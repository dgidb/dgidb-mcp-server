import pandas as pd
import sys, requests
import re, csv, json, math, os

from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from concurrent.futures import ThreadPoolExecutor, as_completed

from host_dgidb_civic_MCP import *

import logging

# Silence urllib3 debug logs
logging.getLogger("urllib3").setLevel(logging.WARNING)
logging.getLogger("urllib3.connectionpool").setLevel(logging.WARNING)

# Set True to see build/normalization diagnostics; False prints only the
# metrics reported in the paper.
VERBOSE = False
def vlog(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)

evidence_query = """
query evidenceItems($molecularProfileName: String, $diseaseName: String, $therapyName: String, $significance: EvidenceSignificance) {
evidenceItems( molecularProfileName: $molecularProfileName, diseaseName: $diseaseName, therapyName: $therapyName, significance: $significance) {
    nodes { 
        status 
        evidenceType
        evidenceDirection 
        significance
        molecularProfile{
            variants{
                name
                feature{
                    name
                }
            }
        }
        disease{
            displayName
        }
        therapies{
            name
        }
        variantOrigin 
        description
        evidenceLevel
        evidenceRating
        id 
    }
}
}
"""

labels_filename = 'data/per_gene_ranked_joint_task_drug_lists_dataset.csv'

if not os.path.exists(labels_filename):
    print('First use: python create_joint_task_dataset.py')
    sys.exit()

labels_df = pd.read_csv(labels_filename)

# -------------------------------------------------------
# ----------------- Normalization -----------------------
# -------------------------------------------------------

GENE_BASE = "https://normalize.cancervariants.org/gene/normalize"
THERAPY_BASE = "https://normalize.cancervariants.org/therapy/normalize"

def _make_session(timeout=20, pool_maxsize=50):
    s = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=0.3,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET"]),
    )
    adapter = HTTPAdapter(
        max_retries=retries,
        pool_connections=pool_maxsize,
        pool_maxsize=pool_maxsize,
    )
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers.update({"User-Agent": "bulk-normalizer/1.0"})
    s.request_timeout = timeout
    return s

def _canonical_gene_id(payload):
    g = payload.get("gene")
    if not isinstance(g, dict):
        return None
    pc = g.get("primaryCode")
    if isinstance(pc, str) and pc.strip():
        return pc.strip().lower()
    mappings = g.get("mappings")
    if isinstance(mappings, list) and mappings:
        coding = (mappings[0] or {}).get("coding") or {}
        cid = coding.get("id")
        if isinstance(cid, str) and cid.strip():
            return cid.strip().lower()
    gid = g.get("id")
    if isinstance(gid, str) and gid.strip():
        return gid.strip().lower()
    return None

def _canonical_therapy_id(payload):
    t = payload.get("therapy")
    if not isinstance(t, dict):
        return None
    pc = t.get("primaryCode")
    if isinstance(pc, str) and pc.strip():
        return pc.strip().lower()
    mappings = t.get("mappings")
    if isinstance(mappings, list) and mappings:
        coding = (mappings[0] or {}).get("coding") or {}
        cid = coding.get("id")
        if isinstance(cid, str) and cid.strip():
            return cid.strip().lower()
    tid = t.get("id")
    if isinstance(tid, str) and tid.strip():
        return tid.strip().lower()
    return None

def parse_list_field(field):
    if field is None or (isinstance(field, float) and math.isnan(field)):
        return []
    s = str(field).strip()
    if not s:
        return []
    return [tok.strip() for tok in s.split(",") if tok.strip()]

def to_canonical_set_fast(items, lookup):
    out = set()
    for x in items:
        if not x:
            continue
        cid = lookup.get(x, None)
        if cid:
            out.add(cid)
    return out

def canonicalize_list(items, lookup):
    """
    Convert a list of raw terms to a list of canonical IDs,
    preserving order and dropping duplicates while keeping
    the first occurrence.
    """
    out = []
    seen = set()
    for x in items:
        if not x:
            continue
        cid = lookup.get(x, None)
        if cid and cid not in seen:
            seen.add(cid)
            out.append(cid)
    return out

def _norm_one(session, base_url, q):
    try:
        r = session.get(base_url, params={"q": q}, timeout=session.request_timeout)
        r.raise_for_status()
        data = r.json()
        if base_url == GENE_BASE:
            return _canonical_gene_id(data)
        else:
            return _canonical_therapy_id(data)
    except Exception as e:
        print(f"Normalization error for {q} @ {base_url}: {e}")
        return None

def bulk_normalize(terms, base_url, max_workers=16):
    uniques = sorted({t for t in terms if isinstance(t, str) and t.strip()})
    session = _make_session()
    out = {}

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = {ex.submit(_norm_one, session, base_url, q): q for q in uniques}
        for fut in as_completed(futs):
            q = futs[fut]
            out[q] = fut.result()

    return out

# ----------------- Relevance scoring for genes -------------------

def evidence_level_score(level):
    """
    Map CIViC evidenceLevel (A strongest -> E weakest) to a numeric score.
    Unknown / null -> 0.
    """
    if isinstance(level, float) and math.isnan(level):
        return 0
    if level is None:
        return 0
    s = str(level).strip().upper()
    mapping = {"A": 5, "B": 4, "C": 3, "D": 2, "E": 1}
    return mapping.get(s, 0)

def evidence_rating_score(rating):
    """
    Map CIViC evidenceRating (5 strongest -> 1 weakest, null lowest) to numeric.
    """
    if isinstance(rating, float) and math.isnan(rating):
        return 0
    if rating is None or rating == "":
        return 0
    try:
        r = int(rating)
    except Exception:
        return 0
    if r < 1 or r > 5:
        return 0
    return r

def gene_relevance(level, rating):
    """
    Combine evidenceLevel and evidenceRating into a single relevance score.
    Higher score = stronger evidence.
    Level is weighted more heavily than rating to enforce A > B > C > D > E
    regardless of rating.
    """
    lvl = evidence_level_score(level)
    rat = evidence_rating_score(rating)
    return lvl * 10 + rat  # e.g., B5 > B4, and all Bx > Cx

# ----------------- NDCG utilities -------------------

def ndcg_at_k_binary(pred_list, gold_set, k=20):
    """
    NDCG@k with binary relevance (1 if item in gold_set else 0),
    used for drug lists.

    Special rule for no-gold-drug cases:
    - If len(gold_set) == 0:
        * return 1.0 if pred_list is empty
        * return 0.0 otherwise
    """
    if len(gold_set) == 0:
        return 1.0 if len(pred_list) == 0 else 0.0

    dcg = 0.0
    for i, item in enumerate(pred_list[:k], start=1):
        rel = 1.0 if item in gold_set else 0.0
        if rel > 0.0:
            dcg += rel / math.log2(i + 1)

    ideal_len = min(k, len(gold_set))
    ideal_dcg = 0.0
    for i in range(1, ideal_len + 1):
        ideal_dcg += 1.0 / math.log2(i + 1)

    if ideal_dcg == 0.0:
        return 0.0
    return dcg / ideal_dcg

def ndcg_from_relevance(pred_list, gold_rels, k=None):
    """
    NDCG for gene ranking with graded relevance.

    gold_rels: dict[item] -> relevance (>= 0)
    pred_list: ranked list of items
    k: cutoff; if None, use len(gold_rels)

    Special rule for no-gold-gene cases:
    - If len(gold_rels) == 0:
        * return 1.0 if pred_list is empty
        * return 0.0 otherwise
    """
    if not gold_rels:
        return 1.0 if len(pred_list) == 0 else 0.0

    if k is None:
        k = len(gold_rels)

    dcg = 0.0
    for i, item in enumerate(pred_list[:k], start=1):
        rel = gold_rels.get(item, 0.0)
        if rel > 0.0:
            dcg += rel / math.log2(i + 1)

    # Ideal DCG: sort gold items by relevance descending
    sorted_rels = sorted(gold_rels.values(), reverse=True)
    ideal_rels = sorted_rels[:k]
    ideal_dcg = 0.0
    for i, rel in enumerate(ideal_rels, start=1):
        if rel > 0.0:
            ideal_dcg += rel / math.log2(i + 1)

    if ideal_dcg == 0.0:
        return 0.0
    return dcg / ideal_dcg

# -------------------------------------------------------
# ----------------- BUILD GOLD STRUCTURES ---------------
# -------------------------------------------------------

key_cols = [
    "molecularProfile_name",
    "disease_name",
    "therapies",
    "evidenceType",
    "significance",
]

# -----------------------------------------------
# Precompute normalization for GOLD ONLY (once)
# -----------------------------------------------

gold_gene_terms = labels_df["label_CIViC_genes"].dropna().astype(str).tolist()

gold_drug_terms = []
for s in labels_df["label_DGIdb_drugs"]:
    gold_drug_terms.extend(parse_list_field(s))

unique_gold_genes = set(gold_gene_terms)
unique_gold_drugs = set(gold_drug_terms)

vlog("Normalizing GOLD terms only once...")
vlog("Unique gold genes:", len(unique_gold_genes))
vlog("Unique gold drugs:", len(unique_gold_drugs))

gold_gene_map = bulk_normalize(unique_gold_genes, GENE_BASE, max_workers=12)
gold_drug_map = bulk_normalize(unique_gold_drugs, THERAPY_BASE, max_workers=12)

# GOLD:
# - gold_gene_rels_by_key: key -> {gene_canon -> relevance}
# - gold_gene_pairs: set of (key..., gene_canon)
# - gold_drug_lists: (key..., gene_canon) -> ordered list of canonical drugs
# - gold_drug_sets: (key..., gene_canon) -> set of canonical drugs
# - gold_pairs: set of (key..., gene_canon, drug_canon)

gold_gene_rels_by_key = {}
gold_gene_pairs = set()

gold_drug_lists = {}   # (key..., gene) -> ordered list
gold_drug_sets = {}    # (key..., gene) -> set
gold_pairs = set()

for _, row in labels_df.iterrows():
    key = tuple(row[c] for c in key_cols)

    gene_raw = str(row["label_CIViC_genes"]).strip()
    gene_canon = gold_gene_map.get(gene_raw, None)
    if not gene_canon:
        continue

    # compute relevance for this gene from evidence_level/rating
    level = row.get("evidence_level", None)
    rating = row.get("evidence_rating", None)
    rel = gene_relevance(level, rating)

    # track strongest evidence per (key, gene)
    gene_rels = gold_gene_rels_by_key.setdefault(key, {})
    prev_rel = gene_rels.get(gene_canon, -1)
    if rel > prev_rel:
        gene_rels[gene_canon] = rel

    gold_gene_pairs.add(key + (gene_canon,))

    # drugs
    drugs_raw = parse_list_field(row["label_DGIdb_drugs"])
    drug_canon_list = canonicalize_list(drugs_raw, gold_drug_map)
    drug_canon_set = set(drug_canon_list)

    key_gene = key + (gene_canon,)
    gold_drug_lists[key_gene] = drug_canon_list
    gold_drug_sets[key_gene] = drug_canon_set

    for drug_canon in drug_canon_set:
        gold_pairs.add(key_gene + (drug_canon,))

vlog("Total gold canonical gene-drug pairs:", len(gold_pairs))
vlog("Total gold (key, gene) pairs:", len(gold_gene_pairs))
vlog("Total gold contexts (key) for genes:", len(gold_gene_rels_by_key))
vlog("--------------------------------------------------\n")


# ======================================================
#   TIE-AWARE GOLD (boundary ties at the top-20 cutoff)
# ======================================================
# At the top-20 cutoff, drugs that share the boundary interaction score within
# the same FDA tier are interchangeable: the model should get credit for listing
# ANY of them for the remaining slot(s), must not be penalised for listing a
# different tied drug, and must not be required to list all of them. We build,
# per gene, from live DGIdb interaction scores:
#   core      : canonical drugs ranked ABOVE the boundary tie group (required)
#   tie       : canonical drugs sharing the boundary score (interchangeable)
#   slots     : how many tie members fit inside the top-K (= K - |core|)
#   required  : the truncated gold size G (recall denominator)
#   relevant  : core | tie (every drug that counts as correct)
# Genes with no boundary tie reduce exactly to the strict top-K gold, so only
# the boundary-tie genes are affected.
#
# This re-derives the gold ranking from current DGIdb. It needs network access
# to DGIDB_GRAPHQL_URL and the VICC therapy normaliser, and regenerates the gold
# deterministically (also correcting the small DGIdb drift / non-strict-20 cap
# seen earlier). No LLM output is regenerated.

import json

DGIDB_GRAPHQL_URL = "https://dgidb.org/api/graphql"
TOP_K = 20
SCORE_ROUND = 6
DGIDB_CACHE = "data/dgidb_interactions_cache.json"

GENE_INTERACTIONS_QUERY = """
query genes($names: [String!]!) {
  genes(names: $names) {
    nodes {
      name
      interactions { drug { name approved } interactionScore }
    }
  }
}
"""


def fetch_dgidb_interactions(gene_symbols, batch=25):
    """raw gene symbol -> [(drug_upper, approved_bool, score_or_None), ...].
    ADAPTER: replace the network call with your host_dgidb_civic_MCP fetcher to
    match the exact data your gold pipeline used. Results are cached to disk."""
    cache = {}
    if os.path.exists(DGIDB_CACHE):
        try:
            with open(DGIDB_CACHE) as f:
                cache = json.load(f)
        except Exception:
            cache = {}

    need = sorted({g for g in gene_symbols if g and g not in cache})
    if need:
        sess = _make_session()
        for i in range(0, len(need), batch):
            chunk = need[i:i + batch]
            r = sess.post(
                DGIDB_GRAPHQL_URL,
                json={"query": GENE_INTERACTIONS_QUERY, "variables": {"names": chunk}},
                headers={"Content-Type": "application/json"},
                timeout=60,
            )
            r.raise_for_status()
            nodes = (((r.json() or {}).get("data") or {}).get("genes") or {}).get("nodes") or []
            for node in nodes:
                gname = (node.get("name") or "").strip()
                inters = []
                for it in node.get("interactions") or []:
                    d = it.get("drug") or {}
                    dname = (d.get("name") or "").strip().upper()
                    if not dname:
                        continue
                    sc = it.get("interactionScore")
                    sc = round(float(sc), SCORE_ROUND) if sc is not None else None
                    inters.append([dname, bool(d.get("approved")), sc])
                cache[gname] = inters
            for g in chunk:  # genes DGIdb returned nothing for
                cache.setdefault(g, [])
        try:
            with open(DGIDB_CACHE, "w") as f:
                json.dump(cache, f)
        except Exception:
            pass

    return {g: [tuple(x) for x in cache.get(g, [])] for g in gene_symbols}


def _rank_raw(inters):
    """Two-tier ranking: FDA-approved first, then interaction score desc."""
    def key(t):
        _, ap, sc = t
        return (0 if ap else 1, -(sc if sc is not None else -math.inf))
    return sorted(inters, key=key)


def _core_tie(raw_ranked, top_k=TOP_K):
    """(core_raw_names, tie_raw_names). tie is non-empty ONLY when a score tie
    straddles the top_k cutoff (the only case the allowance changes anything)."""
    n = len(raw_ranked)
    if n <= top_k:
        return [d[0] for d in raw_ranked], []
    b_ap, b_sc = raw_ranked[top_k - 1][1], raw_ranked[top_k - 1][2]
    if raw_ranked[top_k][1] == b_ap and raw_ranked[top_k][2] == b_sc:
        first_idx = next(i for i, (_, ap, sc) in enumerate(raw_ranked)
                         if ap == b_ap and sc == b_sc)
        core = [d[0] for d in raw_ranked[:first_idx]]
        tie = [d[0] for d in raw_ranked if d[1] == b_ap and d[2] == b_sc]
        return core, tie
    return [d[0] for d in raw_ranked[:top_k]], []


# strict canonical gold set per gene_canon (from the frozen gold already built)
gold_strict_set_by_gene = {}
for key_gene, dset in gold_drug_sets.items():
    gold_strict_set_by_gene.setdefault(key_gene[-1], dset)  # identical across keys

# gene_canon -> a representative raw gene symbol to query DGIdb
raw_symbol_by_gene_canon = {}
for raw_gene, gcanon in gold_gene_map.items():
    if gcanon and gcanon not in raw_symbol_by_gene_canon:
        raw_symbol_by_gene_canon[gcanon] = raw_gene

vlog("Fetching DGIdb interaction scores for tie-aware gold...")
raw_inter = fetch_dgidb_interactions(sorted(set(raw_symbol_by_gene_canon.values())))

# canonicalise only the names we need (top-K plus any boundary tie group)
names_to_norm = set()
ranked_by_gene_canon = {}
for gcanon, raw_gene in raw_symbol_by_gene_canon.items():
    ranked = _rank_raw(raw_inter.get(raw_gene, []))
    ranked_by_gene_canon[gcanon] = ranked
    core_raw, tie_raw = _core_tie(ranked)
    names_to_norm.update(core_raw)
    names_to_norm.update(tie_raw)

drug_norm = bulk_normalize(names_to_norm, THERAPY_BASE, max_workers=12)


def _canon_set(names):
    return {drug_norm[n] for n in names if drug_norm.get(n)}


# Build two gold maps from the SAME regenerated ranking:
#   strict_gold_fetched : truncate at 20, no tie allowance
#   tieaware_gold       : boundary-tie allowance
# (tie-aware - strict_fetched) therefore isolates the tie effect from any drift.
strict_gold_fetched = {}
tieaware_gold = {}
n_boundary = 0
n_no_inter = 0
for gcanon, ranked in ranked_by_gene_canon.items():
    if not ranked:
        n_no_inter += 1
    core_raw, tie_raw = _core_tie(ranked)
    strict_relevant = _canon_set([d[0] for d in ranked[:TOP_K]])
    G = len(strict_relevant)

    strict_gold_fetched[gcanon] = dict(
        core=set(strict_relevant), tie=set(), slots=0,
        required=G, relevant=set(strict_relevant),
    )

    core_c = _canon_set(core_raw)
    tie_c = _canon_set(tie_raw) - core_c
    slots = max(0, min(G - len(core_c), len(tie_c)))
    tieaware_gold[gcanon] = dict(
        core=core_c, tie=tie_c, slots=slots,
        required=G, relevant=core_c | tie_c,
    )
    if tie_c:
        n_boundary += 1

vlog(f"Tie-aware gold built for {len(tieaware_gold)} genes; "
      f"{n_boundary} have a boundary tie group; "
      f"{n_no_inter} returned no interactions.")
vlog("--------------------------------------------------\n")


def score_drug_lists(items, gold_map, top_k=TOP_K):
    """items: iterable of (gene_canon, predicted_ordered_drug_list).
    gold_map: gene_canon -> {core, tie, slots, required, relevant}.

    Capacity-aware micro precision / recall / F1 + mean NDCG@k:
      - precision numerator = predicted drugs that are in `relevant`
      - recall   numerator  = |P n core| + min(|P n tie|, slots)   (capped)
      - recall   denominator= required (= strict gold size G)
    So any boundary-tied drug counts toward the remaining slot(s), listing extra
    tied drugs is not penalised, and unambiguous (core) drugs remain required."""
    sum_corr = sum_pred = sum_tp = sum_req = 0
    ndcgs = []
    for gene_canon, pred_list in items:
        g = gold_map.get(gene_canon)
        if g is None:
            core = tie = relevant = set()
            slots = required = 0
        else:
            core, tie = g["core"], g["tie"]
            slots, required, relevant = g["slots"], g["required"], g["relevant"]
        P = set(pred_list)
        sum_corr += len(P & relevant)
        sum_pred += len(P)
        sum_tp += len(P & core) + min(len(P & tie), slots)
        sum_req += required
        ndcgs.append(ndcg_at_k_binary(pred_list, relevant, k=top_k))
    precision = sum_corr / sum_pred if sum_pred else 0.0
    recall = sum_tp / sum_req if sum_req else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    ndcg = sum(ndcgs) / len(ndcgs) if ndcgs else 0.0
    return dict(contexts=len(ndcgs), precision=precision, recall=recall, f1=f1, ndcg=ndcg)


# ======================================================
#   Evaluate MCP and No-MCP **without redoing gold norms**
# ======================================================

for LLM_filename in [
    "data/per_gene_GPT_MCP_joint_rank.csv",
    "data/per_gene_GPT_no_MCP_joint_rank.csv"
]:

    print("\n=========== Evaluating:", LLM_filename, "============\n")

    LLM_df = pd.read_csv(LLM_filename)

    # Collect predicted terms (for normalization)
    pred_gene_terms = LLM_df["LLM_CIViC_genes"].dropna().astype(str).tolist()

    pred_drug_terms = []
    for s in LLM_df["LLM_DGIdb_drugs"]:
        pred_drug_terms.extend(parse_list_field(s))

    # Combine gold + pred for THIS evaluation run
    unique_pred_genes = unique_gold_genes | set(pred_gene_terms)
    unique_pred_drugs = unique_gold_drugs | set(pred_drug_terms)

    vlog("Unique genes to normalize for this model:", len(unique_pred_genes))
    vlog("Unique drugs to normalize for this model:", len(unique_pred_drugs))

    # Normalize only prediction-side deltas
    pred_gene_map = bulk_normalize(unique_pred_genes, GENE_BASE, max_workers=12)
    pred_drug_map = bulk_normalize(unique_pred_drugs, THERAPY_BASE, max_workers=12)

    # Build predicted structures:
    # - pred_gene_pairs: set of (key..., gene)
    # - pred_gene_lists_by_key: key -> ordered list of canonical predicted genes
    # - pred_pairs: set of (key..., gene, drug)
    # - pred_drug_lists: (key..., gene) -> ordered list of canonical predicted drugs
    # - pred_drug_sets: (key..., gene) -> set of canonical predicted drugs

    pred_gene_pairs = set()
    pred_gene_lists_by_key = {}
    pred_pairs = set()
    pred_drug_lists = {}
    pred_drug_sets = {}

    # iterate rows in file order so that ordering of genes per key matches LLM ranking
    for _, row in LLM_df.iterrows():
        key = tuple(row[c] for c in key_cols)

        # gene prediction
        val_gene = row.get("LLM_CIViC_genes", None)
        if isinstance(val_gene, float) and math.isnan(val_gene):
            gene_raw = ""
        else:
            gene_raw = str(val_gene).strip()

        gene_canon = pred_gene_map.get(gene_raw, None)
        if not gene_canon:
            # no usable gene prediction; skip this row
            continue

        # maintain gene set and ranked list per key
        pred_gene_pairs.add(key + (gene_canon,))
        gene_list = pred_gene_lists_by_key.setdefault(key, [])
        if gene_canon not in gene_list:
            gene_list.append(gene_canon)

        # drug prediction for this (key, gene)
        val_drugs = row.get("LLM_DGIdb_drugs", None)
        drugs_raw = parse_list_field(val_drugs)
        drug_canon_list = canonicalize_list(drugs_raw, pred_drug_map)
        drug_canon_set = set(drug_canon_list)

        key_gene = key + (gene_canon,)
        pred_drug_lists[key_gene] = drug_canon_list
        pred_drug_sets[key_gene] = drug_canon_set

        for drug_canon in drug_canon_set:
            pred_pairs.add(key_gene + (drug_canon,))

    # =====================================================================
    # STEP 1 — CIViC gene identification (P / R / F1 + NDCG)
    # =====================================================================
    gene_tp_count = len(gold_gene_pairs & pred_gene_pairs)
    gene_fp_count = len(pred_gene_pairs - gold_gene_pairs)
    gene_fn_count = len(gold_gene_pairs - pred_gene_pairs)

    gene_precision = gene_tp_count / (gene_tp_count + gene_fp_count) if (gene_tp_count + gene_fp_count) else 0.0
    gene_recall = gene_tp_count / (gene_tp_count + gene_fn_count) if (gene_tp_count + gene_fn_count) else 0.0
    gene_f1 = (2 * gene_precision * gene_recall / (gene_precision + gene_recall)) if (gene_precision + gene_recall) else 0.0

    gene_ndcg_values = []
    for key in set(gold_gene_rels_by_key) | set(pred_gene_lists_by_key):
        gene_ndcg_values.append(
            ndcg_from_relevance(pred_gene_lists_by_key.get(key, []),
                                gold_gene_rels_by_key.get(key, {}), k=None))
    mean_gene_ndcg = sum(gene_ndcg_values) / len(gene_ndcg_values) if gene_ndcg_values else 0.0

    # =====================================================================
    # STEP 2 — DGIdb drug retrieval (PRIMARY = unique-gene, tie-aware)
    # =====================================================================
    # Collapse predictions to one entry per unique canonical gene (first
    # occurrence) so recurring genes are not over-weighted, then score with the
    # tie-aware gold (any drug tied at the top-20 boundary counts toward the
    # remaining slot[s]). The strict / per-occurrence variants are computed too,
    # for the robustness table only.
    def _collapse_to_unique_gene(drug_lists, drug_sets):
        out_list, out_set = {}, {}
        for key_gene, lst in drug_lists.items():
            g = key_gene[-1]
            if g not in out_list:
                out_list[g] = lst
                out_set[g] = drug_sets.get(key_gene, set())
        return out_list, out_set

    pred_drug_lists_g, _ = _collapse_to_unique_gene(pred_drug_lists, pred_drug_sets)

    po_keys = set(gold_drug_sets) | set(pred_drug_lists)
    po_items = [(kg[-1], pred_drug_lists.get(kg, [])) for kg in po_keys]
    ug_keys = set(gold_strict_set_by_gene) | set(pred_drug_lists_g)
    ug_items = [(g, pred_drug_lists_g.get(g, [])) for g in ug_keys]

    po_strict = score_drug_lists(po_items, strict_gold_fetched)
    po_tie = score_drug_lists(po_items, tieaware_gold)
    ug_strict = score_drug_lists(ug_items, strict_gold_fetched)
    ug_tie = score_drug_lists(ug_items, tieaware_gold)  # <-- PRIMARY

    # ----------------------------- report --------------------------------
    model = "GPT-5 + MCP" if "no_MCP" not in LLM_filename else "GPT-5 (no MCP)"
    print(f"\n================= {model} =================")

    print("\nStep 1 — CIViC gene identification")
    print(f"  Precision {gene_precision:.2f}   Recall {gene_recall:.2f}   "
          f"F1 {gene_f1:.2f}   NDCG {mean_gene_ndcg:.2f}")

    print("\nStep 2 — DGIdb drug retrieval  (unique-gene, tie-aware)")
    print(f"  Precision {ug_tie['precision']:.2f}   Recall {ug_tie['recall']:.2f}   "
          f"F1 {ug_tie['f1']:.2f}   NDCG {ug_tie['ndcg']:.2f}")

    print("\nSupplementary — drug-step robustness "
          "(gene recurrence x tie handling)")
    print(f"  {'variant':<20}{'P':>7}{'R':>7}{'F1':>7}{'NDCG':>7}")
    for lab, d in [("per-occ strict", po_strict), ("per-occ tie-aware", po_tie),
                   ("unique  strict", ug_strict), ("unique  tie-aware", ug_tie)]:
        print(f"  {lab:<20}{d['precision']:>7.2f}{d['recall']:>7.2f}"
              f"{d['f1']:>7.2f}{d['ndcg']:>7.2f}")
    print()