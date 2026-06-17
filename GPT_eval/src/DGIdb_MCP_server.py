from mcp.server.fastmcp import FastMCP
import requests, json, sys
from difflib import SequenceMatcher
import logging
logging.basicConfig(level=logging.DEBUG)
import re
import unicodedata
from collections import Counter
import functools, sys
print = functools.partial(print, file=sys.stderr, flush=True) 
import os


gene_query = """
query genes($names: [String!]) {
    genes(names: $names) {
        nodes {
            name
            interactions {
                drug {
                    name
                    approved
                }
                interactionScore
                interactionTypes {
                    type
                    directionality
                }
            }
        }
    }
}
"""

drug_query = """
query drugs($names: [String!]) {
    drugs(names: $names) {
        nodes {
            name
            interactions {
                gene {
                    name
                }
                interactionScore
                interactionTypes {
                    type
                    directionality
                }
            }
        }
    }
}
"""

drug_info_query = """
query drugs_info($names: [String!]) {
    drugs(names: $names) {
        nodes {
            name
            conceptId
            approved
            immunotherapy
            antiNeoplastic
            drugAttributes {
                name
                value
                }
            drugApprovalRatings {
                rating
                source {
                    sourceDbName
                    sourceTrustLevel {
                        level
                    }
                }
            }
        }
    }
}
"""

gene_category_query = """
query gene_category($names: [String!]){
    genes(names: $names) {
        nodes {
            name
            longName
            conceptId
            geneCategoriesWithSources {
            name
            sourceNames
            }
        }
    }
}
"""


def filter_and_sort(interactions, N=20): #N should depend on how many genes/drugs were input
    # Sort key:
    # - bool converts to int (False=0, True=1), so negate it (not approved=0, approved=1 → sort approved first)
    # - secondary sort is interactionScore
    sorted_interactions = sorted(
        interactions,
        key=lambda x: (
            not x.get("drug", {}).get("approved", False),   # False if approved → comes first
            -x.get("interactionScore", 0)                   # negative for descending
        )
    )
    return sorted_interactions[:N]

def combine_matching_nodes(nodes, term):
    """combine nodes that contain the substring {term}"""
    term_node = {'interactions': []}

    for node in nodes:
        node_name = node['name']

        if term in node_name:
            for interaction in node['interactions']:
                term_node['interactions'].append(interaction)

    return term_node

def distribute_nodes(d_node_lens, total = 100):
    """
    d_node_lens has something like {'drug1': 25, 'drug2': 75}

    sort the nodes by len in ascending order
    Iteration:
        compute the average size to  len / total
        allow up to this much for the current node
        len-=1
        total-=size(i)

    """

    sorted_d_node_lens = {k: v for k, v in sorted(d_node_lens.items(), key=lambda item: item[1])}
    i = 0 #iteration count
    num_nodes = len(sorted_d_node_lens)

    updated_d_node_lens = {}

    for key in sorted_d_node_lens:
        val = sorted_d_node_lens[key]
        allowed_count = int(total / (num_nodes - i))

        curr_count = min(val, allowed_count)
        total -= curr_count

        updated_d_node_lens[key] = curr_count

        i+=1
    
    return updated_d_node_lens

def select_nodes(nodes, normalized_names, N=100, min_N=40):
    #Count how many interactions are present for each gene to create d_node_lens for distribute_nodes()
    d_node_counts = {}
    d_name_nodes = {} #contains all the interactions that match the gene/drug name (node)

    for name in normalized_names:
        selected_node = combine_matching_nodes(nodes, name)
        d_name_nodes[name] = selected_node
        d_node_counts[name] = len(selected_node['interactions'])

    total = N

    if len(nodes) == 1:
        total = min_N

    d_node_counts = distribute_nodes(d_node_counts, total=total)

    for name in d_node_counts:
        d_name_nodes[name] = filter_and_sort(d_name_nodes[name]['interactions'], N=d_node_counts[name])

    return d_name_nodes
    

mcp = FastMCP("DGIdb Query App")

HEADERS = {"Content-Type": "application/json"}

with open("data/VICC_drug_alias_map.json", "r", encoding='utf-8') as f:
    drug_map = json.load(f)

with open("data/VICC_gene_alias_map.json", "r", encoding='utf-8') as f:
    gene_map = json.load(f)

def normalize_str(s: str) -> str:
    """
    Decompose accents, strip diacritics, lowercase, remove non-alphanumerics,
    collapse whitespace, and trim.
    """
    # 1. Decompose accents (NFKD) and strip combining marks
    s = unicodedata.normalize('NFKD', s)
    s = ''.join(ch for ch in s if not unicodedata.combining(ch))
    # 2. Lowercase
    s = s.lower()
    # 3. Keep only letters, digits, spaces
    s = re.sub(r'[^a-z0-9\s]', ' ', s)
    # 4. Collapse runs of whitespace
    s = re.sub(r'\s+', ' ', s)
    return s.strip()

def dice_coefficient(a: str, b: str) -> float:
    """
    Compute Sørensen–Dice coefficient between two strings based on bigrams.
    Returns 1.0 for exact match on <2‑char strings, 0.0 if either is empty.
    """
    if not a or not b:
        return 0.0
    if len(a) < 2 or len(b) < 2:
        return 1.0 if a == b else 0.0

    # Build bigram multisets
    ba = [a[i:i+2] for i in range(len(a) - 1)]
    bb = [b[i:i+2] for i in range(len(b) - 1)]
    ca, cb = Counter(ba), Counter(bb)

    # Intersection size
    intersection = sum(min(ca[gram], cb[gram]) for gram in ca)
    total = sum(ca.values()) + sum(cb.values())

    return (2.0 * intersection) / total if total > 0 else 0.0


def normalize_entity(name, lookup, threshold=0.7):
    """
    Map a free-form name to its primary alias via fuzzy matching (Dice/bigrams).
    Tie-break: prefer the primary alias over non-primary aliases.
    Also handle alias collisions that normalize to the same string.
    """
    if not name:
        return None

    q_norm = normalize_str(name)

    # Collect candidates: (primary_str, normalized_alias_str, is_primary, alias_order)
    candidates = []
    primaries_norm = {}  # normalized primary -> primary (for exact primary match)
    for primary, aliases in lookup.items():
        p_norm = normalize_str(primary)
        primaries_norm[p_norm] = primary
        candidates.append((primary, p_norm, True, -1))  # primary has alias_order -1
        for idx, alias in enumerate(aliases):
            candidates.append((primary, normalize_str(alias), False, idx))

    # 1) Hard preference: exact match to a PRIMARY alias
    if q_norm in primaries_norm:
        return primaries_norm[q_norm].upper()

    # 2) Score all candidates; prefer primaries on ties
    best = None
    best_key = None  # tuple used for tie-breaking
    EPS = 1e-12

    for primary, a_norm, is_primary, alias_order in candidates:
        score = dice_coefficient(q_norm, a_norm)

        # Sort key (higher is better):
        #   1) score
        #   2) is_primary (True > False)
        #   3) exact string equality (rare tie-breaker if scores are equal)
        #   4) shorter alias length (slight nudge toward tighter match)
        #   5) smaller alias_order (earlier-listed alias wins)
        exact = (a_norm == q_norm)
        key = (
            round(score, 12),
            1 if is_primary else 0,
            1 if exact else 0,
            -len(a_norm),
            -alias_order,
        )

        if best is None or key > best_key:
            best = (primary, score)
            best_key = key

    if best and best[1] >= threshold:
        return best[0].upper()
    return None


def run_query(query, variables: dict) -> dict:
    print("Running GraphQL Query with variables:", variables, file=sys.stderr)
    r = requests.post(
        "https://dgidb.org/api/graphql",
        json={"query": query, "variables": variables},
        headers=HEADERS,
        timeout=30,
    )
    print("GraphQL Response:", r.text, file=sys.stderr)
    r.raise_for_status()
    return r.json()

def parse_list(s: str):
    # split on comma or whitespace (one or more), strip empties
    return [d.strip() for d in re.split(r'[,\s]+', s.strip()) if d.strip()]


@mcp.tool(title="Gets the drugs that interact with a list of genes. Input like 'gene1, gene2'.")
def get_drug_interactions_for_gene_list(gene_names):
    normalized_gene_names = []
    gene_names = parse_list(gene_names)

    for gene in gene_names:
        gene = normalize_entity(gene, gene_map)
        normalized_gene_names.append(gene)

    resp = run_query(gene_query, {'names': normalized_gene_names})

    nodes = resp['data']['genes']['nodes']

    if not nodes:
        return resp
    
    return select_nodes(nodes, normalized_gene_names)


@mcp.tool(title="Gets drug info including approval, if used in immunotherapy, and other drug attributes for a list of drugs. Input like 'drug1, drug2'.")
def get_drug_info(drug_names):
    normalized_drug_names = []
    drug_names = parse_list(drug_names)
    print(drug_names)

    for drug in drug_names:
        #drug = normalize_entity(drug, drug_map)
        normalized_drug_names.append(drug)

    print(normalized_drug_names)

    resp = run_query(drug_info_query, {'names': normalized_drug_names})

    nodes = resp['data']['drugs']['nodes']
    return nodes

    
@mcp.tool(title="Gets the genes that interact for a list of drugs. Input like 'drug1, drug2'.")
def get_gene_interactions_for_drug_list(drug_names):

    normalized_drug_names = []
    drug_names = parse_list(drug_names)
    print(drug_names)

    for drug in drug_names:
        drug = normalize_entity(drug, drug_map)
        normalized_drug_names.append(drug)

    print(normalized_drug_names)

    resp = run_query(drug_query, {'names': normalized_drug_names})

    nodes = resp['data']['drugs']['nodes']

    if not nodes:
        return resp
    
    return select_nodes(nodes, normalized_drug_names)

@mcp.tool(title="Gets gene category info from DGIdb for a list of genes. Input like 'gene1, gene2'.")
def get_gene_categories(gene_names):
    normalized_gene_names = []
    gene_names = parse_list(gene_names)
    print(gene_names)

    for gene in gene_names:
        #gene = normalize_entity(gene, gene_map)
        normalized_gene_names.append(gene)

    print(normalized_gene_names)

    resp = run_query(gene_category_query, {'names': normalized_gene_names})
    print(resp)

    nodes = resp['data']['genes']['nodes']
    if not nodes:
        return resp
    #combine all the nodes that represent the same drugs
    #only return the info for the exact string match if possible

    selected_nodes = []
    for gene_norm in normalized_gene_names:
        #first iterate over all of them to find an exact match
        found = False
        for node in nodes:
            gene_name = node['name']
            print(gene_norm.upper(), gene_name.upper(), gene_norm.upper() == gene_name.upper())
            if gene_norm.upper() == gene_name.upper():
                selected_nodes.append(node)
                found = True
                break
        
        #if there is no exact match for any return all of them
        if not found:
            print('not found', gene_norm)
            return nodes
        
    return selected_nodes


if __name__ == "__main__":           # <-- only runs when *you* call the file
    #get_gene_interactions_for_drug_list('CRISABOROLE')
    print(get_drug_info('FARUDODSTAT,[D-TYR8]CYN 154806'))
    #print(get_gene_categories('CRAT,GRIN2A'))
    #mcp.run(transport="stdio")       # recommended pattern in the quick-start 