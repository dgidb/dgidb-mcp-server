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
from typing import Optional, List

from DGIdb_MCP_server import select_nodes

mcp = FastMCP("CIViC+DGIdb Query App")

HEADERS = {"Content-Type": "application/json"}

# with open("../src/data/disease_name_map.json", "r") as f:
#     d_disease_map = json.load(f)

# with open("../src/data/therapy_name_map.json", "r") as f:
#     d_therapy_map = json.load(f)

# with open("../src/data/molecular_profile_map.json", "r") as f:
#     MP_map = json.load(f)

with open("data/VICC_gene_alias_map.json", "r", encoding='utf-8') as f:
    gene_map = json.load(f)

def load_drug_list_from_txt(path: str):
    with open(path, encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]
    
dgidb_drug_list = load_drug_list_from_txt('data/primary_dgidb_drug_names.txt')

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

def normalize_entity(name, lookup, threshold = 0.7):
    """
    Map a free‑form name to its primary alias via fuzzy matching (Dice/bigrams).

    :param name:    Input string (or None)
    :param lookup:  Dict mapping primary → list of aliases
    :param threshold: Minimum similarity (0–1) to accept
    :return:        The matched primary string, or None if no good match
    """
    if not name:
        return None

    q_norm = normalize_str(name)

    if name.lower() in lookup.keys():
        return name

    # Build normalized‑alias → primary map
    alias_to_primary = {}
    for primary, aliases in lookup.items():
        # also include the primary itself as an alias
        alias_to_primary[normalize_str(primary)] = primary
        for alias in aliases:
            alias_to_primary[normalize_str(alias)] = primary

    # Find best candidate by Dice score
    best_target, best_score = None, -1.0
    for norm_alias in alias_to_primary:
        score = dice_coefficient(q_norm, norm_alias)
        if score > best_score:
            best_score, best_target = score, norm_alias

    if best_score >= threshold:
        return alias_to_primary[best_target]
    return None

def normalize_dgidb_drug(name, lookup=dgidb_drug_list, threshold = 0.7):
    """
    Map a free-form drug name to a canonical DGIdb drug name using a list lookup
    and Sørensen–Dice (bigram) similarity.

    :param name:       Input drug string (or None)
    :param lookup:     List of canonical drug names (e.g., from DGIdb)
    :param threshold:  Minimum similarity (0–1) to accept
    :return:           Best-matching canonical name, or None if no good match
    """
    if not name or not lookup:
        return None

    q_norm = normalize_str(name)

    # Fast path: exact normalized match
    for cand in lookup:
        if cand and normalize_str(cand) == q_norm:
            return cand

    # Otherwise, choose the highest Dice coefficient
    best_name, best_score = None, -1.0
    for cand in lookup:
        if not cand:
            continue
        score = dice_coefficient(q_norm, normalize_str(cand))
        if score > best_score:
            best_name, best_score = cand, score

    return best_name if best_score >= threshold else None


def run_civic_query(query, variables: dict) -> dict:
    print("Running GraphQL Query with variables:", variables, file=sys.stderr)
    r = requests.post(
        "https://civicdb.org/api/graphql",
        json={"query": query, "variables": variables},
        headers=HEADERS,
        timeout=30,
    )
    print("GraphQL Response:", r.text, file=sys.stderr)
    r.raise_for_status()
    return r.json()

def civic_rating_rank(rating):
    """
    Convert a CIViC evidenceRating into a sortable rank.
    Higher numeric rating = better.
    None = worst.
    Returned rank is negated so higher ratings sort first.
    """
    if rating is None:
        return float('inf')  # push None to the end
    
    try:
        val = int(rating)
    except (TypeError, ValueError):
        return float('inf')
    
    # return negative so that higher rating sorts earlier
    return -val

def sort_CIViC_entries(evidence_units):
    # evidenceLevel ordering: A > B > C > D > E
    level_order = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4}

    def sort_key(item):
        # primary = evidenceLevel (A first)
        lvl = item.get('evidenceLevel')
        lvl_rank = level_order.get(lvl, 99)

        # secondary = evidenceRating (higher first)
        rating = item.get('evidenceRating')
        rating_rank = civic_rating_rank(rating)

        return (lvl_rank, rating_rank)

    evidence_units.sort(key=sort_key)
    return evidence_units


@mcp.tool(title="CIViC evidence items for a gene variant (optional) + disease (optional) + therapy (optional). Contains clinical significance information for specific publications.")
def get_variant_evidence(molecularProfileName=None, diseaseName=None, therapyName=None, significance=None):
    """Return CIViC JSON for the requested disease / molecular profile."""

    print(f"Tool called with variant: {molecularProfileName}, disease: {diseaseName}", file=sys.stderr, flush=True)

    if all(arg is None for arg in [molecularProfileName, diseaseName, therapyName]):
        return 'One of the following must be specified: molecularProfileName, diseaseName, therapyName'

    #molecularProfileName = normalize_entity(molecularProfileName, MP_map)
    #diseaseName = normalize_entity(diseaseName, d_disease_map)

    if therapyName:
        therapyName = therapyName.split(', ')[0]

    #therapyName = normalize_entity(therapyName, d_therapy_map)

    #print(f"Normalized name: {molecularProfileName}", file=sys.stderr, flush=True)
    
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

    variables = {
        k: v
        for k, v in {
            "diseaseName": diseaseName,
            "molecularProfileName": molecularProfileName,
            "therapyName": therapyName,
            "significance": significance
        }.items()
        if v is not None and v.lower() != 'none' and v.lower() != 'null'           # keep only non-None values
    }

    resp = run_civic_query(evidence_query, variables)

    if 'data' not in resp: #query error
        return resp
    
    evidence_items = resp['data']['evidenceItems']['nodes']

    evidence_items = sort_CIViC_entries(evidence_items)

    for item in evidence_items:
        print(item['evidenceLevel'], item['evidenceRating'])

    for item in evidence_items:
        eid = item['id']
        item['url'] = f'https://identifiers.org/civic.eid:{eid}'
        del item['id']

    return_object = {}

    return_object['Field Descriptions'] = (
        "evidenceType: Category describing the type of clinical or biological evidence (e.g., predictive, diagnostic).\n"
        "evidenceDirection: Indicates whether the evidence supports or refutes the association.\n"
        "significance: The clinical relevance of the evidence.\n"
        "description: Detailed summary of the evidence from CIViC curators.\n"
        "evidenceLevel: Describes the robustness of the study type. A - Validated association, B - Clinical evidence, C - Case study, D - Preclinical evidence, and E - Inferential association\n"
        "evidenceRating: Quality score assigned to the evidence by curators (scored 1-5).\n"
        "url: Direct link to the CIViC record for this evidence item.\n"
        "When returning information to users you MUST cite URLs used for specific information."
    )

    return_object['API Results'] = evidence_items

    return return_object


@mcp.tool(title="CIViC assertions for a gene variant (optional) + disease (optional) + therapy (optional). Contains clinical significance information across multiple publications.")
def get_variant_assertions(molecularProfileName=None, diseaseName=None, therapyName=None, significance=None):
    """Return CIViC JSON for the requested disease / molecular profile."""

    print(f"Tool called with variant: {molecularProfileName}, disease: {diseaseName}", file=sys.stderr, flush=True)

    #molecularProfileName = normalize_entity(molecularProfileName, MP_map)
    #diseaseName = normalize_entity(diseaseName, d_disease_map)

    if therapyName:
        therapyName = therapyName.split(', ')[0]

    #therapyName = normalize_entity(therapyName, d_therapy_map)

    #print(f"Normalized name: {molecularProfileName}", file=sys.stderr, flush=True)
    
    assertions_query = """
    query assertions($molecularProfileName: String!, $diseaseName: String, $therapyName: String, $significance: EvidenceSignificance) {
        assertions(molecularProfileName: $molecularProfileName, diseaseName: $diseaseName, therapyName: $therapyName, significance: $significance) {
            nodes { 
                status 
                assertionType
                assertionDirection 
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
                summary
                id 
            }
        }
    }"""

    variables = {
        k: v
        for k, v in {
            "diseaseName": diseaseName,
            "molecularProfileName": molecularProfileName,
            "therapyName": therapyName,
            "significance": significance
        }.items()
        if v is not None and v.lower() != 'none' and v.lower() != 'null'           # keep only non-None values
    }

    resp = run_civic_query(assertions_query, variables)

    if 'data' not in resp: #query error
        return resp
    
    assertions = resp['data']['assertions']['nodes']

    #assertions = sort_CIViC_entries(assertions)

    #modify to get a clickable link to CIViC
    for item in assertions:
        aid = item['id']
        item['url'] = f'https://identifiers.org/civic.aid:{aid}'
        del item['id']

    return_object = {}

    return_object['Field Descriptions'] = (
        "assertionType: Category describing the type of clinical or biological evidence (e.g., predictive, diagnostic).\n"
        "assertionDirection: Indicates whether the evidence supports or refutes the association.\n"
        "significance: The clinical relevance of the evidence.\n"
        "summary: Detailed summary of the evidence from CIViC curators.\n"
        "url: Direct link to the CIViC record for this evidence item.\n"
        "When returning information to users you MUST cite URLs used for specific information."
    )

    return_object['API Results'] = assertions

    return return_object

def run_dgidb_query(query, variables: dict) -> dict:
    print("Running GraphQL Query with variables:", variables, file=sys.stderr)
    r = requests.post(
        "https://dgidb.org/api/graphql",
        json={"query": query, "variables": variables},
        headers=HEADERS,
        timeout=30,
    )
    #print("GraphQL Response:", r.text, file=sys.stderr)
    r.raise_for_status()
    return r.json()

def parse_list(s: str):
    # split on comma or whitespace (one or more), strip empties
    return [d.strip() for d in re.split(r'[,\s]+', s.strip()) if d.strip()]

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


@mcp.tool(title="Gets the drugs that interact with a list of genes. Input like 'gene1, gene2'.")
def get_drug_interactions_for_gene_list(gene_names):
    normalized_gene_names = []
    gene_names = parse_list(gene_names)

    for gene in gene_names:
        gene = normalize_entity(gene, gene_map)
        if gene:
            normalized_gene_names.append(gene.upper())

    if not normalized_gene_names:
        return None

    print('norm', normalized_gene_names)
    resp = run_dgidb_query(gene_query, {'names': normalized_gene_names})

    nodes = resp['data']['genes']['nodes']

    if not nodes:
        return resp
    
    return select_nodes(nodes, normalized_gene_names, N=20, min_N=20)


@mcp.tool(title="Gets drug info including approval, if used in immunotherapy, and other drug attributes for a list of drugs. Input like 'drug1, drug2'.")
def get_drug_info(drug_names):
    normalized_drug_names = []
    drug_names = parse_list(drug_names)
    print(drug_names)

    for drug in drug_names:
        drug = normalize_dgidb_drug(drug)
        normalized_drug_names.append(drug)

    print(normalized_drug_names)

    resp = run_dgidb_query(drug_info_query, {'names': normalized_drug_names})

    nodes = resp['data']['drugs']['nodes']
    if not nodes:
        return resp
    #combine all the nodes that represent the same drugs
    #only return the info for the exact string match if possible

    selected_nodes = []
    for drug_norm in normalized_drug_names:
        #first iterate over all of them to find an exact match
        found = False
        for node in nodes:
            drug_name = node['name']
            print(drug_norm.upper(), drug_name.upper(), drug_norm.upper() == drug_name.upper())
            if drug_norm.upper() == drug_name.upper():
                selected_nodes.append(node)
                found = True
                break
        
        #if there is no exact match for any return all of them
        if not found:
            print('not found', drug_norm)
            return nodes
        
    return selected_nodes

    
@mcp.tool(title="Gets the genes that interact for a list of drugs. Input like 'drug1, drug2'.")
def get_gene_interactions_for_drug_list(drug_names):

    normalized_drug_names = []
    drug_names = parse_list(drug_names)
    print(drug_names)

    for drug in drug_names:
        drug = normalize_dgidb_drug(drug)
        normalized_drug_names.append(drug)

    print(normalized_drug_names)

    resp = run_dgidb_query(drug_query, {'names': normalized_drug_names})

    nodes = resp['data']['drugs']['nodes']

    if not nodes:
        return resp
    
    return select_nodes(nodes, normalized_drug_names, N=20, min_N=20)

@mcp.tool(title="Gets gene category info from DGIdb for a list of genes. Input like 'gene1, gene2'.")
def get_gene_categories(gene_names):
    normalized_gene_names = []
    gene_names = parse_list(gene_names)
    print(gene_names)

    for gene in gene_names:
        gene = normalize_entity(gene, gene_map)
        normalized_gene_names.append(gene)

    print(normalized_gene_names)

    resp = run_dgidb_query(gene_category_query, {'names': normalized_gene_names})
    #print(resp)

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
    #print(get_gene_categories('ABCC10,ALCAM'))
    print(get_drug_interactions_for_gene_list('BTK'))
    #get_variant_evidence('BRAF V600E')
    #print(get_variant_evidence(diseaseName='Lung Non-small Cell Carcinoma', therapyName='Paclitaxel', significance='RESISTANCE'))
    #mcp.run(transport="stdio")       # recommended pattern in the quick-start 