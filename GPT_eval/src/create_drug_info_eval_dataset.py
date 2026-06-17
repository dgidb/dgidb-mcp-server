import pandas as pd
from datetime import datetime
from openai import AsyncOpenAI
from agents import Agent, Runner, OpenAIChatCompletionsModel, ModelSettings
from agents.mcp import MCPServerStdio
import csv
import random, re, hashlib
import numpy as np
import time, io
import requests, json, sys
#use csv index to keep track of MCP output files.

drugs_df = pd.read_csv('data/drugs.tsv', sep='\t')
print(len(drugs_df))
#can one drug have multiple answers for approved, immunotherapy, anti_neoplastic
#remove any drugs where this is the case

conflicting_drugs = (
    drugs_df
    .groupby("drug_name")[["approved", "immunotherapy", "anti_neoplastic"]]
    .nunique()
    .reset_index()
)

# Identify drugs with >1 unique value in any of the three columns
conflicting_drugs = conflicting_drugs[
    (conflicting_drugs[["approved", "immunotherapy", "anti_neoplastic"]] > 1).any(axis=1)
]

print(conflicting_drugs)

# Filter them out from the dataframe
drugs_df = drugs_df[~drugs_df["drug_name"].isin(conflicting_drugs["drug_name"])]
drugs_df = drugs_df.drop_duplicates(subset=["drug_name"], keep="first").reset_index(drop=True)

print(len(set(drugs_df['drug_name'])))

strata_counts = (
    drugs_df
    .groupby(["approved", "immunotherapy", "anti_neoplastic"])
    .size()
    .reset_index(name="count")
    .sort_values("count")
)

print(strata_counts)
print(f"\nSmallest stratum: {strata_counts['count'].min()} drugs")
print(f"Strata with <15 drugs: {(strata_counts['count'] < 15).sum()} / {len(strata_counts)} populated")
print(f"Missing strata (0 drugs): {8 - len(strata_counts)} of 8 combinations")

# --- Step 1: Create balanced 100-drug evaluation set ---

def balanced_multilabel_subset_fixedN(
    df,
    target_size=100,
    label_cols=("approved", "immunotherapy", "anti_neoplastic"),
    group_col="drug_name",
    random_state=42,
):
    rng = np.random.default_rng(random_state)

    # 1) One row per drug_name
    d = df.drop_duplicates(subset=[group_col], keep="first").copy()

    # 2) Coerce to booleans robustly (handles True/False or "True"/"False")
    for c in label_cols:
        if d[c].dtype != bool:
            d[c] = (
                d[c]
                .astype(str)
                .str.strip()
                .str.lower()
                .map({"true": True, "false": False, "1": True, "0": False, "yes": True, "no": False})
            )
            d[c] = d[c].fillna(False).astype(bool)

    # 3) Feasibility for 50/50 at N
    desired_pos = target_size // 2  # 50
    pos = {c: int(d[c].sum()) for c in label_cols}
    neg = {c: len(d) - pos[c] for c in label_cols}

    # If any label lacks enough positives or negatives to hit 50/50, reduce target for that label
    feasible_pos = {c: min(desired_pos, pos[c], neg[c] if False else desired_pos) for c in label_cols}
    # (the neg[c] check is only relevant if you tried to enforce negatives directly; we enforce via fill)

    # 4) Greedy selection of positive-covering rows up to per-label targets
    deficits = {c: feasible_pos[c] for c in label_cols}  # how many positives still needed per label
    pool = d.sample(frac=1.0, random_state=random_state)  # shuffle once
    selected_idx = []
    remaining = pool.copy()

    while sum(deficits.values()) > 0 and len(selected_idx) < target_size:
        # Candidates that don't overshoot any label (i.e., don't have True for a label already satisfied)
        allowed = remaining.copy()
        for c in label_cols:
            if deficits[c] == 0:
                allowed = allowed[~allowed[c]]

        if allowed.empty:
            # Can't add any more positive rows without overshoot; break and accept near-50/50
            break

        # Score: number of unmet labels a row would satisfy (weighted by remaining deficits)
        score = np.zeros(len(allowed), dtype=int)
        for i, (_, row) in enumerate(allowed.iterrows()):
            s = 0
            for c in label_cols:
                if deficits[c] > 0 and bool(row[c]):
                    s += deficits[c]  # weight by how "needed" this label is
            score[i] = s

        best_idx = allowed.index[np.argmax(score)]
        best = allowed.loc[best_idx]
        selected_idx.append(best_idx)

        # Update deficits
        for c in label_cols:
            if best[c] and deficits[c] > 0:
                deficits[c] -= 1

        # Remove from pool
        remaining = remaining.drop(index=best_idx)

    # 5) Fill to target_size with all-false rows (keeps negatives high to approach 50/50)
    needed = target_size - len(selected_idx)
    if needed > 0:
        all_false_mask = ~remaining[list(label_cols)].any(axis=1)
        all_false = remaining[all_false_mask]
        if len(all_false) >= needed:
            fill = all_false.sample(n=needed, random_state=random_state)
        else:
            # If not enough all-false, fill with whatever remains (will slightly deviate from 50/50)
            rest = remaining.sample(n=needed - len(all_false), random_state=random_state)
            fill = pd.concat([all_false, rest], axis=0)
        subset = pd.concat([d.loc[selected_idx], fill], axis=0)
    else:
        # Already at N rows
        subset = d.loc[selected_idx]

    # 6) Shuffle final subset for neutrality
    subset = subset.sample(frac=1.0, random_state=random_state).reset_index(drop=True)

    # 7) Report achieved balance
    achieved = subset[list(label_cols)].mean().to_frame("positive_rate")
    counts = subset[list(label_cols)].sum().astype(int).to_frame("true_count")
    summary = achieved.join(counts)

    info = {
        "target_size": target_size,
        "selected_rows": len(subset),
        "feasible_pos_per_label": feasible_pos,
        "remaining_deficits_after_selection": deficits,
    }
    return subset, summary, info

# ---- Use it ----
label_cols = ["approved", "immunotherapy", "anti_neoplastic"]
subset_100, balance_summary, info = balanced_multilabel_subset_fixedN(
    drugs_df,
    target_size=100,
    label_cols=label_cols,
    group_col="drug_name",
    random_state=7
)

subset_100.to_csv('data/eval_drugs_balanced_100.csv')

print(info)
print(balance_summary.T)

# --- Step 2: Re-verify labels via live DGIdb GraphQL API ---

HEADERS = {"Content-Type": "application/json"}

#remaking the labels based on the API results
drug_info_query = """
query drugs_info($names: [String!]) {
    drugs(names: $names) {
        nodes {
            name
            conceptId
            approved
            immunotherapy
            antiNeoplastic
        }
    }
}
"""

def run_query(query, variables: dict) -> dict:
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

df = pd.read_csv("data/eval_drugs_balanced_100.csv")

#df = df[df['drug_name']=='OLAPARIB']

updated_dict = {'index': [], 'drug_name':[], 'concept_id':[],'approved':[],'immunotherapy':[],'anti_neoplastic':[],'source_db_name':[],'source_db_version':[]}

for _, row in df.iterrows():
    index, drug_name, concept_id, approved, immunotherapy, anti_neoplastic, source_db_name, source_db_version = row[['index', 'drug_name', 
                                        'concept_id', 'approved', 'immunotherapy', 'anti_neoplastic', 'source_db_name', 'source_db_version']].values
    

    resp = run_query(drug_info_query, {'names':[drug_name]})

    nodes = resp['data']['drugs']['nodes']
    #print(nodes)

    if len(nodes) == 0:
        actual_approved = False
        actual_immunotherapy = False
        actual_anti_neoplastic = False
        actual_concept_id = None

    else:
        found = False
        for node in nodes:
            if drug_name.upper() == node['name'].upper():
                actual_approved = node['approved']
                actual_immunotherapy = node['immunotherapy']
                actual_anti_neoplastic = node['antiNeoplastic']
                actual_concept_id = node['conceptId']
                found = True
                break
        
        if not found:
            print('never found str match')
            actual_approved = nodes[0]['approved']
            actual_immunotherapy = nodes[0]['immunotherapy']
            actual_anti_neoplastic = nodes[0]['antiNeoplastic']
            actual_concept_id = nodes[0]['conceptId']

    print(drug_name, approved, actual_approved)
    #print(drug_name, concept_id, actual_concept_it)

    updated_dict['index'].append(index)
    updated_dict['drug_name'].append(drug_name)
    updated_dict['concept_id'].append(actual_concept_id)
    updated_dict['approved'].append(actual_approved)
    updated_dict['immunotherapy'].append(actual_immunotherapy)
    updated_dict['anti_neoplastic'].append(actual_anti_neoplastic)
    updated_dict['source_db_name'].append(source_db_name)
    updated_dict['source_db_version'].append(source_db_version)

df = pd.DataFrame(updated_dict)
df.to_csv('data/eval_drugs_ground_truth.csv')