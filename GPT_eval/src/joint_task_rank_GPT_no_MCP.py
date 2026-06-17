import os, functools, csv, re
import pandas as pd
from datetime import datetime
from pathlib import Path
from typing import List
import openai

def extract_list_from_file(path: str, marker: str = "### LLM OUTPUT ###") -> List[str]:
    """
    Read `path`, find the last occurrence of `marker`, and parse the first non-empty
    line after it as a comma/whitespace-separated list of gene symbols.

    Returns an ordered, de-duplicated list of uppercase symbols (e.g., ['ASNS', 'BRAF']).
    If the marker isn't found or no genes are present, returns [].
    """
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    i = text.rfind(marker)
    if i == -1:
        return []

    tail = text[i + len(marker):]
    # Take the first non-empty line after the marker
    line = next((ln.strip() for ln in tail.splitlines() if ln.strip()), "")
    if not line:
        return []

    # Split on commas or whitespace, clean tokens, keep order & dedupe
    tokens = re.split(r"[,\s]+", line)
    vals, seen = [], set()
    for tok in tokens:
        sym = re.sub(r"[^A-Za-z0-9\-]", "", tok).upper()
        if sym and sym not in seen:
            vals.append(sym)
            seen.add(sym)
    return vals


# with open('open_ai_key.txt', 'r') as file:
#     OPENAI_API_KEY = file.read().rstrip()

#os.environ["OPENAI_API_KEY"] = OPENAI_API_KEY
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")


# === Main function to run a single query and log it ===
def run_query_with_logging(run_id, messages):
    run_id = run_id.replace(':', '_').replace('*', '').replace('/', '_')
    log_dir = "data/no_MCP_ranked_joint_task_GPT-5"
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, f"{run_id}.txt")

    # Open log file and redirect print/logging to it
    with open(log_path, "w", encoding="utf-8") as log_file:

        # Redirect print
        print = functools.partial(__builtins__.print, file=log_file, flush=True)  # noqa: F821

        print(f"### RUN ID ###\n{run_id}\n")
        print("### MESSAGES ###")
        print(messages)
        print()

        client = openai.OpenAI(api_key=OPENAI_API_KEY)

        response = client.chat.completions.create(
            model="gpt-5-2025-08-07",
            messages=messages,
            temperature=1,
        )

        msg = response.choices[0].message.content

        print("### LLM OUTPUT ###")
        print(msg)

        text = str(msg or "").strip()

        # Use the last non-empty line (in case the model adds extra lines)
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        line = lines[-1] if lines else ""

        # Split into tokens (preserve duplicates in full_list)
        tokens = line.split(',')

        # Full list (no dedup)
        full_list = [tok.strip().upper() for tok in tokens if tok.strip()]

        # Deduplicated list (order-preserving)
        seen = set()
        dedup_list = []
        for item in full_list:
            if item not in seen:
                seen.add(item)
                dedup_list.append(item)

    # Return the parsed list (e.g., genes or drugs) and the log path if needed
    return dedup_list, full_list


# this df has the list of MP, but not the actual per gene
# these are per_gene_ranked_joint_task_drug_lists_dataset.csv
df = pd.read_csv('data/subset_50_ranked_joint_task_drug_lists_dataset.csv')

OUT_PATH = "data/per_gene_GPT_no_MCP_joint_rank.csv"

# If file doesn't exist, create it with the correct header
if not os.path.exists(OUT_PATH):
    header = [
        "molecularProfile_name",
        "disease_name",
        "therapies",
        "evidenceType",
        "significance",
        "LLM_CIViC_genes",
        "LLM_DGIdb_drugs",
    ]
    pd.DataFrame(columns=header).to_csv(OUT_PATH, index=False)

out_df = pd.read_csv(OUT_PATH)
df = df[~df['molecularProfile_name'].isin(out_df['molecularProfile_name'])]

print(len(df))

# save genes as file names to indicate they have already been used
FIELDNAMES = [
    'molecularProfile_name',
    'disease_name',
    'therapies',
    'evidenceType',
    'significance',
    'LLM_CIViC_genes',
    'LLM_DGIdb_drugs'
]

with open(OUT_PATH, "a", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    for _, row in df.iterrows():

        molecularProfile_name, disease_name, therapies, evidenceType, significance = row[
            ['molecularProfile_name', 'disease_name', 'therapies', 'evidenceType', 'significance']
        ]
        print(molecularProfile_name)

        messages = [
            {
                "role": "system",
                "content": (
                    "You are a chatbot for CIViC and the Drug Gene Interaction Database (DGIdb). "
                    "When listing genes from CIViC, you must rank them using CIViC evidence strength."
                    "\n\nGENE RANKING RULES:\n"
                    "1. Rank genes first by evidenceLevel (A strongest → E weakest).\n"
                    "2. Within the same evidenceLevel, rank genes by evidenceRating (5 → 1 → null).\n"
                    "3. Return only HGNC gene symbols, in this ranked order, as a comma-separated list "
                    "like: GENE1,GENE2,GENE3."
                )
            },
            {
                "role": "user",
                "content":
                    f"According to CIViC, list genes whose variants/molecular profiles are associated with resistance to {therapies} in {disease_name} "
                    "(significance = RESISTANCE)?\n"
            }
        ]

        full_gene_list, dedup_gene_list = run_query_with_logging(molecularProfile_name, messages)
        print(dedup_gene_list)
        if not dedup_gene_list:
            print('no gene list', molecularProfile_name)
            writer.writerow([molecularProfile_name, disease_name, therapies, evidenceType, significance, '', ''])
            continue

        # Unsure if this rule is needed
        # "4. If multiple entries share the same drug name, treat them as one drug and use the "
        #     "highest interactionScore for ranking.\n"

        for gene in dedup_gene_list:

            messages = [
                {
                    "role": "system",
                    "content": (
                        "You are a chatbot for the CIViC and Drug Gene Interaction Database (DGIdb). "
                        "You receive a gene and DGIdb interaction for that gene. "
                        "Your task is to return a ranked list of up to 20 drugs."
                        "\n\n"
                        "RANKING RULES:\n"
                        "1. List FDA-approved drugs first.\n"
                        "2. Within the approved drugs, sort by the DGIdb interactionScore in descending order "
                        "(higher interactionScore = higher priority).\n"
                        "3. After all approved drugs, list non-approved drugs, also sorted by interactionScore "
                        "in descending order.\n"
                        "4. Return only drug names (no scores), as a comma-separated list like: DRUG1,DRUG2."
                    )
                },
                {
                    "role": "user",
                    "content":
                        f"What drugs target or modulate this gene product in DGIdb? Gene: {gene} \n"
                        "Return a comma-separated list of up to 20 drug names following the ranking rules."
                }
            ]

            full_drug_list, dedup_drug_list = run_query_with_logging(molecularProfile_name + f'_dgidb_{gene}', messages)
            print(full_drug_list)

            writer.writerow([
                molecularProfile_name,
                disease_name,
                therapies,
                evidenceType,
                significance,
                gene,
                ','.join(full_drug_list)
            ])
        
        #break