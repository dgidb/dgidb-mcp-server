import asyncio, os, sys, functools, logging
import pandas as pd
from datetime import datetime
from openai import AsyncOpenAI
from agents import Agent, Runner, OpenAIChatCompletionsModel, ModelSettings
from agents.mcp import MCPServerStdio
import csv
import random, re, hashlib
import numpy as np
import time, io
from pathlib import Path
from typing import List

class MCPToolFilter(logging.Filter):
    def filter(self, record):
        if record.levelname != "DEBUG":
            return False
        # only look at the openai.agents loggers
        if not record.name.startswith("openai.agents"):
            return False

        msg = record.getMessage()
        # let through both invocation _and_ return messages
        return (
            msg.startswith("Invoking MCP tool")
            or msg.startswith("MCP tool")
        )

class NewlineFormatter(logging.Formatter):
    def format(self, record):
        # first let the base Formatter build the string
        s = super().format(record)
        # then un‐escape any “\n” sequences into real newlines
        return s.replace("\\n", "\n")
    
# with open('open_ai_key.txt', 'r') as file:
#     OPENAI_API_KEY = file.read().rstrip()

# os.environ["OPENAI_API_KEY"] = OPENAI_API_KEY

OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")

# === Set up the MCP server params (only once) ===
server_params = {
    "command": "fastmcp",
    "args": ["run", "host_dgidb_civic_MCP.py:mcp", "--transport", "stdio"],
    "env": {
        "OPENAI_API_KEY": OPENAI_API_KEY,
        "MCP_DEBUG": "1",
    },
    "errlog": sys.stderr,
}

# === Main function to run a single query and log it ===
async def run_query_with_logging(run_id, query):
    run_id = run_id.replace(':', '_').replace('*', '').replace('/', '_')
    log_dir = "data/MCP_ranked_joint_task_GPT-5"
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, f"{run_id}.txt")

    # Open log file and redirect print/logging to it
    with open(log_path, "w", encoding="utf-8") as log_file:

        # Redirect print
        print = functools.partial(__builtins__.print, file=log_file, flush=True)

        # Set up logging
        handler = logging.StreamHandler(log_file)
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(NewlineFormatter("%(message)s"))
        handler.addFilter(MCPToolFilter())
        logging.getLogger().handlers = [handler]
        logging.getLogger().setLevel(logging.DEBUG)

        # Start MCP server
        srv = MCPServerStdio(name="host_dgidb_civic_MCP", params=server_params, client_session_timeout_seconds=60.0)
        await srv.__aenter__()

        agent = Agent(
            name="host_dgidb_civic_MCP",
            instructions = (
                    "Use the tools to answer precision oncology questions for CIViC and drug-gene interaction questions for the Drug Gene Interaction Database (DGIdb)."
                ),            
            model=OpenAIChatCompletionsModel(
                model="gpt-5-2025-08-07",
                openai_client=AsyncOpenAI()
            ),
            model_settings=ModelSettings(temperature=1),
            mcp_servers=[srv],
        )

        result = await Runner.run(agent, query)

        print("### LLM OUTPUT ###")
        print(result.final_output)

        # Build full_list (with duplicates) and dedup_list (order-preserving uniqueness)
        if isinstance(result.final_output, (list, tuple)):
            # Model returned a Python list
            full_list = [str(x).strip().upper() for x in result.final_output if str(x).strip()]

        else:
            text = str(result.final_output or "").strip()
            # Use the last non-empty line (in case the model adds extra lines)
            lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
            line = lines[-1] if lines else ""

            # Split on commas (preserve duplicates)
            tokens = line.split(',')
            full_list = [tok.strip().upper() for tok in tokens if tok.strip()]

        # Make dedup_list from full_list (preserve order, no duplicates)
        seen = set()
        dedup_list = []
        for item in full_list:
            if item not in seen:
                seen.add(item)
                dedup_list.append(item)

        # Example: ['ASNS', 'BRAF']
        #print("Parsed genes:", gene_list)

        # --- Extract the three YES/NO answers ---
        # answers_tuple = _extract_yes_no_tuple(getattr(result, "final_output", "") or "")
        # if answers_tuple is None:
        #     print("ERROR: Failed to parse 2-line CSV (fda_approval,immunotherapy,anti_neoplastic) from LLM output.")
        #     # Optionally echo raw output again for debugging
        #     print("RAW OUTPUT FOLLOWS:")
        #     print(result.final_output)
        #     answers_tuple = (None, None, None)

        await srv.__aexit__(None, None, None)

    # Return both the parsed answers and the log path (for run_index.csv logging)
    return full_list, dedup_list

#this df has the list of MP, but not the actual per gene
#these are per_gene_ranked_joint_task_drug_lists_dataset.csv
df = pd.read_csv('data/subset_50_ranked_joint_task_drug_lists_dataset.csv')

OUT_PATH = "data/per_gene_GPT_MCP_joint_rank.csv"

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

#save genes as file names to indicate they have already been used
FIELDNAMES = ['molecularProfile_name','disease_name','therapies','evidenceType','significance', 'LLM_CIViC_genes', 'LLM_DGIdb_drugs']

with open(OUT_PATH, "a", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    for _, row in df.iterrows():

        molecularProfile_name,disease_name,therapies,evidenceType,significance = row[['molecularProfile_name','disease_name','therapies','evidenceType','significance']]
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
        
        full_gene_list, dedup_gene_list = asyncio.run(run_query_with_logging(molecularProfile_name, messages))
        if not dedup_gene_list:
            print('no gene list', molecularProfile_name)
            writer.writerow([molecularProfile_name,disease_name,therapies,evidenceType,significance,'', ''])
            continue

        #Unsure if this rule is needed
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
        
            full_drug_list, dedup_drug_list = asyncio.run(run_query_with_logging(molecularProfile_name+f'_dgidb_{gene}', messages))
            print(full_drug_list)

            writer.writerow([molecularProfile_name,disease_name,therapies,evidenceType,significance,gene, ','.join(full_drug_list)])