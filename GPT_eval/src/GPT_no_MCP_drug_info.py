import asyncio, os, sys, functools, logging
import pandas as pd
from datetime import datetime
import openai
from agents import Agent, Runner, OpenAIChatCompletionsModel, ModelSettings
from agents.mcp import MCPServerStdio
import csv
import random, re, hashlib
import numpy as np
import time, io

# with open('open_ai_key.txt', 'r') as file:
#     OPENAI_API_KEY = file.read().rstrip()

# os.environ["OPENAI_API_KEY"] = OPENAI_API_KEY

eval_df = pd.read_csv('data/eval_drugs_ground_truth.csv')

out_csv = 'data/predictions_no_mcp_dgidb_prompt.csv'
if os.path.exists(out_csv):
    out_df = pd.read_csv(out_csv)
    eval_df = eval_df[~eval_df['drug_name'].isin(out_df['drug_name'])]
    write_header = False
else:
    write_header = True

print(len(eval_df))

def _normalize_yes_no(val: str):
    s = str(val).strip().strip('"').strip("'").upper()
    if s in ("YES", "Y", "TRUE", "T", "1"):
        return True
    if s in ("NO", "N", "FALSE", "F", "0"):
        return False
    return None

def _extract_yes_no_tuple(text: str):
    """
    Expect a 2-line CSV anywhere in the output:
    fda_approval,immunotherapy,anti_neoplastic
    <YES/NO>,<YES/NO>,<YES/NO>
    Returns a 3-tuple of normalized 'YES'/'NO' strings, or None on failure.
    """
    if not text:
        return None

    # Find the last occurrence of the exact header, then take the next non-empty line
    header_pat = r'(?im)^\s*fda_approval\s*,\s*immunotherapy\s*,\s*anti_neoplastic\s*$'
    matches = list(re.finditer(header_pat, text))
    if not matches:
        return None

    start = matches[-1].end()  # take the last header match
    tail = text[start:]

    # Next non-empty line after the header is the values line
    m_vals = re.search(r'^\s*([^\r\n]+?)\s*$', tail, flags=re.MULTILINE)
    if not m_vals:
        return None

    values_line = m_vals.group(1)
    try:
        row = next(csv.reader(io.StringIO(values_line)))
    except Exception:
        return None

    if len(row) != 3:
        return None

    norm = tuple(_normalize_yes_no(v) for v in row)
    if any(v is None for v in norm):
        return None
    return norm  # ('YES'/'NO', 'YES'/'NO', 'YES'/'NO')


OUT_DIR = 'data/no_MCP_drug_info_GPT-5'

FIELDNAMES = [
    "index", "drug_name", "approved", "immunotherapy", "anti_neoplastic", "pred_approved", "pred_immunotherapy", "pred_anti_neoplastic"
]

with open(out_csv, "a", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
    if write_header:
        writer.writeheader()
    for idx, row in eval_df.iterrows():
        index, drug_name, approved, immunotherapy, anti_neoplastic = row[["index", "drug_name", "approved", "immunotherapy", "anti_neoplastic"]].values
        print(index, drug_name)
        messages = [
            {
                "role": "system",
                "content": (
                    "You are a chatbot for the Drug Gene Interaction Database (DGIdb)"
                )
            },
            {
                "role": "user",
                "content": 
                    f"""
                    Drug: {drug_name}

                    For the drug named above ONLY, state if each property applies (YES/NO) based on DGIdb.

                    Output EXACTLY this 2-line CSV (no extra text):
                    fda_approval,immunotherapy,anti_neoplastic
                    <YES/NO>,<YES/NO>,<YES/NO>
                    """
            }
        ]

        api_key = os.getenv("OPENAI_API_KEY")
        # Initialize client (replace with your own API key)
        client = openai.OpenAI(api_key=api_key)

        response = client.chat.completions.create(
                    model="gpt-5-2025-08-07",                    
                    messages=messages,
                    temperature=1
                )

        msg = response.choices[0].message.content
        print(msg)

        with open(os.path.join(OUT_DIR, drug_name.replace(':', '_')), 'w') as file:
            file.write(msg)

        answers_tuple = _extract_yes_no_tuple(msg)
        if answers_tuple is None:
            print("ERROR: Failed to parse 2-line CSV (fda_approval,immunotherapy,anti_neoplastic) from LLM output.")
            # Optionally echo raw output again for debugging
            print("RAW OUTPUT FOLLOWS:")
            print(msg)
            answers_tuple = (None, None, None)

        pred_approved, pred_immunotherapy, pred_anti_neoplastic = answers_tuple

        row_dict = {
            'index':index,
            'drug_name':drug_name,
            'approved':approved,
            'immunotherapy': immunotherapy,
            "anti_neoplastic": anti_neoplastic, 
            "pred_approved": pred_approved, 
            "pred_immunotherapy": pred_immunotherapy,
            "pred_anti_neoplastic": pred_anti_neoplastic
        }

        writer.writerow(row_dict)
        time.sleep(10)