"""
Inference script for the generic-prompt MCP experiment.
Answers the question: How often does the model use the MCP tool when
the prompt does NOT name or describe the database?

System prompt used here:
    "You are a helpful AI assistant with access to external tools."
compared to the DGIdb-specific prompt used in GPT_MCP_drug_info.py.

Outputs:
  data/MCP_prop_used_drug_info_general_chatbot/<drug_name>.txt  — per-drug logs (MCP usage detectable here)
  data/predictions_mcp_generic_prompt.csv                       — structured predictions for eval_drug_info.py

On each run the script:
  1. Parses any existing log files not yet written to the predictions CSV.
  2. Runs inference for any drugs that have no log file yet.
"""

import asyncio, os, sys, functools, logging
import pandas as pd
from openai import AsyncOpenAI
from agents import Agent, Runner, OpenAIChatCompletionsModel, ModelSettings
from agents.mcp import MCPServerStdio
import csv
import re, io

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
    Returns a 3-tuple of booleans, or None on failure.
    """
    if not text:
        return None

    header_pat = r'(?im)^\s*fda_approval\s*,\s*immunotherapy\s*,\s*anti_neoplastic\s*$'
    matches = list(re.finditer(header_pat, text))
    if not matches:
        return None

    start = matches[-1].end()
    tail = text[start:]

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
    return norm


class MCPToolFilter(logging.Filter):
    def filter(self, record):
        if record.levelname != "DEBUG":
            return False
        if not record.name.startswith("openai.agents"):
            return False
        msg = record.getMessage()
        return (
            msg.startswith("Invoking MCP tool")
            or msg.startswith("MCP tool")
        )

class NewlineFormatter(logging.Formatter):
    def format(self, record):
        s = super().format(record)
        return s.replace("\\n", "\n")

# with open('open_ai_key.txt', 'r') as file:
#     OPENAI_API_KEY = file.read().rstrip()

# os.environ["OPENAI_API_KEY"] = OPENAI_API_KEY

OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")

server_params = {
    "command": "fastmcp",
    "args": ["run", "DGIdb_MCP_server.py:mcp", "--transport", "stdio"],
    "env": {
        "OPENAI_API_KEY": OPENAI_API_KEY,
        "MCP_DEBUG": "1",
    },
    "errlog": sys.stderr,
}

async def run_query_with_logging(run_id, query):
    run_id = run_id.replace(':', '_')
    log_dir = "data/MCP_prop_used_drug_info_general_chatbot"
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, f"{run_id}.txt")

    with open(log_path, "w", encoding="utf-8") as log_file:

        print = functools.partial(__builtins__.print, file=log_file, flush=True)

        handler = logging.StreamHandler(log_file)
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(NewlineFormatter("%(message)s"))
        handler.addFilter(MCPToolFilter())
        logging.getLogger().handlers = [handler]
        logging.getLogger().setLevel(logging.DEBUG)

        srv = MCPServerStdio(name="DGIdb_MCP", params=server_params, client_session_timeout_seconds=60.0)
        await srv.__aenter__()

        agent = Agent(
            name="agent",
            instructions=(
                "You are a helpful AI assistant with access to external tools. Follow the user's instructions carefully."
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

        answers_tuple = _extract_yes_no_tuple(getattr(result, "final_output", "") or "")
        if answers_tuple is None:
            print("ERROR: Failed to parse 2-line CSV (fda_approval,immunotherapy,anti_neoplastic) from LLM output.")
            print("RAW OUTPUT FOLLOWS:")
            print(result.final_output)
            answers_tuple = (None, None, None)

        await srv.__aexit__(None, None, None)

    return answers_tuple, log_path


# ── Setup ─────────────────────────────────────────────────────────────────────

eval_df = pd.read_csv('data/eval_drugs_ground_truth.csv')
eval_df['_key'] = eval_df['drug_name'].astype(str).str.replace(':', '_', regex=False)

log_dir = "data/MCP_prop_used_drug_info_general_chatbot"
os.makedirs(log_dir, exist_ok=True)

out_csv = 'data/predictions_mcp_generic_prompt.csv'

FIELDNAMES = [
    "index", "drug_name", "approved", "immunotherapy", "anti_neoplastic",
    "pred_approved", "pred_immunotherapy", "pred_anti_neoplastic"
]

# Drugs already written to the predictions CSV
if os.path.exists(out_csv):
    already_in_csv = set(pd.read_csv(out_csv)['drug_name'].astype(str).str.replace(':', '_', regex=False))
    write_header = False
else:
    already_in_csv = set()
    write_header = True

# All log files that exist, keyed by sanitised drug name
existing_logs = {
    os.path.splitext(f)[0]: os.path.join(log_dir, f)
    for f in os.listdir(log_dir)
    if f.endswith(".txt")
}

logs_to_parse   = {k: v for k, v in existing_logs.items() if k not in already_in_csv}
needs_inference = eval_df[~eval_df['_key'].isin(existing_logs.keys())]

print(f"Log files to parse into CSV : {len(logs_to_parse)}")
print(f"Drugs needing inference     : {len(needs_inference)}")

with open(out_csv, "a", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
    if write_header:
        writer.writeheader()

    # ── Step 1: Write predictions from existing log files ─────────────────────
    for key, log_path in logs_to_parse.items():
        with open(log_path, 'r', encoding='utf-8', errors='ignore') as lf:
            content = lf.read()

        answers = _extract_yes_no_tuple(content)
        if answers is None:
            print(f"WARNING: could not parse predictions from {log_path} — skipping.")
            continue

        match = eval_df[eval_df['_key'] == key]
        if match.empty:
            print(f"WARNING: no ground-truth row found for '{key}' — skipping.")
            continue

        row = match.iloc[0]
        pred_approved, pred_immunotherapy, pred_anti_neoplastic = answers
        writer.writerow({
            'index':               row['index'],
            'drug_name':           row['drug_name'],
            'approved':            row['approved'],
            'immunotherapy':       row['immunotherapy'],
            'anti_neoplastic':     row['anti_neoplastic'],
            'pred_approved':       pred_approved,
            'pred_immunotherapy':  pred_immunotherapy,
            'pred_anti_neoplastic': pred_anti_neoplastic,
        })
        print(f"  Parsed from log: {row['drug_name']}")

    # ── Step 2: Run inference for drugs with no log file yet ──────────────────
    for idx, row in needs_inference.iterrows():
        index, drug_name, approved, immunotherapy, anti_neoplastic = row[
            ["index", "drug_name", "approved", "immunotherapy", "anti_neoplastic"]
        ].values
        print(index, drug_name)

        messages = [
            {
                "role": "user",
                "content": f"""
                Drug: {drug_name}
                For the drug named above ONLY, state if each property applies (YES/NO).

                Output EXACTLY this 2-line CSV (no extra text):
                fda_approval,immunotherapy,anti_neoplastic
                <YES/NO>,<YES/NO>,<YES/NO>
                """
            }
        ]

        answers, log_path = asyncio.run(run_query_with_logging(drug_name, messages))
        pred_approved, pred_immunotherapy, pred_anti_neoplastic = answers

        writer.writerow({
            'index':               index,
            'drug_name':           drug_name,
            'approved':            approved,
            'immunotherapy':       immunotherapy,
            'anti_neoplastic':     anti_neoplastic,
            'pred_approved':       pred_approved,
            'pred_immunotherapy':  pred_immunotherapy,
            'pred_anti_neoplastic': pred_anti_neoplastic,
        })