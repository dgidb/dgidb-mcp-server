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

def sanitize_filename(name: str, replacement: str = "_", max_length: int = 255, ensure_unique: bool = False) -> str:
    """
    Sanitize a string to be used as a safe filename.

    Args:
        name (str): Original filename (without path or extension).
        replacement (str): Replacement for unsafe characters (default: "-").
        max_length (int): Max length for filename (default: 255).
        ensure_unique (bool): Whether to append hash to ensure uniqueness when truncated.

    Returns:
        str: Safe, valid filename.
    """

    name = name.replace('::', '-')
    # Invalid characters for Windows and general use
    invalid_chars = r'[<>:"/\\|?*\n\r\t]'
    name = re.sub(invalid_chars, replacement, name)

    # Collapse multiple replacement chars
    name = re.sub(re.escape(replacement) + r'{2,}', replacement, name)

    # Strip leading/trailing replacement characters
    name = name.strip(replacement)

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
    "args": ["run", "DGIdb_MCP_server.py:mcp", "--transport", "stdio"],
    "env": {
        "OPENAI_API_KEY": OPENAI_API_KEY,
        "MCP_DEBUG": "1",
    },
    "errlog": sys.stderr,
}

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


# === Main function to run a single query and log it ===
async def run_query_with_logging(run_id, query):
    run_id = run_id.replace(':', '_')
    log_dir = "data/MCP_drug_info_GPT-5"
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
        srv = MCPServerStdio(name="DGIdb_MCP", params=server_params, client_session_timeout_seconds=60.0)
        await srv.__aenter__()

        agent = Agent(
            name="DGIdb_MCP",
            instructions = (
                    "Use the tools to answer drug-gene interaction questions for the Drug Gene Interaction Database (DGIdb)."
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

        # --- Extract the three YES/NO answers ---
        answers_tuple = _extract_yes_no_tuple(getattr(result, "final_output", "") or "")
        if answers_tuple is None:
            print("ERROR: Failed to parse 2-line CSV (fda_approval,immunotherapy,anti_neoplastic) from LLM output.")
            # Optionally echo raw output again for debugging
            print("RAW OUTPUT FOLLOWS:")
            print(result.final_output)
            answers_tuple = (None, None, None)

        await srv.__aexit__(None, None, None)

    # Return both the parsed answers and the log path (for run_index.csv logging)
    return answers_tuple, log_path

eval_df = pd.read_csv('data/eval_drugs_ground_truth.csv')

out_csv = 'data/predictions_mcp_dgidb_prompt.csv'
if os.path.exists(out_csv):
    out_df = pd.read_csv(out_csv)
    eval_df = eval_df[~eval_df['drug_name'].isin(out_df['drug_name'])]
    write_header = False
else:
    write_header = True

print(len(eval_df))

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

        answers, log_path = asyncio.run(run_query_with_logging(drug_name, messages))

        pred_approved, pred_immunotherapy, pred_anti_neoplastic = answers

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