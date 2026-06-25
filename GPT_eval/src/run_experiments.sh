#!/usr/bin/env bash
# Run the full drug classification evaluation pipeline in order.
#
# Data files produced by this pipeline:
#   data/eval_drugs_balanced_100.csv          — 100-drug balanced eval set
#   data/eval_drugs_ground_truth.csv          — labels re-verified via DGIdb API
#   data/predictions_mcp_dgidb_prompt.csv     — Experiment 1: MCP + DGIdb-specific prompt
#   data/predictions_no_mcp_dgidb_prompt.csv  — Experiment 2: No MCP + DGIdb-specific prompt
#   data/predictions_mcp_generic_prompt.csv   — Experiment 3: MCP + generic prompt
#   data/MCP_drug_info_GPT-5/                    — per-drug logs for Experiment 1
#   data/no_MCP_drug_info_GPT-5/                 — per-drug logs for Experiment 2
#   data/MCP_prop_used_drug_info_general_chatbot/ — per-drug logs for Experiment 3
#   data/metrics_mcp_tool_usage_by_prompt.csv — MCP tool-usage analysis for Experiment 3

set -euo pipefail

cd GPT_eval/src

# Pin to the package versions used for the published evaluation
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

# ── Dataset creation ──────────────────────────────────────────────────────────
# Uncomment to regenerate the evaluation set from scratch.
# Warning: this re-queries the DGIdb API for all 100 drugs and will overwrite
# eval_drugs_balanced_100.csv and eval_drugs_ground_truth.csv.

# python create_drug_info_eval_dataset.py
# python create_joint_task_dataset.py

export OPENAI_API_KEY='[INSERT YOUR API KEY]'

# ── Check if API key has been set ─────────────────────────────────────────────
if [[ -n "$OPENAI_API_KEY" && "$OPENAI_API_KEY" != '[INSERT YOUR API KEY]' ]]; then

    # ── Drug classification experiments ──────────────────────────────────────
    python GPT_MCP_drug_info.py           # Experiment 1: MCP + DGIdb-specific prompt
    python GPT_no_MCP_drug_info.py        # Experiment 2: No MCP + DGIdb-specific prompt
    python GPT_MCP_drug_info_prop_used.py # Experiment 3: MCP + generic prompt (measures tool-use rate)

    # ── DGIdb + CIViC joint task: LLM-assisted drug candidate ranking ─────────
    # Requires host_dgidb_civic_MCP.py in the working directory (launched as a
    # subprocess by joint_task_rank_GPT_MCP.py via `fastmcp run`).
    python joint_task_rank_GPT_MCP.py
    python joint_task_rank_GPT_no_MCP.py

else
    echo "No API key specified — skipping experiments and evaluating previous outputs."
fi

# ── Evaluate all experiments ──────────────────────────────────────────────
python eval_drug_info.py              # Reports metrics for all 3 experiments + MCP usage analysis


# ── Evaluate joint task ───────────────────────────────────────────────────────
python eval_rank_joint.py