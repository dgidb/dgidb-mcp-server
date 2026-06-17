"""
Evaluates all three drug classification experiments and reports metrics.

Experiments:
  1. MCP — DGIdb prompt        (predictions_mcp_dgidb_prompt.csv)
  2. No MCP — DGIdb prompt     (predictions_no_mcp_dgidb_prompt.csv)
  3. MCP — Generic prompt      (predictions_mcp_generic_prompt.csv)

Also computes MCP tool-usage analysis for experiment 3:
  - Reads MCP_prop_used_drug_info_general_chatbot/ to detect which calls invoked the MCP tool
  - Reports precision/recall/F1 split by whether the tool was actually used
  - Saves data/metrics_mcp_tool_usage_by_prompt.csv
"""

import numpy as np
import pandas as pd
import os, re, csv, io
from sklearn.metrics import (
    precision_recall_fscore_support,
    precision_score, recall_score, f1_score,
)

LABELS = ['approved', 'immunotherapy', 'anti_neoplastic']

# ── Shared helpers ────────────────────────────────────────────────────────────

def _coerce_bool(s: pd.Series) -> np.ndarray:
    """Coerce booleans or 'True'/'False' strings to a bool numpy array."""
    if s.dtype == bool or s.dtype == np.bool_:
        return s.to_numpy()
    return s.astype(str).str.strip().str.lower().map({'true': True, 'false': False}).to_numpy()

def _normalize_yes_no(val: str):
    s = str(val).strip().strip('"').strip("'").upper()
    if s in ("YES", "Y", "TRUE", "T", "1"):
        return True
    if s in ("NO", "N", "FALSE", "F", "0"):
        return False
    return None

def _extract_yes_no_tuple(text: str):
    """Parse the 2-line CSV block from a log file. Returns a bool 3-tuple or None."""
    if not text:
        return None
    header_pat = r'(?im)^\s*fda_approval\s*,\s*immunotherapy\s*,\s*anti_neoplastic\s*$'
    matches = list(re.finditer(header_pat, text))
    if not matches:
        return None
    tail = text[matches[-1].end():]
    m_vals = re.search(r'^\s*([^\r\n]+?)\s*$', tail, flags=re.MULTILINE)
    if not m_vals:
        return None
    try:
        row = next(csv.reader(io.StringIO(m_vals.group(1))))
    except Exception:
        return None
    if len(row) != 3:
        return None
    norm = tuple(_normalize_yes_no(v) for v in row)
    if any(v is None for v in norm):
        return None
    return norm

# ── Per-experiment metrics ────────────────────────────────────────────────────

def metrics_for_df(df: pd.DataFrame, labels=LABELS):
    rows = []
    y_true_all, y_pred_all = [], []

    for lbl in labels:
        y_true = _coerce_bool(df[lbl])
        y_pred = _coerce_bool(df[f'pred_{lbl}'])

        fp_drugs = df.loc[(~y_true) & (y_pred), 'drug_name'].tolist()
        fn_drugs = df.loc[(y_true) & (~y_pred), 'drug_name'].tolist()
        print(f"  [{lbl}] False positives: {fp_drugs}")
        print(f"  [{lbl}] False negatives: {fn_drugs}")

        tp = int(((y_true == True)  & (y_pred == True)).sum())
        fp = int(((y_true == False) & (y_pred == True)).sum())
        fn = int(((y_true == True)  & (y_pred == False)).sum())
        tn = int(((y_true == False) & (y_pred == False)).sum())

        prec, rec, f1, _ = precision_recall_fscore_support(
            y_true, y_pred, average='binary', pos_label=True, zero_division=0
        )
        support_pos = int((y_true == True).sum())

        if support_pos == 0:
            prec = rec = f1 = float('nan')

        rows.append({
            'label': lbl,
            'TP': tp, 'FP': fp, 'FN': fn, 'TN': tn,
            'precision':        (None if np.isnan(prec) else round(prec, 3)),
            'recall':           (None if np.isnan(rec)  else round(rec,  3)),
            'f1':               (None if np.isnan(f1)   else round(f1,   3)),
            'support_positives': support_pos,
        })
        y_true_all.extend(y_true.tolist())
        y_pred_all.extend(y_pred.tolist())

    per_label = pd.DataFrame(rows).set_index('label')

    micro_p, micro_r, micro_f1, _ = precision_recall_fscore_support(
        y_true_all, y_pred_all, average='micro', zero_division=0
    )
    micro_avg = {
        'precision': round(micro_p,   3),
        'recall':    round(micro_r,   3),
        'f1':        round(micro_f1,  3),
        'support_total_positives': int(per_label['support_positives'].sum()),
    }

    pos_mask = per_label['support_positives'] > 0
    if pos_mask.any():
        macro_avg = {
            'precision':     round(per_label.loc[pos_mask, 'precision'].mean(), 3),
            'recall':        round(per_label.loc[pos_mask, 'recall'].mean(),    3),
            'f1':            round(per_label.loc[pos_mask, 'f1'].mean(),        3),
            'num_labels_used': int(pos_mask.sum()),
        }
    else:
        macro_avg = {'precision': None, 'recall': None, 'f1': None, 'num_labels_used': 0}

    return per_label, micro_avg, macro_avg

# ── Run all three experiments ─────────────────────────────────────────────────

EXPERIMENTS = [
    ('MCP — DGIdb prompt',    'data/predictions_mcp_dgidb_prompt.csv'),
    ('No MCP — DGIdb prompt', 'data/predictions_no_mcp_dgidb_prompt.csv'),
    ('MCP — Generic prompt',  'data/predictions_mcp_generic_prompt.csv'),
]

for name, csv_path in EXPERIMENTS:
    print(f"\n{'='*60}")
    print(f"Experiment: {name}")
    print('='*60)
    df = pd.read_csv(csv_path)
    per_label, micro_avg, macro_avg = metrics_for_df(df)
    print("\nPer-label metrics:")
    print(per_label)
    print("\nMicro average across all decisions:")
    print(micro_avg)
    print("\nMacro average over labels with positive support:")
    print(macro_avg)

# ── MCP tool-usage analysis (generic prompt experiment only) ──────────────────

print(f"\n{'='*60}")
print("MCP Tool Usage Analysis — Generic Prompt Experiment")
print("(split by whether the model actually invoked the MCP tool)")
print('='*60)

ground_truth_df = pd.read_csv('data/eval_drugs_ground_truth.csv')
log_dir = 'data/MCP_prop_used_drug_info_general_chatbot'

rows = []
for fname in os.listdir(log_dir):
    if not fname.endswith('.txt'):
        continue
    with open(os.path.join(log_dir, fname), 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    mcp_used = 'Invoking MCP tool' in content
    result = _extract_yes_no_tuple(content)
    if result is None:
        continue
    pred_approved, pred_immunotherapy, pred_anti_neoplastic = result

    key = fname.replace('.txt', '')
    match = ground_truth_df[
        ground_truth_df['drug_name'].str.replace(':', '_', regex=False) == key.replace(':', '_')
    ]
    if match.empty:
        continue
    gt = match.iloc[0]
    approved, immunotherapy, anti_neoplastic = gt[['approved', 'immunotherapy', 'anti_neoplastic']].tolist()

    rows.extend([
        {'label': 'approved',        'y_true': bool(approved),        'y_pred': bool(pred_approved),        'mcp_used': mcp_used},
        {'label': 'immunotherapy',   'y_true': bool(immunotherapy),   'y_pred': bool(pred_immunotherapy),   'mcp_used': mcp_used},
        {'label': 'anti_neoplastic', 'y_true': bool(anti_neoplastic), 'y_pred': bool(pred_anti_neoplastic), 'mcp_used': mcp_used},
    ])

usage_df = pd.DataFrame(rows)
if usage_df.empty:
    print(f"No log files found in {log_dir}; skipping tool-usage analysis.")
else:
    # Count total drugs per subset (same for all labels within a subset)
    n_mcp    = int((usage_df['mcp_used'] == True).sum()  // len(LABELS))
    n_no_mcp = int((usage_df['mcp_used'] == False).sum() // len(LABELS))

    metrics_rows = []
    for (mcp_used, label), sub in usage_df.groupby(['mcp_used', 'label']):
        y_true = sub['y_true'].astype(bool)
        y_pred = sub['y_pred'].astype(bool)
        n_sub  = n_mcp if mcp_used else n_no_mcp
        n_pos  = int(y_true.sum())
        metrics_rows.append({
            'label':      label,
            'mcp_used':   'MCP' if mcp_used else 'No MCP',
            'n_positive': n_pos,
            'n_subset':   n_sub,
            # Fraction string for display: "10/45" makes the denominator self-evident
            'N':          f"{n_pos}/{n_sub}",
            'precision':  precision_score(y_true, y_pred, pos_label=True, zero_division=0),
            'recall':     recall_score(y_true, y_pred,    pos_label=True, zero_division=0),
            'f1':         f1_score(y_true, y_pred,        pos_label=True, zero_division=0),
        })

    summary = pd.DataFrame(metrics_rows).sort_values(['label', 'mcp_used'])

    # Numeric pivot saved to CSV (keeps n_positive and n_subset as separate columns)
    numeric_pivot = (
        summary
        .pivot(index='label', columns='mcp_used',
               values=['n_positive', 'n_subset', 'precision', 'recall', 'f1'])
        .reindex(columns=pd.MultiIndex.from_product(
            [['n_positive', 'n_subset', 'precision', 'recall', 'f1'], ['MCP', 'No MCP']]
        ))
    )
    numeric_pivot.index = numeric_pivot.index.map({
        'approved':        'FDA-Approved',
        'immunotherapy':   'Immunotherapy',
        'anti_neoplastic': 'Antineoplastic',
    })
    numeric_pivot.round(6).to_csv('data/metrics_mcp_tool_usage_by_prompt.csv')

    # Display pivot uses the fraction string N = "n_positive/n_subset" for readability
    display_pivot = (
        summary
        .pivot(index='label', columns='mcp_used',
               values=['N', 'precision', 'recall', 'f1'])
        .reindex(columns=pd.MultiIndex.from_product(
            [['N', 'precision', 'recall', 'f1'], ['MCP', 'No MCP']]
        ))
    )
    display_pivot.index = display_pivot.index.map({
        'approved':        'FDA-Approved',
        'immunotherapy':   'Immunotherapy',
        'anti_neoplastic': 'Antineoplastic',
    })

    print(f"\n  MCP Used (N={n_mcp} drugs)  vs  No MCP Usage (N={n_no_mcp} drugs)")
    print("  N = positives/subset-size (numerators sum to total positives across all 100 drugs)\n")
    print(display_pivot.round(4))
    print("\nSaved to data/metrics_mcp_tool_usage_by_prompt.csv")