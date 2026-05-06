"""
build_comparison_table.py — Compare QSAR ranking (input Excel row order) with
AutoDock Vina ranking (best = most negative score) for the top 10 CIDs in each
of the three Degree columns, per protein.

Inputs:
  - INPUT_DIR/BMP1_drug_screening_results.xlsx     (Degree 1, Degree 2, Degree 3)
  - INPUT_DIR/MMP9_drug_screening_results.xlsx     (Degree 1, Degree 2, Degree 3)
  - DOCK_DIR/runs/<TARGET>/<TARGET>_docking_results.xlsx
        (must contain columns: cid, vina_score, status; produced by
         dock_top_drugs.py)

Output for each protein (under DOCK_DIR/comparison/):
  <TARGET>_qsar_vs_vina.xlsx with 9 columns:
      Deg1_QSAR  Deg1_Vina  Deg1_Vina_Score
      Deg2_QSAR  Deg2_Vina  Deg2_Vina_Score
      Deg3_QSAR  Deg3_Vina  Deg3_Vina_Score
  10 rows per file (top-10 CIDs per degree).

  DegN_QSAR        : the CID at row k in the input Excel (1-indexed, 1=best)
  DegN_Vina        : the same set of CIDs re-sorted by Vina score (most
                     negative first); CIDs without a usable score keep their
                     CID written here but are placed at the bottom of the rank.
  DegN_Vina_Score  : the Vina score (kcal/mol) of the CID in DegN_Vina at the
                     same row, or 'NO_SMILES' / 'FAILED' / 'NOT_DOCKED' if
                     the docking did not produce a score.

A summary sheet 'rank_changes' in each workbook shows, per CID, how its rank
shifted between QSAR and Vina (positive = improved, negative = worsened).

Usage:
  conda activate sbdd-env
  python build_comparison_table.py
"""

import os
import sys
from pathlib import Path
import pandas as pd


# =========================================================================
# CONFIG
# =========================================================================
INPUT_DIR = Path(os.environ.get("T2D_INPUT_DIR", "./drug_screening_results"))
DOCK_DIR  = Path(os.environ.get("T2D_DOCK_DIR",  "./drug_screening_results/docking_out"))
OUT_DIR   = DOCK_DIR / "comparison"

TARGETS = {
    "BMP1": INPUT_DIR / "BMP1_drug_screening_results.xlsx",
    "MMP9": INPUT_DIR / "MMP9_drug_screening_results.xlsx",
}

DEGREE_COLUMNS = ["Degree 1", "Degree 2", "Degree 3"]
TOP_N = 10


# =========================================================================
# HELPERS
# =========================================================================
def load_dock_results(target):
    """Load the docking results produced by dock_top_drugs.py.
    Returns a dict: cid (int) -> (vina_score_or_None, status)."""
    path = DOCK_DIR / "runs" / target / f"{target}_docking_results.xlsx"
    if not path.exists():
        sys.exit(f"ERROR: docking results not found for {target}: {path}")
    df = pd.read_excel(path)
    # Normalise column names just in case
    cols = {c.lower(): c for c in df.columns}
    cid_col   = cols.get("cid")
    score_col = cols.get("vina_score")
    status_col = cols.get("status")
    if cid_col is None or score_col is None:
        sys.exit(f"ERROR: {path} missing 'cid' or 'vina_score' column")
    out = {}
    for _, row in df.iterrows():
        try:
            cid = int(row[cid_col])
        except (ValueError, TypeError):
            continue
        score = row[score_col]
        status = row[status_col] if status_col else "unknown"
        if pd.isna(score):
            score = None
        else:
            score = float(score)
        out[cid] = (score, str(status) if not pd.isna(status) else "unknown")
    return out


def vina_label_for(cid, score, status):
    """Return what to write in the Vina_Score cell for a CID with no number.
    Used only when score is None (the dock produced no result)."""
    if status is None or status == "unknown":
        return "NOT_DOCKED"
    s = status.lower()
    if "ligand_prep" in s or "embed" in s:
        return "NO_SMILES"   # close enough — couldn't make a 3D structure
    if "vina_failed" in s or "parse_failed" in s or "timeout" in s or "exception" in s:
        return "FAILED"
    return "NOT_DOCKED"


def build_degree_pair(qsar_cids, dock_results):
    """For one degree column, build aligned (qsar_list, vina_list, score_list)
    of length TOP_N.

    qsar_list  : the CIDs in input-file row order (1=best)
    vina_list  : the SAME CIDs re-sorted by Vina score (most negative = best);
                 CIDs without a numeric score sink to the bottom in their
                 original QSAR order
    score_list : the Vina score (or text label) for each CID in vina_list
    """
    # Cap to TOP_N (defensive)
    qsar_cids = list(qsar_cids)[:TOP_N]

    # Pad to TOP_N if input file has fewer rows than expected
    while len(qsar_cids) < TOP_N:
        qsar_cids.append(None)

    # Build (cid, score, label) for re-sorting by Vina
    enriched = []
    for cid in qsar_cids:
        if cid is None:
            enriched.append((None, None, ""))
            continue
        score, status = dock_results.get(cid, (None, "not_in_results"))
        if score is None:
            enriched.append((cid, None, vina_label_for(cid, score, status)))
        else:
            enriched.append((cid, score, score))

    # Sort: numeric scores first by ascending score (most negative = best),
    # then non-numeric (failed) entries in their original QSAR order at the end
    scored = [e for e in enriched if isinstance(e[1], float)]
    failed = [e for e in enriched if not isinstance(e[1], float)]
    scored.sort(key=lambda e: e[1])
    vina_sorted = scored + failed

    qsar_list  = [c if c is not None else "" for c in qsar_cids]
    vina_list  = [(c if c is not None else "") for c, _, _ in vina_sorted]
    score_list = []
    for c, s, label in vina_sorted:
        if c is None:
            score_list.append("")
        elif isinstance(s, float):
            score_list.append(round(s, 3))
        else:
            score_list.append(label)
    return qsar_list, vina_list, score_list


def build_rank_changes(qsar_cids_by_degree, vina_cids_by_degree):
    """Build a long-format summary of rank shifts.
    For each (degree, cid), record qsar_rank, vina_rank, rank_diff.
    rank_diff = qsar_rank - vina_rank   (positive = climbed, negative = fell)"""
    rows = []
    for deg, qsar_list in qsar_cids_by_degree.items():
        vina_list = vina_cids_by_degree[deg]
        # Map cid -> 1-based rank in each list (skip empty cells)
        qsar_rank = {c: i + 1 for i, c in enumerate(qsar_list) if c != ""}
        vina_rank = {c: i + 1 for i, c in enumerate(vina_list) if c != ""}
        for cid, qrank in qsar_rank.items():
            vrank = vina_rank.get(cid)
            rows.append({
                "degree":   deg,
                "cid":      cid,
                "qsar_rank": qrank,
                "vina_rank": vrank,
                "rank_diff": (qrank - vrank) if vrank is not None else None,
            })
    return pd.DataFrame(rows)


# =========================================================================
# MAIN
# =========================================================================
def process_target(target, xlsx_path):
    print(f"\n=== {target} ===")
    if not xlsx_path.exists():
        sys.exit(f"ERROR: input not found: {xlsx_path}")

    df_in = pd.read_excel(xlsx_path)

    dock_results = load_dock_results(target)
    print(f"  loaded {len(dock_results)} docking results")

    # Per-degree QSAR top-10 (preserve input order)
    qsar_by_deg = {}
    vina_by_deg = {}
    score_by_deg = {}
    for col in DEGREE_COLUMNS:
        if col not in df_in.columns:
            print(f"  ! column '{col}' missing; skipping")
            qsar_by_deg[col] = [""] * TOP_N
            vina_by_deg[col] = [""] * TOP_N
            score_by_deg[col] = [""] * TOP_N
            continue
        cids = []
        for v in df_in[col].dropna().head(TOP_N):
            try:
                cids.append(int(float(v)))
            except (ValueError, TypeError):
                continue
        q, v, s = build_degree_pair(cids, dock_results)
        qsar_by_deg[col] = q
        vina_by_deg[col] = v
        score_by_deg[col] = s

    # Wide table: 9 columns, TOP_N rows
    out = {}
    for i, col in enumerate(DEGREE_COLUMNS, 1):
        out[f"Deg{i}_QSAR"]       = qsar_by_deg[col]
        out[f"Deg{i}_Vina"]       = vina_by_deg[col]
        out[f"Deg{i}_Vina_Score"] = score_by_deg[col]
    df_wide = pd.DataFrame(out)
    df_wide.insert(0, "rank", range(1, TOP_N + 1))

    # Rank-changes long table
    df_changes = build_rank_changes(qsar_by_deg, vina_by_deg)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / f"{target}_qsar_vs_vina.xlsx"
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        df_wide.to_excel(writer, sheet_name="qsar_vs_vina", index=False)
        df_changes.to_excel(writer, sheet_name="rank_changes", index=False)
    print(f"  wrote: {out_path}")

    # Console preview
    print(f"  preview:")
    with pd.option_context("display.max_columns", None,
                           "display.width", 200):
        print(df_wide.to_string(index=False))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for target, xlsx in TARGETS.items():
        process_target(target, xlsx)
    print(f"\nDone. Outputs in: {OUT_DIR}")


if __name__ == "__main__":
    main()