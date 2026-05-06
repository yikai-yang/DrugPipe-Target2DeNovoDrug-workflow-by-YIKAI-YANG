"""
07_analyze.py — Per-target and cross-target analysis.

For each target:
  - Loads all per-pocket final.csv files
  - Applies filters (drug-like MW, valid SMILES, successful docking)
  - Computes per-pocket statistics (best score, median, spread)
  - Identifies top candidates with their chemotypes
  - Flags artifacts (tiny molecules, known promiscuous binders)

Also writes a cross-target comparison CSV.

Run:
  python 07_analyze.py                    # analyzes all targets in runs/
  python 07_analyze.py IGFBP7 BGN         # only specific targets
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, Lipinski

RDLogger.DisableLog("rdApp.*")

# ==============================================================
# CONFIG
# ==============================================================
BASE_DIR   = Path(os.environ.get("DRUGPIPE_BASE_DIR", Path(__file__).resolve().parent.parent))
RUNS_DIR   = BASE_DIR / "runs"
ANALYSIS_DIR = BASE_DIR / "analysis"
ANALYSIS_DIR.mkdir(exist_ok=True)

# Drug-likeness filters — drugs failing these get flagged, not removed
MIN_MW          = 150     # below this, likely metabolite/fragment not a drug
MAX_MW          = 600     # Lipinski's upper bound; very large molecules dock unreliably
MAX_ROTATABLE   = 10      # flexible molecules give unreliable Vina scores

# Binding strength thresholds (kcal/mol)
STRONG_BA       = -8.0    # very likely hit
MEDIUM_BA       = -7.0    # plausible hit
WEAK_BA         = -6.0    # borderline

# Known promiscuous / "pan-assay interference" scaffolds worth flagging
# (not a complete list — just common ones that show up in these screens)
PAINS_SMARTS = [
    ("catechol",     "c1cc(O)c(O)cc1"),
    ("rhodanine",    "S=C1SC(=O)NC1"),
    ("quinone",      "O=C1C=CC(=O)C=C1"),
]


def pocket_num_from_dir(p):
    """Parse '5' from 'pocket5'."""
    name = p.name
    if not name.startswith("pocket"):
        return None
    try:
        return int(name.replace("pocket", ""))
    except ValueError:
        return None


def annotate_molecule(smi):
    """Compute drug-likeness descriptors. Returns dict of flags/values."""
    if not isinstance(smi, str):
        return {"mw": None, "logp": None, "rotatable": None,
                "hbd": None, "hba": None, "too_small": True,
                "too_large": False, "too_flexible": False,
                "lipinski_ok": False, "pains_flags": ""}
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return {"mw": None, "logp": None, "rotatable": None,
                "hbd": None, "hba": None, "too_small": True,
                "too_large": False, "too_flexible": False,
                "lipinski_ok": False, "pains_flags": "parse_fail"}

    mw     = Descriptors.MolWt(mol)
    logp   = Descriptors.MolLogP(mol)
    rot    = Lipinski.NumRotatableBonds(mol)
    hbd    = Lipinski.NumHDonors(mol)
    hba    = Lipinski.NumHAcceptors(mol)

    too_small    = mw < MIN_MW
    too_large    = mw > MAX_MW
    too_flexible = rot > MAX_ROTATABLE
    # Lipinski Ro5: MW<=500, logP<=5, HBD<=5, HBA<=10
    lipinski_ok  = mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10

    pains_flags = []
    for name, smarts in PAINS_SMARTS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None and mol.HasSubstructMatch(patt):
            pains_flags.append(name)

    return {"mw": round(mw, 1), "logp": round(logp, 2), "rotatable": rot,
            "hbd": hbd, "hba": hba, "too_small": too_small,
            "too_large": too_large, "too_flexible": too_flexible,
            "lipinski_ok": lipinski_ok, "pains_flags": ",".join(pains_flags)}


def chemotype_cluster(smiles_list, threshold=0.6):
    """Cheap chemotype clustering via Tanimoto. Returns cluster assignments.
    Drugs with pairwise Tanimoto >= threshold are clustered together.
    """
    fps = []
    for s in smiles_list:
        mol = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if mol is None:
            fps.append(None)
            continue
        fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))

    clusters = [-1] * len(fps)
    current_cluster = 0
    for i, fp_i in enumerate(fps):
        if fp_i is None or clusters[i] != -1:
            continue
        clusters[i] = current_cluster
        for j in range(i + 1, len(fps)):
            if fps[j] is None or clusters[j] != -1:
                continue
            sim = DataStructs.TanimotoSimilarity(fp_i, fps[j])
            if sim >= threshold:
                clusters[j] = current_cluster
        current_cluster += 1

    return clusters


def load_target_results(target):
    """Read all pocket final.csv files and combine into one dataframe."""
    target_dir = RUNS_DIR / target
    if not target_dir.exists():
        return None

    pockets_summary = target_dir / f"{target}_pockets_summary.csv"
    druggability = {}
    if pockets_summary.exists():
        pdf = pd.read_csv(pockets_summary)
        druggability = dict(zip(pdf["pocket_num"].astype(int), pdf["druggability"]))

    parts = []
    for pdir in sorted(target_dir.glob("pocket*")):
        n = pocket_num_from_dir(pdir)
        if n is None:
            continue
        final_csv = pdir / f"{target}_pocket{n}_final.csv"
        if not final_csv.exists():
            continue
        df = pd.read_csv(final_csv)
        df["pocket_num"] = n
        df["pocket_druggability"] = druggability.get(n, np.nan)
        parts.append(df)

    if not parts:
        return None
    return pd.concat(parts, ignore_index=True)


def analyze_target(target):
    print(f"\n{'='*72}")
    print(f"Analyzing: {target}")
    print(f"{'='*72}")

    df = load_target_results(target)
    if df is None or len(df) == 0:
        print(f"  no results for {target}, skipping")
        return None

    # Keep only successfully docked rows
    n_total = len(df)
    df = df[df["vina_score"].notna()].copy()
    n_docked = len(df)
    print(f"  {n_docked}/{n_total} rows have valid vina_score")

    # Annotate drug-likeness
    print("  annotating drug-likeness...")
    annotations = df["smiles"].apply(annotate_molecule).apply(pd.Series)
    df = pd.concat([df.reset_index(drop=True), annotations.reset_index(drop=True)], axis=1)

    # Per-pocket statistics
    print("\n  Per-pocket summary:")
    print(f"  {'pocket':<8} {'drugg':<8} {'best':<8} {'median':<8} {'count':<8}")
    per_pocket = []
    for n, grp in df.groupby("pocket_num"):
        scores = grp["vina_score"].dropna()
        if len(scores) == 0:
            continue
        drugg = grp["pocket_druggability"].iloc[0]
        row = {
            "target": target,
            "pocket_num": n,
            "druggability": round(drugg, 3) if pd.notna(drugg) else None,
            "n_docked": len(scores),
            "best_ba": round(scores.min(), 2),
            "median_ba": round(scores.median(), 2),
            "frac_strong": round((scores <= STRONG_BA).mean(), 3),
            "frac_medium": round((scores <= MEDIUM_BA).mean(), 3),
        }
        per_pocket.append(row)
        print(f"  pocket{n:<2d} {drugg:<8.3f} {row['best_ba']:<8} "
              f"{row['median_ba']:<8} {row['n_docked']:<8}")

    # Best pocket per drug
    best_per_drug = df.loc[df.groupby("Drug id")["vina_score"].idxmin()].copy()
    best_per_drug = best_per_drug.sort_values("vina_score").reset_index(drop=True)
    best_per_drug["overall_rank"] = range(1, len(best_per_drug) + 1)

    # Chemotype clustering on top 30
    top30 = best_per_drug.head(30).copy()
    if len(top30) >= 2:
        top30["chemotype_cluster"] = chemotype_cluster(top30["smiles"].tolist(), threshold=0.5)

    # Build "clean" shortlist: drug-like + non-PAINS
    clean = best_per_drug[
        (~best_per_drug["too_small"]) &
        (~best_per_drug["too_large"]) &
        (~best_per_drug["too_flexible"]) &
        (best_per_drug["pains_flags"] == "") &
        (best_per_drug["vina_score"] <= WEAK_BA)
    ].copy()
    clean["clean_rank"] = range(1, len(clean) + 1)

    # Print top 10 of clean list
    print(f"\n  Top 10 drug-like candidates (MW {MIN_MW}-{MAX_MW}, rotbonds≤{MAX_ROTATABLE}, non-PAINS):")
    show_cols = ["clean_rank", "Drug name", "pocket_num", "vina_score",
                 "mw", "tanimoto_score"]
    show_cols = [c for c in show_cols if c in clean.columns]
    if len(clean) == 0:
        print("    no candidates passed filters — likely all hits were tiny or non-drug-like")
    else:
        print(clean[show_cols].head(10).to_string(index=False))

    # Save outputs
    target_analysis_dir = ANALYSIS_DIR / target
    target_analysis_dir.mkdir(exist_ok=True)
    best_per_drug.to_csv(target_analysis_dir / f"{target}_all_ranked.csv", index=False)
    clean.to_csv(target_analysis_dir / f"{target}_clean_shortlist.csv", index=False)
    top30.to_csv(target_analysis_dir / f"{target}_top30_with_chemotypes.csv", index=False)
    if per_pocket:
        pd.DataFrame(per_pocket).to_csv(
            target_analysis_dir / f"{target}_pocket_stats.csv", index=False)

    # Return summary for cross-target table
    best_row = best_per_drug.iloc[0] if len(best_per_drug) else None
    return {
        "target": target,
        "n_pockets": df["pocket_num"].nunique(),
        "n_drugs_docked": df["Drug id"].nunique(),
        "best_vina": round(best_per_drug["vina_score"].min(), 2) if len(best_per_drug) else None,
        "best_drug": best_row["Drug name"] if best_row is not None else None,
        "best_pocket": int(best_row["pocket_num"]) if best_row is not None else None,
        "n_strong_ba": int((best_per_drug["vina_score"] <= STRONG_BA).sum()),
        "n_drug_like": len(clean),
        "clean_best_vina": round(clean["vina_score"].min(), 2) if len(clean) else None,
        "clean_best_drug": clean.iloc[0]["Drug name"] if len(clean) else None,
    }


def main():
    # Which targets to process
    if len(sys.argv) > 1:
        targets = sys.argv[1:]
    else:
        targets = sorted(d.name for d in RUNS_DIR.iterdir()
                         if d.is_dir() and (d / f"{d.name}_pockets_summary.csv").exists())

    if not targets:
        sys.exit("No targets found in runs/")

    print(f"Analyzing {len(targets)} target(s): {', '.join(targets)}")

    summaries = []
    for t in targets:
        s = analyze_target(t)
        if s is not None:
            summaries.append(s)

    # Cross-target comparison
    if summaries:
        summary_df = pd.DataFrame(summaries)
        summary_df = summary_df.sort_values("best_vina", na_position="last").reset_index(drop=True)
        out = ANALYSIS_DIR / "cross_target_summary.csv"
        summary_df.to_csv(out, index=False)

        print(f"\n{'='*72}")
        print("CROSS-TARGET SUMMARY")
        print(f"{'='*72}")
        cols = ["target", "best_vina", "best_drug", "best_pocket",
                "n_strong_ba", "n_drug_like", "clean_best_vina", "clean_best_drug"]
        cols = [c for c in cols if c in summary_df.columns]
        print(summary_df[cols].to_string(index=False))
        print(f"\nWritten to: {out}")


if __name__ == "__main__":
    main()