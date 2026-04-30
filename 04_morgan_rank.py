"""
04_morgan_rank.py — Rank DrugBank drugs per pocket using Morgan FP + inherited ba.

For each pocket:
  - Load generated_docking.csv (ligands + ba from Step 3)
  - Compute Morgan FPs for all generated ligands and all DrugBank drugs
  - For each DrugBank drug, find best-matching generated ligand → inherit its ba
  - Rank drugs by multi-objective score (similarity + ba)

Outputs per pocket:
  runs/<T>/pocket<N>/<T>_pocket<N>_morgan.csv

Runtime: ~2-3 min per pocket (dominated by fingerprinting all 9,716 DrugBank entries).
"""

import os
import sys
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

from config import (
    TARGET_NAME, MORGAN_RADIUS, MORGAN_NBITS,
    DRUGBANK_CSV, POCKETS_CSV,
    pocket_file,
)

RDLogger.DisableLog("rdApp.*")


def smiles_to_fp(smi):
    """Convert SMILES string to Morgan fingerprint, or None if invalid."""
    if not isinstance(smi, str) or not smi.strip():
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=MORGAN_NBITS)


def rank_sum_score(df):
    """Multi-objective score combining structural similarity and binding energy.

    Following the DrugPipe paper: rank each criterion independently, then combine
    with equal weighting. Higher score = better candidate.

    f(x) = (tanimoto_norm) - (ba_norm)     [ba is negative, lower is better]

    We normalize each metric 0-1 (higher=better), then sum.
    """
    s = df["tanimoto_score"].to_numpy()
    b = df["inherited_ba"].to_numpy()

    # Similarity: already 0-1, higher=better
    # Fill NaN with 0 (no info = lowest)
    s_norm = np.nan_to_num(s, nan=0.0)

    # Binding energy: invert so higher=better, then min-max normalize
    # Typical ba range is -10 to -3, so -b gives 3 to 10
    if np.all(np.isnan(b)):
        b_norm = np.zeros_like(b)
    else:
        b_inv = -np.nan_to_num(b, nan=0.0)  # flip sign: higher = more negative = better
        b_min, b_max = np.nanmin(b_inv), np.nanmax(b_inv)
        if b_max > b_min:
            b_norm = (b_inv - b_min) / (b_max - b_min)
        else:
            b_norm = np.zeros_like(b_inv)

    return 0.5 * s_norm + 0.5 * b_norm


def process_pocket(n):
    """Morgan rank DrugBank against pocket N's generated ligands."""
    gen_csv  = pocket_file(n, "generated_docking.csv")
    out_csv  = pocket_file(n, "morgan.csv")

    if os.path.exists(out_csv) and os.path.getsize(out_csv) > 0:
        print(f"[skip]  pocket{n}: morgan ranking already done")
        return pd.read_csv(out_csv)

    if not os.path.exists(gen_csv):
        print(f"[skip]  pocket{n}: no generated_docking.csv (Step 3 failed?)")
        return None

    gen = pd.read_csv(gen_csv)
    gen = gen[gen["ba"].notna() & gen["smiles"].notna()].reset_index(drop=True)
    if len(gen) == 0:
        print(f"[skip]  pocket{n}: no valid docked ligands to compare against")
        return None

    print(f"  pocket{n}: {len(gen)} generated ligands with docking scores")

    # Build fingerprints for generated ligands
    gen_fps, gen_ba, gen_smi = [], [], []
    for _, row in gen.iterrows():
        fp = smiles_to_fp(row["smiles"])
        if fp is None:
            continue
        gen_fps.append(fp)
        gen_ba.append(row["ba"])
        gen_smi.append(row["smiles"])

    if not gen_fps:
        print(f"  pocket{n}: no fingerprintable generated ligands")
        return None

    gen_ba = np.array(gen_ba)

    # Load DrugBank
    drugbank = pd.read_csv(DRUGBANK_CSV)
    print(f"  comparing against {len(drugbank)} DrugBank drugs...")

    # Fingerprint each DrugBank drug and find its best match
    best_sim   = np.zeros(len(drugbank))
    best_idx   = np.zeros(len(drugbank), dtype=int)
    valid_mask = np.ones(len(drugbank), dtype=bool)

    for i, smi in enumerate(drugbank["smiles"]):
        fp = smiles_to_fp(smi)
        if fp is None:
            valid_mask[i] = False
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp, gen_fps)
        best_sim[i] = max(sims)
        best_idx[i] = int(np.argmax(sims))

        if (i + 1) % 2000 == 0:
            print(f"    {i+1}/{len(drugbank)}")

    # Build output
    drugbank = drugbank[valid_mask].copy().reset_index(drop=True)
    best_sim  = best_sim[valid_mask]
    best_idx  = best_idx[valid_mask]

    drugbank["tanimoto_score"]      = best_sim
    drugbank["inherited_ba"]        = gen_ba[best_idx]
    drugbank["matched_ligand_idx"]  = best_idx
    drugbank["matched_ligand_smi"]  = [gen_smi[i] for i in best_idx]
    drugbank["pocket_num"]          = n

    # Multi-objective score + rank
    drugbank["combined_score"] = rank_sum_score(drugbank)
    drugbank = drugbank.sort_values("combined_score", ascending=False).reset_index(drop=True)
    drugbank["rank_in_pocket"] = range(1, len(drugbank) + 1)

    # Save
    cols = ["rank_in_pocket", "pocket_num", "Drug name", "Drug id", "smiles",
            "tanimoto_score", "inherited_ba", "combined_score",
            "matched_ligand_idx", "matched_ligand_smi"]
    cols = [c for c in cols if c in drugbank.columns]
    drugbank[cols].to_csv(out_csv, index=False)
    print(f"  wrote {out_csv}")

    # Quick peek at top 5
    print(f"\n  Top 5 for pocket{n}:")
    peek = drugbank[["rank_in_pocket", "Drug name", "tanimoto_score", "inherited_ba"]].head(5)
    print(peek.to_string(index=False))
    return drugbank


def main():
    print(f"\n### Step 4: Morgan FP ranking for {TARGET_NAME} ###\n")

    if not os.path.exists(DRUGBANK_CSV):
        sys.exit(f"ERROR: DrugBank CSV missing: {DRUGBANK_CSV}")
    if not os.path.exists(POCKETS_CSV):
        sys.exit(f"ERROR: {POCKETS_CSV} missing. Run 01_prepare.py first.")

    df_pockets = pd.read_csv(POCKETS_CSV)
    selected = df_pockets[df_pockets["selected"] == True].reset_index(drop=True)
    print(f"Processing {len(selected)} pocket(s)\n")

    for _, row in selected.iterrows():
        n = int(row["pocket_num"])
        print(f"--- pocket{n} ---")
        process_pocket(n)
        print()

    print("=" * 60)
    print("Morgan FP ranking done for all pockets.")
    print("Next: python 05_dock_drugbank.py")


if __name__ == "__main__":
    main()