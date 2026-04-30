"""
05_dock_drugbank.py — Dock top-N DrugBank candidates per pocket.

For each selected pocket:
  - Read <T>_pocket<N>_morgan.csv (from Step 4)
  - Take top TOP_N_DOCK drugs
  - Convert SMILES → 3D PDBQT
  - Dock against receptor at the pocket center
  - Extract Vina score → merge into <T>_pocket<N>_final.csv

Reuses prepared ligand PDBQTs across pockets (DrugBank drug is the same
molecule regardless of which pocket it's being docked at).

Runtime: TOP_N_DOCK ligands × 15s / N_WORKERS. For 100 × 3 pockets with 8
workers: ~20-30 min.
"""

import os
import sys
import subprocess
import multiprocessing as mp
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from config import (
    TARGET_NAME, N_WORKERS, BOX_SIZE, EXHAUSTIVENESS, N_POSES,
    DOCKING_TIMEOUT, TOP_N_DOCK,
    RECEPTOR_PDBQT, POCKETS_CSV, TARGET_DIR,
    pocket_dir, pocket_file,
)

RDLogger.DisableLog("rdApp.*")

# Shared across pockets — one PDBQT per DrugBank drug, reused
DB_LIGAND_CACHE = os.path.join(TARGET_DIR, "drugbank_ligands_pdbqt")
os.makedirs(DB_LIGAND_CACHE, exist_ok=True)


def smiles_to_pdbqt(smi, out_path):
    """SMILES → 3D → PDBQT via RDKit + obabel. Returns True on success."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        return False
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass
    pdb_tmp = out_path.replace(".pdbqt", ".pdb")
    Chem.MolToPDBFile(mol, pdb_tmp)
    try:
        subprocess.check_call(
            ["obabel", pdb_tmp, "-O", out_path, "--partialcharge", "gasteiger"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        os.remove(pdb_tmp)
        return os.path.exists(out_path) and os.path.getsize(out_path) > 0
    except Exception:
        return False


def extract_best_score(pose_file):
    """Read best Vina score from pose .pdbqt output."""
    if not os.path.exists(pose_file):
        return None
    with open(pose_file) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                return float(line.split()[3])
    return None


def dock_one(args):
    """Parallel worker: dock one ligand at one pocket."""
    dbid, lig_pdbqt, pose_out, cx, cy, cz = args
    if os.path.exists(pose_out):
        return (dbid, extract_best_score(pose_out), "ok")
    cmd = [
        "vina",
        "--receptor", RECEPTOR_PDBQT,
        "--ligand",   lig_pdbqt,
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x",   str(BOX_SIZE),
        "--size_y",   str(BOX_SIZE),
        "--size_z",   str(BOX_SIZE),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(N_POSES),
        "--out", pose_out,
        "--cpu", "1",
    ]
    seed = os.environ.get("DRUGPIPE_VINA_SEED")
    if seed:
        cmd += ["--seed", seed]
    log = pose_out.replace("_out.pdbqt", ".log")
    try:
        with open(log, "w") as lf:
            subprocess.check_call(cmd, stdout=lf, stderr=subprocess.STDOUT,
                                  timeout=DOCKING_TIMEOUT)
        return (dbid, extract_best_score(pose_out), "ok")
    except subprocess.TimeoutExpired:
        return (dbid, None, "timeout")
    except Exception as e:
        return (dbid, None, f"err:{type(e).__name__}")


def prepare_drugbank_pdbqts(top_drugs):
    """Prepare PDBQT for each top drug (reused across pockets)."""
    prepared = {}
    for _, row in top_drugs.iterrows():
        dbid = row["Drug id"]
        smi  = row["smiles"]
        lig_pdbqt = os.path.join(DB_LIGAND_CACHE, f"{dbid}.pdbqt")
        if not os.path.exists(lig_pdbqt):
            if not smiles_to_pdbqt(smi, lig_pdbqt):
                continue
        prepared[dbid] = lig_pdbqt
    return prepared


def process_pocket(n, cx, cy, cz):
    """Dock top-N DrugBank candidates at pocket N."""
    morgan_csv = pocket_file(n, "morgan.csv")
    final_csv  = pocket_file(n, "final.csv")

    if os.path.exists(final_csv) and os.path.getsize(final_csv) > 0:
        print(f"[skip]  pocket{n}: final CSV already exists")
        return pd.read_csv(final_csv)

    if not os.path.exists(morgan_csv):
        print(f"[skip]  pocket{n}: morgan.csv missing (Step 4 failed?)")
        return None

    df = pd.read_csv(morgan_csv)
    top = df.head(TOP_N_DOCK).copy()
    print(f"  pocket{n}: taking top {len(top)} DrugBank candidates")

    # Prepare PDBQTs (cached in shared folder)
    print(f"  preparing ligand PDBQTs (reusing cache)...")
    prepared = prepare_drugbank_pdbqts(top)
    print(f"  {len(prepared)}/{len(top)} ligands ready")

    # Per-pocket pose directory
    pose_dir = os.path.join(pocket_dir(n), "drugbank_poses")
    os.makedirs(pose_dir, exist_ok=True)

    # Build docking tasks
    tasks = []
    for _, row in top.iterrows():
        dbid = row["Drug id"]
        if dbid not in prepared:
            continue
        pose_out = os.path.join(pose_dir, f"{dbid}_out.pdbqt")
        tasks.append((dbid, prepared[dbid], pose_out, cx, cy, cz))

    # Run parallel
    print(f"  docking {len(tasks)} ligands with {N_WORKERS} workers...")
    results = []
    with mp.Pool(N_WORKERS) as pool:
        for i, r in enumerate(pool.imap_unordered(dock_one, tasks), 1):
            results.append(r)
            if i % 10 == 0 or i == len(tasks):
                print(f"    {i}/{len(tasks)}")

    # Merge vina scores back into the ranked table
    res_df = pd.DataFrame(results, columns=["Drug id", "vina_score", "dock_status"])
    merged = top.merge(res_df, on="Drug id", how="left")

    # Re-rank by actual Vina score (more negative = better)
    merged = merged.sort_values("vina_score", ascending=True, na_position="last").reset_index(drop=True)
    merged["final_rank"] = range(1, len(merged) + 1)

    # Save
    cols = ["final_rank", "pocket_num", "Drug name", "Drug id", "smiles",
            "tanimoto_score", "inherited_ba", "vina_score", "combined_score",
            "dock_status"]
    cols = [c for c in cols if c in merged.columns]
    merged[cols].to_csv(final_csv, index=False)
    print(f"  wrote {final_csv}")

    # Quick peek
    ok = merged["vina_score"].notna().sum()
    best = merged["vina_score"].min() if ok else None
    best_str = f"{best:.2f}" if best is not None else "n/a"
    print(f"\n  pocket{n}: {ok}/{len(tasks)} docked, best vina={best_str} kcal/mol")
    peek_cols = [c for c in ["final_rank", "Drug name", "tanimoto_score", "vina_score"]
                 if c in merged.columns]
    print(f"  Top 5 for pocket{n}:")
    print(merged[peek_cols].head(5).to_string(index=False))
    return merged


def main():
    print(f"\n### Step 5: Dock top-{TOP_N_DOCK} DrugBank per pocket for {TARGET_NAME} ###\n")

    if not os.path.exists(RECEPTOR_PDBQT):
        sys.exit(f"ERROR: receptor PDBQT missing: {RECEPTOR_PDBQT}")
    if not os.path.exists(POCKETS_CSV):
        sys.exit(f"ERROR: {POCKETS_CSV} missing. Run 01_prepare.py first.")

    df_pockets = pd.read_csv(POCKETS_CSV)
    selected = df_pockets[df_pockets["selected"] == True].reset_index(drop=True)
    print(f"Processing {len(selected)} pocket(s)")
    print(f"DrugBank PDBQT cache: {DB_LIGAND_CACHE}\n")

    summary = []
    for _, row in selected.iterrows():
        n = int(row["pocket_num"])
        print(f"--- pocket{n} ---")
        df = process_pocket(n, row["center_x"], row["center_y"], row["center_z"])
        if df is not None:
            ok = df["vina_score"].notna().sum()
            best = df["vina_score"].min() if ok else None
            summary.append((n, ok, best))
        print()

    print("=" * 60)
    print("DrugBank docking summary:")
    for n, ok, best in summary:
        best_str = f"{best:.2f}" if best is not None else "n/a"
        print(f"  pocket{n}: {ok} docked, best vina={best_str} kcal/mol")
    print()
    print("Next: python 06_aggregate.py")


if __name__ == "__main__":
    main()