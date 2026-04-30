"""
03_dock_generated.py — Dock DiffSBDD-generated ligands to get binding energies.

For each selected pocket:
  - Read <T>_pocket<N>_ligands_clean.sdf (from Step 2)
  - Convert each ligand to PDBQT
  - Dock against receptor with a pocket-centered box
  - Extract Vina score → write <T>_pocket<N>_generated_docking.csv

Runtime: ~10-20 sec per ligand with 8 parallel workers.
  50 ligands × 3 pockets = ~5-10 min total.
"""

import os
import sys
import subprocess
import multiprocessing as mp
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from config import (
    TARGET_NAME, N_WORKERS, BOX_SIZE, EXHAUSTIVENESS, N_POSES, DOCKING_TIMEOUT,
    RECEPTOR_PDBQT, POCKETS_CSV,
    pocket_dir, pocket_file,
)

RDLogger.DisableLog("rdApp.*")


def sdf_to_pdbqt_single(mol, out_path):
    """Convert one RDKit mol (with 3D coords from DiffSBDD) to PDBQT."""
    # Ensure 3D coords exist and add hydrogens
    mol = Chem.AddHs(mol, addCoords=True)
    # If no 3D conformer, generate one
    if mol.GetNumConformers() == 0:
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
    """Read Vina output and return the best (first) binding energy."""
    if not os.path.exists(pose_file):
        return None
    with open(pose_file) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                return float(line.split()[3])
    return None


def dock_one(args):
    """Worker: dock a single ligand.
    args = (ligand_id, pdbqt_path, pose_path, cx, cy, cz)
    """
    lig_id, lig_pdbqt, pose_out, cx, cy, cz = args
    if os.path.exists(pose_out):
        return (lig_id, extract_best_score(pose_out), "ok")
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
        return (lig_id, extract_best_score(pose_out), "ok")
    except subprocess.TimeoutExpired:
        return (lig_id, None, "timeout")
    except Exception as e:
        return (lig_id, None, f"err:{type(e).__name__}")


def process_pocket(n, cx, cy, cz):
    """Dock all valid ligands for one pocket."""
    pdir = pocket_dir(n)
    clean_sdf = pocket_file(n, "ligands_clean.sdf")
    out_csv   = pocket_file(n, "generated_docking.csv")

    if os.path.exists(out_csv) and os.path.getsize(out_csv) > 0:
        print(f"[skip]  pocket{n}: generated docking already done")
        return pd.read_csv(out_csv)

    if not os.path.exists(clean_sdf):
        print(f"[skip]  pocket{n}: no clean SDF (Step 2 may have failed)")
        return None

    # Prep directories
    lig_dir  = os.path.join(pdir, "generated_ligands_pdbqt")
    pose_dir = os.path.join(pdir, "generated_poses")
    os.makedirs(lig_dir, exist_ok=True)
    os.makedirs(pose_dir, exist_ok=True)

    # Convert each SDF entry → PDBQT + collect SMILES
    print(f"  preparing PDBQTs for pocket{n}...")
    suppl = Chem.SDMolSupplier(clean_sdf, sanitize=True, removeHs=False)
    tasks = []
    smiles_map = {}
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        lig_id = f"gen_{i:03d}"
        lig_pdbqt = os.path.join(lig_dir, f"{lig_id}.pdbqt")
        pose_out  = os.path.join(pose_dir, f"{lig_id}_out.pdbqt")
        if not os.path.exists(lig_pdbqt):
            if not sdf_to_pdbqt_single(mol, lig_pdbqt):
                continue
        smiles_map[lig_id] = Chem.MolToSmiles(mol, True)
        tasks.append((lig_id, lig_pdbqt, pose_out, cx, cy, cz))

    if not tasks:
        print(f"  no valid ligands to dock for pocket{n}")
        return None

    # Dock in parallel
    print(f"  docking {len(tasks)} ligands for pocket{n} with {N_WORKERS} workers...")
    results = []
    with mp.Pool(N_WORKERS) as pool:
        for i, r in enumerate(pool.imap_unordered(dock_one, tasks), 1):
            results.append(r)
            if i % 10 == 0 or i == len(tasks):
                print(f"    {i}/{len(tasks)}")

    # Build output CSV
    df = pd.DataFrame(results, columns=["ligand_id", "ba", "status"])
    df["smiles"] = df["ligand_id"].map(smiles_map)
    df["x"] = cx
    df["y"] = cy
    df["z"] = cz
    df["pocket_num"] = n
    df = df[["pocket_num", "ligand_id", "smiles", "ba", "x", "y", "z", "status"]]
    df = df.sort_values("ba", ascending=True, na_position="last").reset_index(drop=True)
    df.to_csv(out_csv, index=False)

    ok = df["ba"].notna().sum()
    best = df["ba"].min() if ok else None
    print(f"  pocket{n}: docked {ok}/{len(tasks)}, best ba={best}")
    return df


def main():
    print(f"\n### Step 3: Dock generated ligands for {TARGET_NAME} ###\n")

    if not os.path.exists(POCKETS_CSV):
        sys.exit(f"ERROR: {POCKETS_CSV} missing. Run 01_prepare.py first.")
    if not os.path.exists(RECEPTOR_PDBQT):
        sys.exit(f"ERROR: receptor PDBQT missing: {RECEPTOR_PDBQT}")

    df_pockets = pd.read_csv(POCKETS_CSV)
    selected = df_pockets[df_pockets["selected"] == True].reset_index(drop=True)
    print(f"Processing {len(selected)} pocket(s)\n")

    summary = []
    for _, row in selected.iterrows():
        n = int(row["pocket_num"])
        print(f"--- pocket{n} ---")
        df = process_pocket(n, row["center_x"], row["center_y"], row["center_z"])
        if df is not None:
            valid = df["ba"].notna()
            summary.append((n, valid.sum(), df.loc[valid, "ba"].min() if valid.any() else None))
        print()

    print("=" * 60)
    print("Generated-ligand docking summary:")
    for n, count, best in summary:
        best_str = f"{best:.2f}" if best is not None else "n/a"
        print(f"  pocket{n}: {count} ligands docked, best ba={best_str} kcal/mol")
    print()
    print("Next: python 04_morgan_rank.py")


if __name__ == "__main__":
    main()