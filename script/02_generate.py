"""
02_generate.py — Run DiffSBDD to generate candidate ligands per pocket.

Reads from config.py:
    TARGET_NAME, N_LIGANDS, DIFFSBDD_DIR, DIFFSBDD_CKPT, DIFFSBDD_TIMESTEPS
    POCKETS_CSV, INPUT_PDB
    pocket_dir(n), pocket_file(n, suffix)

Inputs:
    runs/<T>/input/<T>.pdb
    runs/<T>/<T>_pockets_summary.csv       (from Step 1)

Outputs (per selected pocket N):
    runs/<T>/pocket<N>/<T>_pocket<N>_marker.sdf        (pocket center marker)
    runs/<T>/pocket<N>/<T>_pocket<N>_ligands.sdf       (raw DiffSBDD output)
    runs/<T>/pocket<N>/<T>_pocket<N>_ligands_clean.sdf (sanitized with RDKit)

Runtime: ~10 sec per ligand on CPU → ~15 min per pocket at N_LIGANDS=100.
"""

import os
import sys
import subprocess
import time
import pandas as pd
from rdkit import Chem, RDLogger

from config import (
    TARGET_NAME, N_LIGANDS, DIFFSBDD_TIMESTEPS,
    DIFFSBDD_DIR, DIFFSBDD_CKPT,
    INPUT_PDB, POCKETS_CSV,
    pocket_dir, pocket_file,
)

# Silence RDKit's loud warnings about unparseable molecules
RDLogger.DisableLog("rdApp.*")


def write_pocket_marker(x, y, z, out_path):
    """Write a minimal valid SDF with a single carbon at the pocket center.
    DiffSBDD uses the position of this 'fake ligand' to decide where to
    place generated atoms.
    """
    sdf = (
        "fake_ligand\n"
        "     RDKit          3D\n"
        "\n"
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        f"{x:10.4f}{y:10.4f}{z:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n"
        "$$$$\n"
    )
    with open(out_path, "w") as f:
        f.write(sdf)

    # Validate the marker parses
    if Chem.SDMolSupplier(out_path)[0] is None:
        sys.exit(f"ERROR: wrote invalid marker SDF: {out_path}")


def run_diffsbdd(pdb_in, marker_sdf, out_sdf):
    """Invoke DiffSBDD's generate_ligands.py for one pocket."""
    ckpt_full = os.path.join(DIFFSBDD_DIR, DIFFSBDD_CKPT)
    if not os.path.exists(ckpt_full):
        sys.exit(f"ERROR: checkpoint missing: {ckpt_full}")

    script_full = os.path.join(DIFFSBDD_DIR, "generate_ligands.py")
    if not os.path.exists(script_full):
        sys.exit(f"ERROR: generate_ligands.py missing: {script_full}")

    cmd = [
        sys.executable,
        "-u",  # unbuffered child stdout/stderr
        "generate_ligands.py",
        DIFFSBDD_CKPT,
        "--pdbfile", pdb_in,
        "--outfile", out_sdf,
        "--ref_ligand", marker_sdf,
        "--n_samples", str(N_LIGANDS),
        "--timesteps", str(DIFFSBDD_TIMESTEPS),
        "--sanitize",
    ]

    print(f"    → running DiffSBDD")
    print("      cwd:", DIFFSBDD_DIR)
    print("      cmd:", " ".join(cmd))
    print("      pdb:", pdb_in)
    print("      marker:", marker_sdf)
    print("      out:", out_sdf)

    t0 = time.time()
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    try:
        proc = subprocess.Popen(
            cmd,
            cwd=DIFFSBDD_DIR,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        last_output = time.time()

        while True:
            line = proc.stdout.readline()
            if line:
                last_output = time.time()
                print(f"      | {line.rstrip()}")
            elif proc.poll() is not None:
                break
            else:
                if time.time() - last_output > 30:
                    print("      | ...still running, no new output in 30s...")
                    last_output = time.time()
                time.sleep(1)

        rc = proc.wait()
        if rc != 0:
            sys.exit(f"ERROR: DiffSBDD exited with code {rc}")

    except FileNotFoundError:
        sys.exit(f"ERROR: couldn't launch DiffSBDD from {DIFFSBDD_DIR}")

    elapsed = time.time() - t0
    print(f"    → done in {elapsed/60:.1f} min")


def sanitize_sdf(in_sdf, out_sdf):
    """Read DiffSBDD output, keep only RDKit-parseable molecules, write clean SDF.
    Returns (total_in, total_kept).
    """
    if not os.path.exists(in_sdf):
        sys.exit(f"ERROR: DiffSBDD produced no output: {in_sdf}")

    suppl = Chem.SDMolSupplier(in_sdf, sanitize=True, removeHs=False)
    writer = Chem.SDWriter(out_sdf)

    total, kept = 0, 0
    for mol in suppl:
        total += 1
        if mol is None:
            continue
        try:
            # Canonicalize via SMILES round-trip to catch any remaining issues
            smi = Chem.MolToSmiles(mol, True)
            if Chem.MolFromSmiles(smi) is None:
                continue
            writer.write(mol)
            kept += 1
        except Exception:
            continue

    writer.close()
    return total, kept


def process_pocket(pocket_num, cx, cy, cz):
    """Full per-pocket generation flow."""
    pdir = pocket_dir(pocket_num)
    os.makedirs(pdir, exist_ok=True)

    marker_path = pocket_file(pocket_num, "marker.sdf")
    raw_path    = pocket_file(pocket_num, "ligands.sdf")
    clean_path  = pocket_file(pocket_num, "ligands_clean.sdf")

    if os.path.exists(clean_path) and os.path.getsize(clean_path) > 0:
        # Count entries in existing file
        existing = sum(1 for _ in Chem.SDMolSupplier(clean_path) if _ is not None)
        print(f"[skip]  pocket{pocket_num}: {existing} ligands already generated")
        return existing

    # Write marker
    if not os.path.exists(marker_path):
        write_pocket_marker(cx, cy, cz, marker_path)
        print(f"    → wrote pocket marker at ({cx:.2f}, {cy:.2f}, {cz:.2f})")

    # Run DiffSBDD if raw output doesn't exist yet
    if not os.path.exists(raw_path):
        run_diffsbdd(INPUT_PDB, marker_path, raw_path)

    # Sanitize
    total, kept = sanitize_sdf(raw_path, clean_path)
    pct = 100 * kept / total if total else 0
    print(f"    → sanitized: {kept}/{total} kept ({pct:.0f}%)")

    if kept == 0:
        print(f"    ⚠️  no valid ligands from pocket{pocket_num} — skipping downstream")

    return kept


def main():
    print(f"\n### Step 2: Generate ligands for {TARGET_NAME} ###\n")

    if not os.path.exists(POCKETS_CSV):
        sys.exit(f"ERROR: {POCKETS_CSV} missing. Run 01_prepare.py first.")

    df = pd.read_csv(POCKETS_CSV)
    selected = df[df["selected"] == True].reset_index(drop=True)

    if len(selected) == 0:
        sys.exit("ERROR: no pockets marked 'selected' in summary CSV.")

    print(f"Processing {len(selected)} selected pocket(s):")
    for _, row in selected.iterrows():
        print(f"  pocket{int(row['pocket_num'])}  "
              f"druggability={row['druggability']:.3f}  "
              f"center=({row['center_x']:.2f}, {row['center_y']:.2f}, {row['center_z']:.2f})")
    print()

    results = []
    for _, row in selected.iterrows():
        n = int(row["pocket_num"])
        print(f"--- pocket{n} ---")
        kept = process_pocket(n, row["center_x"], row["center_y"], row["center_z"])
        results.append((n, kept))
        print()

    print("=" * 60)
    print("Generation summary:")
    for n, kept in results:
        print(f"  pocket{n}: {kept} valid ligands")
    print()
    print("Next: python 03_dock_generated.py")


if __name__ == "__main__":
    main()