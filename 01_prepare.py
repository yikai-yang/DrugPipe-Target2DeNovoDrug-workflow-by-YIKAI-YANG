"""
01_prepare.py — Prepare protein receptor and find druggable pockets.

Reads from config.py:
    - PDB_FILE, TARGET_NAME, N_POCKETS
    - INPUT_PDB, RECEPTOR_PDBQT, FPOCKET_DIR, POCKETS_CSV

Outputs:
    - runs/<TARGET>/input/<TARGET>.pdb                  (copy of input)
    - runs/<TARGET>/receptor/<TARGET>.pdbqt             (prepared receptor)
    - runs/<TARGET>/fpocket/<TARGET>_out/               (fpocket output)
    - runs/<TARGET>/<TARGET>_pockets_summary.csv        (ranked pockets)

Run once per target. Safe to re-run: skips work that's already done.
"""

import os
import re
import shutil
import subprocess
import sys
import pandas as pd

from config import (
    TARGET_NAME, PDB_FILE, N_POCKETS, MIN_DRUGGABILITY,
    INPUT_PDB, RECEPTOR_DIR, RECEPTOR_PDBQT,
    TARGET_DIR, FPOCKET_DIR, POCKETS_CSV,
)


def copy_input_pdb():
    """Put a copy of the input PDB inside the target folder."""
    if os.path.exists(INPUT_PDB):
        print(f"[skip]   Input PDB already copied: {INPUT_PDB}")
        return
    if not os.path.exists(PDB_FILE):
        sys.exit(f"ERROR: input PDB not found: {PDB_FILE}")
    shutil.copy(PDB_FILE, INPUT_PDB)
    print(f"[done]   Copied {PDB_FILE} → {INPUT_PDB}")


def strip_waters_and_hetatm(pdb_in, pdb_out):
    """Write a cleaned PDB with HOH and non-standard HETATMs removed.
    Keeps standard residues and cofactors like metals if present.
    """
    keep_hetatm = {"MG", "ZN", "CA", "FE", "MN", "CU", "NA", "K"}  # common metal ions
    with open(pdb_in) as f_in, open(pdb_out, "w") as f_out:
        for line in f_in:
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname == "HOH":
                    continue
                if resname not in keep_hetatm:
                    continue
            f_out.write(line)


def prepare_receptor():
    """Convert PDB to PDBQT using Open Babel."""
    if os.path.exists(RECEPTOR_PDBQT):
        print(f"[skip]   Receptor PDBQT already exists: {RECEPTOR_PDBQT}")
        return

    os.makedirs(RECEPTOR_DIR, exist_ok=True)
    cleaned = os.path.join(RECEPTOR_DIR, f"{TARGET_NAME}_cleaned.pdb")
    strip_waters_and_hetatm(INPUT_PDB, cleaned)

    # -xr = rigid receptor; -p 7.4 = protonation at physiological pH;
    # --partialcharge gasteiger = add partial charges
    cmd = [
        "obabel", cleaned,
        "-O", RECEPTOR_PDBQT,
        "-xr",
        "-p", "7.4",
        "--partialcharge", "gasteiger",
    ]
    print(f"[run]    obabel receptor prep...")
    try:
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        sys.exit("ERROR: obabel failed. Check that Open Babel is installed: `which obabel`")

    if not os.path.exists(RECEPTOR_PDBQT) or os.path.getsize(RECEPTOR_PDBQT) < 1000:
        sys.exit(f"ERROR: receptor PDBQT is empty or tiny: {RECEPTOR_PDBQT}")

    os.remove(cleaned)
    print(f"[done]   Receptor PDBQT: {RECEPTOR_PDBQT}")


def run_fpocket():
    """Run fpocket on the prepared PDB.
    Fpocket creates a sibling folder named <pdb_basename>_out next to the input.
    We move/symlink it into our structured runs/ layout.
    """
    # fpocket will write into the same directory as the input pdb.
    # We run it inside INPUT_DIR so the output lands next to the copy.
    input_dir = os.path.dirname(INPUT_PDB)
    fpocket_raw_out = os.path.join(input_dir, f"{TARGET_NAME}_out")

    target_fpocket_dir = os.path.join(FPOCKET_DIR, f"{TARGET_NAME}_out")

    if os.path.exists(target_fpocket_dir):
        print(f"[skip]   Fpocket output already exists: {target_fpocket_dir}")
        return target_fpocket_dir

    os.makedirs(FPOCKET_DIR, exist_ok=True)

    print(f"[run]    fpocket -f {INPUT_PDB}")
    try:
        subprocess.check_call(
            ["fpocket", "-f", INPUT_PDB],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError:
        sys.exit("ERROR: fpocket failed. Check installation: `which fpocket`")

    if not os.path.exists(fpocket_raw_out):
        sys.exit(f"ERROR: fpocket didn't produce expected output at {fpocket_raw_out}")

    # Move fpocket output into our structured location
    shutil.move(fpocket_raw_out, target_fpocket_dir)
    print(f"[done]   Fpocket output: {target_fpocket_dir}")
    return target_fpocket_dir


def parse_pocket_info(info_file):
    """Parse fpocket's <target>_info.txt file into a list of dicts.
    One dict per pocket with keys: pocket_num, score, druggability, volume,
    hydrophobicity, polarity_score, num_alpha_spheres.
    """
    if not os.path.exists(info_file):
        sys.exit(f"ERROR: fpocket info file missing: {info_file}")

    with open(info_file) as f:
        content = f.read()

    # Split by "Pocket N :" headers
    pocket_blocks = re.split(r"Pocket (\d+) :\s*\n", content)[1:]
    # Format: [num1, block1, num2, block2, ...]

    pockets = []
    for i in range(0, len(pocket_blocks), 2):
        num = int(pocket_blocks[i])
        block = pocket_blocks[i + 1]

        def grab(key):
            m = re.search(rf"{re.escape(key)}\s*:\s*(-?\d+\.?\d*)", block)
            return float(m.group(1)) if m else None

        pockets.append({
            "pocket_num":         num,
            "score":              grab("Score"),
            "druggability":       grab("Druggability Score"),
            "num_alpha_spheres":  grab("Number of Alpha Spheres"),
            "total_sasa":         grab("Total SASA"),
            "volume":              grab("Volume"),
            "hydrophobicity":     grab("Mean local hydrophobic density"),
            "polarity":           grab("Polarity score"),
        })
    return pockets


def compute_pocket_center(pocket_atm_pdb):
    """Compute centroid of all atoms in a pocket atm.pdb file."""
    xs, ys, zs = [], [], []
    with open(pocket_atm_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
    if not xs:
        return None, None, None
    return sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)


def build_pockets_summary(fpocket_out_dir):
    """Parse pocket info + compute centers, write ranked CSV."""
    info_file = os.path.join(fpocket_out_dir, f"{TARGET_NAME}_info.txt")
    pockets = parse_pocket_info(info_file)

    # Add center coordinates for each pocket
    pockets_dir = os.path.join(fpocket_out_dir, "pockets")
    for p in pockets:
        n = p["pocket_num"]
        atm = os.path.join(pockets_dir, f"pocket{n}_atm.pdb")
        cx, cy, cz = compute_pocket_center(atm)
        p["center_x"] = round(cx, 4) if cx is not None else None
        p["center_y"] = round(cy, 4) if cy is not None else None
        p["center_z"] = round(cz, 4) if cz is not None else None

    df = pd.DataFrame(pockets)

    # Rank by druggability (primary) then score (tiebreaker)
    df = df.sort_values(
        by=["druggability", "score"],
        ascending=[False, False],
        na_position="last"
    ).reset_index(drop=True)
    df["druggability_rank"] = range(1, len(df) + 1)

    # Flag which pockets we'll actually use downstream
    df["selected"] = False
    top = df[df["druggability"] >= MIN_DRUGGABILITY].head(N_POCKETS)
    df.loc[top.index, "selected"] = True

    # Write summary
    cols = [
        "druggability_rank", "pocket_num", "selected",
        "druggability", "score", "volume",
        "num_alpha_spheres", "hydrophobicity", "polarity",
        "center_x", "center_y", "center_z",
    ]
    df[cols].to_csv(POCKETS_CSV, index=False)
    return df


def print_summary(df):
    """Human-readable printout of the ranking."""
    print()
    print("=" * 80)
    print(f"Pocket ranking for {TARGET_NAME} (top {N_POCKETS} selected)")
    print("=" * 80)
    shown_cols = ["druggability_rank", "pocket_num", "selected",
                  "druggability", "volume",
                  "center_x", "center_y", "center_z"]
    print(df[shown_cols].head(10).to_string(index=False))
    print()
    selected_count = df["selected"].sum()
    print(f"Selected {selected_count} pocket(s) for downstream DiffSBDD + docking:")
    for _, row in df[df["selected"]].iterrows():
        print(f"  pocket{int(row['pocket_num']):<3d}  "
              f"druggability={row['druggability']:.3f}  "
              f"center=({row['center_x']:.2f}, {row['center_y']:.2f}, {row['center_z']:.2f})")

    # Warn if druggability is low across the board
    best = df["druggability"].max()
    if best < 0.3:
        print()
        print(f"⚠️  WARNING: best druggability is only {best:.3f}.")
        print(f"    This protein may have no highly druggable sites.")
        print(f"    Results will likely be modest. Consider:")
        print(f"      - Checking the structure is correct / properly protonated")
        print(f"      - Using MIN_DRUGGABILITY=0 (current: {MIN_DRUGGABILITY}) to process anyway")


def main():
    print(f"\n### Step 1: Prepare {TARGET_NAME} ###\n")

    copy_input_pdb()
    prepare_receptor()
    fpocket_out_dir = run_fpocket()
    df = build_pockets_summary(fpocket_out_dir)
    print_summary(df)

    print()
    print(f"Pocket summary written to: {POCKETS_CSV}")
    print(f"Review the ranking, then run:  python 02_generate.py")


if __name__ == "__main__":
    main()