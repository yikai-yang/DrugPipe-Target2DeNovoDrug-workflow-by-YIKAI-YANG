"""
rerun_SPARC.py — Isolated rerun of all pipeline stages for SPARC_monomer.

Purpose: diagnose the 0.0 kcal/mol docking anomaly seen in the original run.

Identical parameters to the first run, EXCEPT:
  - Output tree is in runs_SPARC_rerun/ (sibling to runs/), not touching originals
  - AutoDock Vina is invoked with --seed 42 for deterministic diagnosis
  - Extra validation/diagnostics at each stage, written to DIAGNOSTICS.md

Run:
  conda activate sbdd-env
  cd <repo_root>/drugpipe
  caffeinate -i python rerun_SPARC.py
"""

import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
import pandas as pd


# ==============================================================
# ISOLATED PATHS — everything is under RERUN_ROOT, nothing touches original
# ==============================================================
BASE_DIR       = Path(os.environ.get("DRUGPIPE_BASE_DIR", Path(__file__).resolve().parent.parent))
SCRIPTS_DIR    = Path(__file__).resolve().parent
SHARED_DIR     = BASE_DIR / "shared"

# Original (untouched — used only as the PDB source)
ORIGINAL_PDB   = BASE_DIR / "PDB" / "SPARC_1BMO_monomer_energy_minimized copy.pdb"

# Rerun tree
RERUN_ROOT     = BASE_DIR / "runs_SPARC_rerun"
TARGET_NAME    = "SPARC_monomer"
TARGET_DIR     = RERUN_ROOT / TARGET_NAME
INPUT_DIR      = TARGET_DIR / "input"
INPUT_PDB      = INPUT_DIR / f"{TARGET_NAME}.pdb"
RECEPTOR_DIR   = TARGET_DIR / "receptor"
RECEPTOR_PDBQT = RECEPTOR_DIR / f"{TARGET_NAME}.pdbqt"
FPOCKET_DIR    = TARGET_DIR / "fpocket"
POCKETS_CSV    = TARGET_DIR / f"{TARGET_NAME}_pockets_summary.csv"
SUMMARY_DIR    = TARGET_DIR / "summary"
DIAG_MD        = TARGET_DIR / "DIAGNOSTICS.md"
RERUN_LOG      = TARGET_DIR / "rerun.log"

# ==============================================================
# PARAMETERS — IDENTICAL to original run
# ==============================================================
N_POCKETS         = 2
MIN_DRUGGABILITY  = 0.0
N_LIGANDS         = 30
DIFFSBDD_TIMESTEPS = 500
BOX_SIZE          = 22.0
EXHAUSTIVENESS    = 8
N_POSES           = 5
DOCKING_TIMEOUT   = 600
TOP_N_DOCK        = 100
MORGAN_RADIUS     = 2
MORGAN_NBITS      = 2048
N_WORKERS         = 8

# Diagnostic-only change
VINA_SEED = 42

# Shared assets (read-only)
DRUGBANK_CSV   = SHARED_DIR / "drugbank.csv"
DIFFSBDD_DIR   = SHARED_DIR / "DiffSBDD"
DIFFSBDD_CKPT  = "checkpoints/crossdocked_fullatom_cond.ckpt"


# ==============================================================
# LOGGING
# ==============================================================
DIAG_LINES = []


def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    RERUN_LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(RERUN_LOG, "a") as f:
        f.write(line + "\n")


def diag(section, msg):
    """Append to the diagnostics markdown being built up."""
    DIAG_LINES.append(f"### {section}\n\n{msg}\n")
    log(f"[diag] {section}: {msg.splitlines()[0][:100]}")


def write_diagnostics():
    DIAG_MD.parent.mkdir(parents=True, exist_ok=True)
    with open(DIAG_MD, "w") as f:
        f.write(f"# SPARC_monomer Rerun Diagnostics\n\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Original PDB: `{ORIGINAL_PDB}`\n\n")
        f.write(f"Rerun tree: `{RERUN_ROOT}`\n\n")
        f.write("---\n\n")
        for entry in DIAG_LINES:
            f.write(entry + "\n")
    log(f"Diagnostics written: {DIAG_MD}")


# ==============================================================
# HELPERS
# ==============================================================
def run_cmd(cmd, **kwargs):
    """Run a subprocess, raising on failure, capturing output."""
    r = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    return r


def receptor_atom_stats(pdbqt_path):
    """Return (n_atoms, min_xyz, max_xyz, centroid)."""
    xs, ys, zs = [], [], []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except ValueError:
                    continue
    if not xs:
        return 0, None, None, None
    n = len(xs)
    mn = (min(xs), min(ys), min(zs))
    mx = (max(xs), max(ys), max(zs))
    ctr = (sum(xs)/n, sum(ys)/n, sum(zs)/n)
    return n, mn, mx, ctr


# ==============================================================
# STAGE 1 — Protein preparation + pocket detection
# ==============================================================
def stage1_prepare():
    log("")
    log("=" * 70)
    log("STAGE 1 — Protein prep + fpocket")
    log("=" * 70)

    # Create isolated tree
    for d in [TARGET_DIR, INPUT_DIR, RECEPTOR_DIR, SUMMARY_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    # Copy PDB
    if not ORIGINAL_PDB.exists():
        sys.exit(f"ERROR: source PDB not found: {ORIGINAL_PDB}")
    shutil.copy(ORIGINAL_PDB, INPUT_PDB)
    log(f"  copied PDB: {ORIGINAL_PDB.name} -> {INPUT_PDB}")

    # Clean receptor for docking:
    # Keep protein ATOM records only. Remove glycans/ligands/metals as HETATM,
    # plus LINK/CONECT/SSBOND/MODRES/etc that can confuse Open Babel.
    cleaned = RECEPTOR_DIR / f"{TARGET_NAME}_cleaned.pdb"
    n_in, n_out, n_het_removed = 0, 0, 0

    with open(INPUT_PDB) as fi, open(cleaned, "w") as fo:
        for line in fi:
            n_in += 1

            if line.startswith("ATOM  "):
                fo.write(line)
                n_out += 1

            elif line.startswith("TER"):
                fo.write(line)
                n_out += 1

            elif line.startswith("HETATM"):
                n_het_removed += 1

        fo.write("END\n")
        n_out += 1

    log(f"  cleaned PDB: kept {n_out}/{n_in} lines "
        f"(removed {n_het_removed} HETATM; removed LINK/CONECT/SSBOND/etc)")

    # Convert to PDBQT
    cmd = [
        "obabel", "-ipdb", str(cleaned),
        "-opdbqt", "-O", str(RECEPTOR_PDBQT),
        "-xr", "-p", "7.4", "--partialcharge", "gasteiger",
    ]
    log(f"  running obabel...")
    r = run_cmd(cmd)
    if r.returncode != 0 or not RECEPTOR_PDBQT.exists():
        sys.exit(f"ERROR: obabel failed.\n  stderr: {r.stderr[:500]}")
    # cleaned.unlink()  # keep for debugging

    # DIAGNOSTIC — receptor validity
    n_atoms, mn, mx, ctr = receptor_atom_stats(RECEPTOR_PDBQT)
    size = RECEPTOR_PDBQT.stat().st_size
    if n_atoms == 0:
        diag("receptor_pdbqt", f"**FAIL** — PDBQT has 0 atoms (size={size} bytes)")
        sys.exit("ERROR: receptor PDBQT has no atoms — obabel produced an empty output")
    extents = (mx[0]-mn[0], mx[1]-mn[1], mx[2]-mn[2])
    diag("receptor_pdbqt",
         f"- atoms: **{n_atoms}**\n"
         f"- file size: {size:,} bytes\n"
         f"- bounding box: x∈[{mn[0]:.1f}, {mx[0]:.1f}], "
         f"y∈[{mn[1]:.1f}, {mx[1]:.1f}], z∈[{mn[2]:.1f}, {mx[2]:.1f}]\n"
         f"- extents: ({extents[0]:.1f}, {extents[1]:.1f}, {extents[2]:.1f}) Å\n"
         f"- centroid: ({ctr[0]:.2f}, {ctr[1]:.2f}, {ctr[2]:.2f})")

    # Run fpocket
    log(f"  running fpocket...")
    r = run_cmd(["fpocket", "-f", str(INPUT_PDB)])
    if r.returncode != 0:
        sys.exit(f"ERROR: fpocket failed.\n  stderr: {r.stderr[:500]}")

    fpocket_raw = INPUT_DIR / f"{TARGET_NAME}_out"
    if not fpocket_raw.exists():
        sys.exit(f"ERROR: fpocket output missing: {fpocket_raw}")

    FPOCKET_DIR.mkdir(parents=True, exist_ok=True)
    dest = FPOCKET_DIR / f"{TARGET_NAME}_out"
    if dest.exists():
        shutil.rmtree(dest)
    shutil.move(str(fpocket_raw), str(dest))
    log(f"  fpocket output moved to: {dest}")

    # Parse pocket info
    info_file = dest / f"{TARGET_NAME}_info.txt"
    if not info_file.exists():
        sys.exit(f"ERROR: fpocket info file missing: {info_file}")
    with open(info_file) as f:
        content = f.read()

    pocket_blocks = re.split(r"Pocket (\d+) :\s*\n", content)[1:]
    pockets = []
    for i in range(0, len(pocket_blocks), 2):
        num = int(pocket_blocks[i])
        block = pocket_blocks[i + 1]
        def grab(key):
            m = re.search(rf"{re.escape(key)}\s*:\s*(-?\d+\.?\d*)", block)
            return float(m.group(1)) if m else None
        pockets.append({
            "pocket_num":        num,
            "score":             grab("Score"),
            "druggability":      grab("Druggability Score"),
            "num_alpha_spheres": grab("Number of Alpha Spheres"),
            "total_sasa":        grab("Total SASA"),
            "volume":            grab("Volume"),
            "hydrophobicity":    grab("Mean local hydrophobic density"),
            "polarity":          grab("Polarity score"),
        })

    # Compute centroids from atm.pdb files
    pockets_dir = dest / "pockets"
    for p in pockets:
        atm = pockets_dir / f"pocket{p['pocket_num']}_atm.pdb"
        xs, ys, zs = [], [], []
        if atm.exists():
            with open(atm) as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        try:
                            xs.append(float(line[30:38]))
                            ys.append(float(line[38:46]))
                            zs.append(float(line[46:54]))
                        except ValueError:
                            continue
        if xs:
            p["center_x"] = round(sum(xs)/len(xs), 4)
            p["center_y"] = round(sum(ys)/len(ys), 4)
            p["center_z"] = round(sum(zs)/len(zs), 4)
        else:
            p["center_x"] = p["center_y"] = p["center_z"] = None

    df = pd.DataFrame(pockets).sort_values(
        by=["druggability", "score"], ascending=[False, False], na_position="last"
    ).reset_index(drop=True)
    df["druggability_rank"] = range(1, len(df) + 1)
    df["selected"] = False
    top = df[df["druggability"] >= MIN_DRUGGABILITY].head(N_POCKETS)
    df.loc[top.index, "selected"] = True

    cols = ["druggability_rank", "pocket_num", "selected",
            "druggability", "score", "volume", "num_alpha_spheres",
            "hydrophobicity", "polarity", "center_x", "center_y", "center_z"]
    df[cols].to_csv(POCKETS_CSV, index=False)

    # DIAGNOSTIC — pocket centroid vs receptor bounding box
    sel = df[df["selected"]].reset_index(drop=True)
    lines = [f"Selected {len(sel)} pockets (top {N_POCKETS} by druggability):\n"]
    anomalies = []
    for _, row in sel.iterrows():
        n = int(row['pocket_num'])
        cx, cy, cz = row['center_x'], row['center_y'], row['center_z']
        drg = row['druggability']
        half = BOX_SIZE / 2
        # Check: does the docking box (22 Å cube around centroid) overlap the receptor?
        box_min = (cx - half, cy - half, cz - half)
        box_max = (cx + half, cy + half, cz + half)
        overlaps = all(box_min[i] <= mx[i] and box_max[i] >= mn[i] for i in range(3))
        # Check: is centroid inside receptor bounding box?
        inside = all(mn[i] <= (cx, cy, cz)[i] <= mx[i] for i in range(3))
        status = "✓" if overlaps else "✗ BOX OUTSIDE RECEPTOR"
        lines.append(f"- pocket{n}: druggability={drg:.3f}, "
                     f"centroid=({cx:.2f}, {cy:.2f}, {cz:.2f}), "
                     f"inside_receptor={inside}, box_overlaps={overlaps} {status}")
        if not overlaps:
            anomalies.append(f"pocket{n} docking box does not overlap receptor")
        elif not inside:
            anomalies.append(f"pocket{n} centroid outside receptor bounding box (but box overlaps)")

    diag("pocket_selection", "\n".join(lines))
    if anomalies:
        diag("POCKET_ANOMALIES_FOUND",
             "⚠ The following issues may explain the 0.0 Vina scores:\n- " +
             "\n- ".join(anomalies))

    log(f"  wrote pocket summary: {POCKETS_CSV}")
    return df


# ==============================================================
# STAGE 2–5 — Reuse existing scripts via subprocess with env overrides
# ==============================================================
def run_pipeline_step(script_name, step_label):
    """Run one of the numbered step scripts against our isolated config.
    
    We dynamically rewrite the paths by setting DRUGPIPE_TARGET, DRUGPIPE_PDB,
    AND patching the run directory via environment.
    """
    log("")
    log("=" * 70)
    log(f"{step_label} — {script_name}")
    log("=" * 70)

    env = os.environ.copy()
    env["DRUGPIPE_TARGET"] = TARGET_NAME
    env["DRUGPIPE_PDB"]    = str(INPUT_PDB)
    # Critical: redirect the runs/ base so step scripts write into our isolated tree
    env["DRUGPIPE_RUNS_DIR"] = str(RERUN_ROOT)
    # Tell step 3 and 5 to seed Vina
    env["DRUGPIPE_VINA_SEED"] = str(VINA_SEED)

    t0 = time.time()
    r = subprocess.run(
        [sys.executable, script_name],
        cwd=SCRIPTS_DIR, env=env,
        capture_output=True, text=True,
        timeout=7200,
    )
    dur = time.time() - t0

    # Dump output
    step_log = TARGET_DIR / f"{step_label}_{script_name}.log"
    with open(step_log, "w") as f:
        f.write(f"=== {script_name} ===\n")
        f.write(f"=== duration: {dur:.1f}s ===\n\n")
        f.write(f"=== STDOUT ===\n{r.stdout}\n\n=== STDERR ===\n{r.stderr}\n")

    if r.returncode != 0:
        log(f"  FAILED (exit {r.returncode}). See {step_log.name}")
        for line in r.stderr.strip().splitlines()[-5:]:
            log(f"    | {line}")
        return False

    log(f"  done in {dur:.0f}s")
    return True


# ==============================================================
# STAGE 3 POST-VALIDATION — Catch the 0.0 anomaly immediately
# ==============================================================
def validate_stage3():
    """After step 3 runs, verify that Vina produced real scores, not 0.0."""
    log("")
    log("Validating Stage 3 outputs...")

    df_pockets = pd.read_csv(POCKETS_CSV)
    sel = df_pockets[df_pockets["selected"] == True]

    lines = []
    anomaly = False
    for _, row in sel.iterrows():
        n = int(row["pocket_num"])
        csv = TARGET_DIR / f"pocket{n}" / f"{TARGET_NAME}_pocket{n}_generated_docking.csv"
        if not csv.exists():
            lines.append(f"- pocket{n}: **NO OUTPUT CSV** ({csv})")
            anomaly = True
            continue
        df = pd.read_csv(csv)
        if len(df) == 0:
            lines.append(f"- pocket{n}: CSV exists but empty")
            anomaly = True
            continue
        valid = df["ba"].dropna()
        if len(valid) == 0:
            lines.append(f"- pocket{n}: 0 valid scores, all NaN")
            anomaly = True
            continue
        n_zero = (valid == 0.0).sum()
        n_neg  = (valid < 0).sum()
        lines.append(f"- pocket{n}: {len(valid)} scored, best={valid.min():.3f}, "
                     f"median={valid.median():.3f}, n_zero={n_zero}, n_negative={n_neg}")
        if n_neg == 0:
            anomaly = True
            lines.append(f"  ⚠ **ANOMALY**: pocket{n} has no negative Vina scores — "
                         f"same bug as original run")

    diag("stage3_validation", "\n".join(lines))
    if anomaly:
        diag("STAGE3_ANOMALY",
             "The 0.0 anomaly is reproducing. This points to either:\n"
             "1. Receptor PDBQT geometry problem (check receptor_pdbqt diagnostic above)\n"
             "2. Pocket centroid off-structure (check pocket_selection diagnostic)\n"
             "3. A deeper issue with how Vina interprets the SPARC receptor\n\n"
             "Inspect the per-pocket Vina log files in pocket*/generated_poses/*.log for clues.")
        return False
    return True


# ==============================================================
# STAGE 7 — Run analysis on rerun only
# ==============================================================
def run_analysis():
    log("")
    log("=" * 70)
    log("STAGE 7 — Analysis")
    log("=" * 70)

    env = os.environ.copy()
    env["DRUGPIPE_RUNS_DIR"] = str(RERUN_ROOT)

    r = subprocess.run(
        [sys.executable, "07_analyze.py", TARGET_NAME],
        cwd=SCRIPTS_DIR, env=env,
        capture_output=True, text=True,
    )
    log_file = TARGET_DIR / "stage7_07_analyze.log"
    with open(log_file, "w") as f:
        f.write(f"=== STDOUT ===\n{r.stdout}\n\n=== STDERR ===\n{r.stderr}\n")

    if r.returncode != 0:
        log(f"  analysis failed: see {log_file.name}")
        return False

    log(f"  analysis done")
    return True


# ==============================================================
# MAIN
# ==============================================================
def main():
    if RERUN_LOG.exists():
        RERUN_LOG.unlink()

    log("#" * 70)
    log(f"SPARC_monomer ISOLATED RERUN — {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"  source PDB: {ORIGINAL_PDB}")
    log(f"  output tree: {RERUN_ROOT} (separate from runs/)")
    log(f"  Vina seed: {VINA_SEED}")
    log(f"  all other parameters identical to original run")
    log("#" * 70)

    t0 = time.time()

    # Stage 1 — prepare and validate receptor + pockets
    try:
        stage1_prepare()
    except SystemExit as e:
        log(f"Stage 1 ABORTED: {e}")
        write_diagnostics()
        return

    # Stage 2 — generate ligands
    if not run_pipeline_step("02_generate.py", "stage2"):
        write_diagnostics()
        return

    # Stage 3 — dock generated ligands
    if not run_pipeline_step("03_dock_generated.py", "stage3"):
        write_diagnostics()
        return

    # Stage 3 post-validation — catch the 0.0 anomaly early
    stage3_ok = validate_stage3()
    if not stage3_ok:
        log("")
        log("⚠ Stage 3 validation found the 0.0 anomaly reproducing.")
        log("  Continuing with Stage 4-5 anyway so diagnostics are complete.")

    # Stage 4 — Morgan rank
    if not run_pipeline_step("04_morgan_rank.py", "stage4"):
        write_diagnostics()
        return

    # Stage 5 — dock DrugBank
    if not run_pipeline_step("05_dock_drugbank.py", "stage5"):
        write_diagnostics()
        return

    # Stage 6 — aggregate
    if (SCRIPTS_DIR / "06_aggregate.py").exists():
        run_pipeline_step("06_aggregate.py", "stage6")

    # Stage 7 — analyze
    run_analysis()

    elapsed_min = (time.time() - t0) / 60
    log("")
    log("#" * 70)
    log(f"RERUN COMPLETE — {elapsed_min:.1f} min elapsed")
    log(f"  diagnostics: {DIAG_MD}")
    log(f"  outputs in: {TARGET_DIR}")
    log("#" * 70)

    write_diagnostics()


if __name__ == "__main__":
    main()