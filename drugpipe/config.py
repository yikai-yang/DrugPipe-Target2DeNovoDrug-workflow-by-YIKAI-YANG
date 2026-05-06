"""
config.py — DrugPipe workflow configuration

This file is imported by every step (01-06). Edit the TARGET section for
each new protein you run; leave everything else alone unless you're
changing machines.

Usage from other scripts:
    from config import *
"""

import os
from pathlib import Path

# ==============================================================
# TARGET — EDIT THIS SECTION FOR EACH NEW PROTEIN
# ==============================================================
BASE_DIR       = os.environ.get(
    "DRUGPIPE_BASE_DIR",
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)

TARGET_NAME = os.environ.get("DRUGPIPE_TARGET", "SPARC_monomer")  # e.g. "IGFBP7", "BGN", etc.
PDB_FILE    = os.environ.get("DRUGPIPE_PDB", "")  # set DRUGPIPE_PDB or edit _sanity_check will warn
# ==============================================================
# PIPELINE PARAMETERS — adjust based on time budget
# ==============================================================
# Pocket selection
N_POCKETS            = 2        # top-N druggable pockets to process (3-10)
MIN_DRUGGABILITY     = 0.0      # set >0 to skip low-quality pockets; 0 = take top-N regardless

# DiffSBDD generation
N_LIGANDS            = 30      # samples per pocket (expect 30-70% valid after sanitization)
DIFFSBDD_TIMESTEPS   = 500      # more = slower but better quality; 500 is DiffSBDD default
DIFFSBDD_BATCH_SIZE  = 10       # increase if you have GPU memory

# Docking
BOX_SIZE             = 22.0     # Å, cubic box side (22 = good default)
EXHAUSTIVENESS       = 8        # Vina default; lower (4) = faster, (16) = more thorough
N_POSES              = 5        # poses Vina reports per ligand
DOCKING_TIMEOUT      = 600      # seconds per dock; kill runaway jobs

# DrugBank screening
TOP_N_DOCK           = 100      # top Morgan FP hits per pocket to dock
MORGAN_RADIUS        = 2        # ECFP4 (radius 2) is standard for drug similarity
MORGAN_NBITS         = 2048

# Parallelization
N_WORKERS            = 8        # parallel docking jobs (your MacBook has 10 cores)

# ==============================================================
# SHARED ASSETS — edit once per machine, not per target
# ==============================================================
SHARED_DIR     = f"{BASE_DIR}/shared"

DRUGBANK_CSV   = f"{SHARED_DIR}/drugbank.csv"
DIFFSBDD_DIR   = f"{SHARED_DIR}/DiffSBDD"
DIFFSBDD_CKPT  = "checkpoints/crossdocked_fullatom_cond.ckpt"  # relative to DIFFSBDD_DIR

# ==============================================================
# DERIVED PATHS — auto-generated, don't edit
# ==============================================================
RUNS_DIR = os.environ.get("DRUGPIPE_RUNS_DIR", f"{BASE_DIR}/runs")
TARGET_DIR      = f"{RUNS_DIR}/{TARGET_NAME}"

INPUT_DIR       = f"{TARGET_DIR}/input"
INPUT_PDB       = f"{INPUT_DIR}/{TARGET_NAME}.pdb"

RECEPTOR_DIR    = f"{TARGET_DIR}/receptor"
RECEPTOR_PDBQT  = f"{RECEPTOR_DIR}/{TARGET_NAME}.pdbqt"

FPOCKET_DIR     = f"{TARGET_DIR}/fpocket"
POCKETS_CSV     = f"{TARGET_DIR}/{TARGET_NAME}_pockets_summary.csv"

SUMMARY_DIR     = f"{TARGET_DIR}/summary"
FINAL_RANKING   = f"{SUMMARY_DIR}/{TARGET_NAME}_final_ranking.csv"

def pocket_dir(n):
    """Return path to pocket N's folder: runs/IGFBP7/pocket3/ etc."""
    return f"{TARGET_DIR}/pocket{n}"

def pocket_file(n, suffix):
    """Return a target-and-pocket-tagged file path.
    Example: pocket_file(3, 'ligands.sdf') -> runs/IGFBP7/pocket3/IGFBP7_pocket3_ligands.sdf
    """
    return f"{pocket_dir(n)}/{TARGET_NAME}_pocket{n}_{suffix}"

# ==============================================================
# INITIALIZATION — runs on every import
# ==============================================================
def _ensure_dirs():
    """Create target folder structure if missing (safe to re-run)."""
    for d in [TARGET_DIR, INPUT_DIR, RECEPTOR_DIR, SUMMARY_DIR]:
        os.makedirs(d, exist_ok=True)

def _sanity_check():
    """Verify shared assets exist; warn if target inputs are missing."""
    missing = []
    if not os.path.exists(DRUGBANK_CSV):
        missing.append(f"DrugBank CSV: {DRUGBANK_CSV}")
    if not os.path.exists(DIFFSBDD_DIR):
        missing.append(f"DiffSBDD repo: {DIFFSBDD_DIR}")
    ckpt_full = os.path.join(DIFFSBDD_DIR, DIFFSBDD_CKPT)
    if not os.path.exists(ckpt_full):
        missing.append(f"DiffSBDD checkpoint: {ckpt_full}")
    if not os.path.exists(PDB_FILE):
        missing.append(f"Target PDB: {PDB_FILE}")
    if missing:
        print("⚠️  WARNING: missing files/folders:")
        for m in missing:
            print(f"    - {m}")

_ensure_dirs()

if __name__ == "__main__":
    # Run `python config.py` to print current settings and sanity-check.
    print(f"\n{'='*60}")
    print(f"DrugPipe configuration — target: {TARGET_NAME}")
    print(f"{'='*60}")
    print(f"Input PDB:       {PDB_FILE}")
    print(f"Output folder:   {TARGET_DIR}")
    print(f"")
    print(f"Pockets:         top {N_POCKETS} by druggability")
    print(f"Ligands/pocket:  {N_LIGANDS}")
    print(f"Docking box:     {BOX_SIZE}³ Å, exhaustiveness {EXHAUSTIVENESS}")
    print(f"DrugBank top-N:  {TOP_N_DOCK} per pocket")
    print(f"Parallel workers:{N_WORKERS}")
    print(f"")
    print(f"Shared assets:")
    print(f"  DrugBank CSV:  {DRUGBANK_CSV}")
    print(f"  DiffSBDD:      {DIFFSBDD_DIR}")
    print(f"  Checkpoint:    {DIFFSBDD_CKPT}")
    print()
    _sanity_check()