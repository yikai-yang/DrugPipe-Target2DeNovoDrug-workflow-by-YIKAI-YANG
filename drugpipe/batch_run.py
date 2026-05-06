"""
batch_run.py — Run Steps 3-5 for all targets, robust to per-target failures.

Assumes:
  - config.py reads DRUGPIPE_TARGET and DRUGPIPE_PDB env vars (see edit)
  - Steps 01 and 02 have already been run for every target listed below

Behaviour:
  - Hard 2-hour timeout per step (kills stuck runs, continues to next target)
  - Per-target exception handling — one failure doesn't abort others
  - Full stdout/stderr captured to per-target log files
  - Detects and deletes truncated output CSVs (< 200 bytes) so re-runs are clean
  - Skips targets with missing runs/ or no selected pockets
  - Writes batch_summary.csv with status of every target × step

Run:     python batch_run.py
Monitor: tail -f ../batch_run.log
"""

import os
import subprocess
import sys
import time
from pathlib import Path
import pandas as pd

# ==============================================================
# CONFIGURATION
# ==============================================================
BASE_DIR    = Path(os.environ.get("DRUGPIPE_BASE_DIR", Path(__file__).resolve().parent.parent))
SCRIPTS_DIR = Path(__file__).resolve().parent
RUNS_DIR    = BASE_DIR / "runs"

# Targets to process. IGFBP7 omitted — already completed.
TARGETS = [
    "IGFBP7",                       # resume from where Step 5 was interrupted
    "BGN",
    "BMP1",
    "INHBA_with_propeptide",
    "INHBA_without_propeptide",
    "MMP9",
    "SPARC_monomer",
    "TANGL",
    "TNC",
]

# (script_name, expected per-pocket output suffix for truncation detection)
STEPS = [
    ("03_dock_generated.py", "generated_docking.csv"),
    ("04_morgan_rank.py",    "morgan.csv"),
    ("05_dock_drugbank.py",  "final.csv"),
]

STEP_TIMEOUT_SEC = 7200    # 2 hours per step
MIN_CSV_BYTES    = 200     # smaller than this → treat as corrupt partial

BATCH_LOG     = BASE_DIR / "batch_run.log"
BATCH_SUMMARY = BASE_DIR / "batch_summary.csv"

# ==============================================================


def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(BATCH_LOG, "a") as f:
        f.write(line + "\n")


def get_pdb_path(target):
    """Step 01 copied the original PDB here — use it as canonical."""
    return RUNS_DIR / target / "input" / f"{target}.pdb"


def get_selected_pockets(target):
    """Read pockets_summary.csv from Step 01 output."""
    csv = RUNS_DIR / target / f"{target}_pockets_summary.csv"
    if not csv.exists():
        return None
    df = pd.read_csv(csv)
    sel = df[df["selected"] == True]
    return sel["pocket_num"].astype(int).tolist()


def cleanup_truncated_csvs(target, suffix):
    """Delete per-pocket CSVs for this step that are too small to be real.
    Prevents the idempotent-skip logic in step scripts from trusting a
    corrupt CSV on re-run.
    """
    pockets = get_selected_pockets(target) or []
    deleted = []
    for n in pockets:
        csv = RUNS_DIR / target / f"pocket{n}" / f"{target}_pocket{n}_{suffix}"
        if csv.exists() and csv.stat().st_size < MIN_CSV_BYTES:
            try:
                csv.unlink()
                deleted.append(csv.name)
            except Exception:
                pass
    if deleted:
        log(f"    cleaned truncated CSVs: {deleted}")


def run_step(target, script_name, out_suffix):
    """Invoke one step for one target in a subprocess with env vars set."""
    env = os.environ.copy()
    env["DRUGPIPE_TARGET"] = target
    env["DRUGPIPE_PDB"]    = str(get_pdb_path(target))

    step_log = RUNS_DIR / target / f"batch_{script_name}.log"
    step_log.parent.mkdir(parents=True, exist_ok=True)

    log(f"  → {script_name}")
    t0 = time.time()
    try:
        r = subprocess.run(
            [sys.executable, script_name],
            cwd=SCRIPTS_DIR,
            env=env,
            capture_output=True,
            text=True,
            timeout=STEP_TIMEOUT_SEC,
        )
    except subprocess.TimeoutExpired:
        log(f"    TIMEOUT after {STEP_TIMEOUT_SEC}s")
        with open(step_log, "w") as f:
            f.write(f"TIMEOUT after {STEP_TIMEOUT_SEC}s\n")
        cleanup_truncated_csvs(target, out_suffix)
        return False, "timeout"

    # Persist full log for post-mortem
    with open(step_log, "w") as f:
        f.write(f"=== python {script_name} ===\n")
        f.write(f"=== DRUGPIPE_TARGET={target} ===\n")
        f.write(f"=== duration: {time.time()-t0:.1f}s ===\n\n")
        f.write(f"=== STDOUT ===\n{r.stdout}\n")
        f.write(f"\n=== STDERR ===\n{r.stderr}\n")

    if r.returncode != 0:
        log(f"    FAILED (exit {r.returncode}) — log: runs/{target}/{step_log.name}")
        for line in r.stderr.strip().splitlines()[-3:]:
            log(f"      | {line}")
        cleanup_truncated_csvs(target, out_suffix)
        return False, f"exit_{r.returncode}"

    log(f"    done in {time.time()-t0:.0f}s")
    return True, "ok"


def process_target(target):
    log("")
    log("=" * 60)
    log(f"TARGET: {target}")
    log("=" * 60)

    target_dir = RUNS_DIR / target
    if not target_dir.exists():
        log(f"  SKIP: {target_dir.relative_to(BASE_DIR)} missing (Steps 0-2 not done?)")
        return {"target": target, "status": "skip_missing_dir"}

    pdb = get_pdb_path(target)
    if not pdb.exists():
        log(f"  SKIP: {pdb.relative_to(BASE_DIR)} missing")
        return {"target": target, "status": "skip_missing_pdb"}

    pockets = get_selected_pockets(target)
    if pockets is None:
        log(f"  SKIP: pockets_summary.csv missing")
        return {"target": target, "status": "skip_no_summary"}
    if not pockets:
        log(f"  SKIP: no selected pockets")
        return {"target": target, "status": "skip_no_pockets"}
    log(f"  pockets to process: {pockets}")

    result = {"target": target, "pockets": str(pockets)}
    for script_name, out_suffix in STEPS:
        step_key = script_name.split("_")[0]   # "03", "04", "05"
        ok, msg = run_step(target, script_name, out_suffix)
        result[f"step_{step_key}"] = msg
        if not ok:
            log(f"  abort {target} at {script_name}")
            result["status"] = f"failed_at_{step_key}"
            return result

    result["status"] = "ok"
    log(f"  {target}: COMPLETED all steps")
    return result


def main():
    log("")
    log("#" * 60)
    log(f"BATCH RUN STARTED — {len(TARGETS)} targets")
    log(f"  log file: {BATCH_LOG.relative_to(BASE_DIR)}")
    log("#" * 60)

    batch_t0 = time.time()
    results = []
    for target in TARGETS:
        try:
            r = process_target(target)
        except KeyboardInterrupt:
            log("INTERRUPTED by user — stopping batch")
            break
        except Exception as e:
            log(f"  UNEXPECTED EXCEPTION for {target}: {type(e).__name__}: {e}")
            r = {"target": target, "status": f"unexpected_{type(e).__name__}"}
        results.append(r)

    elapsed_min = (time.time() - batch_t0) / 60
    log("")
    log("#" * 60)
    log(f"BATCH COMPLETE — {elapsed_min:.1f} min elapsed")
    log("#" * 60)

    df = pd.DataFrame(results)
    df.to_csv(BATCH_SUMMARY, index=False)
    log(f"Summary written to: {BATCH_SUMMARY.relative_to(BASE_DIR)}")
    log("")
    log("Final status:")
    for r in results:
        log(f"  {r['target']:<30s} {r['status']}")


if __name__ == "__main__":
    main()