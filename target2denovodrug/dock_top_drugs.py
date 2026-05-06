"""
dock_top_drugs.py — Blind-dock the top-10 CIDs per Degree column from
BMP1_drug_screening_results.xlsx and MMP9_drug_screening_results.xlsx
against their respective proteins.

Inputs (expected in INPUT_DIR):
  - BMP1_drug_screening_results.xlsx        (cols: Degree 1, Degree 2, Degree 3 -> CIDs)
  - BMP1 copy.pdb
  - MMP9_drug_screening_results.xlsx
  - MMP9_energy_minimized copy.pdb

Outputs (under OUTPUT_DIR):
  runs/
    BMP1/
      receptor/BMP1.pdbqt
      ligands/
        CID<XXXX>.sdf
        CID<XXXX>.pdbqt
      poses/
        CID<XXXX>_out.pdbqt
        CID<XXXX>_vina.log
      BMP1_docking_results.xlsx
    MMP9/
      (same structure)
  all_targets_docking_summary.xlsx   <- combined workbook (one sheet per target)

Usage:
  conda activate sbdd-env
  python dock_top_drugs.py
"""

import os
import re
import sys
import time
import shutil
import subprocess
import urllib.parse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from rdkit import Chem
from rdkit.Chem import AllChem


# =========================================================================
# CONFIG — edit INPUT_DIR / OUTPUT_DIR to match your layout
# =========================================================================
INPUT_DIR  = Path(os.environ.get("T2D_INPUT_DIR",  "./drug_screening_results"))
OUTPUT_DIR = Path(os.environ.get("T2D_OUTPUT_DIR", "./drug_screening_results/docking_out"))

TARGETS = {
    "BMP1": {
        "pdb":   INPUT_DIR / "BMP1 copy.pdb",
        "xlsx":  INPUT_DIR / "BMP1_drug_screening_results.xlsx",
    },
    "MMP9": {
        "pdb":   INPUT_DIR / "MMP9_energy_minimized copy.pdb",
        "xlsx":  INPUT_DIR / "MMP9_drug_screening_results.xlsx",
    },
}

TOP_N_PER_DEGREE = 10
DEGREE_COLUMNS   = ["Degree 1", "Degree 2", "Degree 3"]

# Blind-dock box — sized to cover full protein, see compute_blind_box()
BOX_PADDING     = 5.0        # Å padding around receptor bounding box
MAX_BOX_EDGE    = 80.0       # Å cap per axis — Vina's practical limit with
                             # default 0.375 Å grid; larger boxes fail with
                             # 'Cube size too large' or similar at setup

# Vina parameters (consistent with main DrugPipe pipeline)
EXHAUSTIVENESS   = 8
N_POSES          = 5
DOCKING_TIMEOUT  = 600
VINA_SEED        = 42
N_WORKERS        = 8

# Ligand preparation
MMFF_VARIANT     = "MMFF94s"
EMBED_MAX_TRIES  = 5
PH_FOR_PROT      = 7.4


# =========================================================================
# HELPERS
# =========================================================================
def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def build_session():
    """Requests session with retry for PubChem calls."""
    s = requests.Session()
    retry = Retry(
        total=5, backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    s.mount("https://", HTTPAdapter(max_retries=retry))
    return s


def fetch_cid_info(cid, session):
    """Fetch SMILES and preferred name for a PubChem CID.
       Returns dict with keys: cid, smiles, name (name may be None).

       PubChem deprecated 'CanonicalSMILES' / 'IsomericSMILES' in favour of
       'ConnectivitySMILES' / 'SMILES'. We try modern names first, then fall
       back to legacy names so this works against any PubChem version.
       For docking we prefer stereochemistry-aware SMILES."""
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid"
    out = {"cid": cid, "smiles": None, "name": None}

    # Try SMILES property names in priority order.
    # 'SMILES'              = modern, stereochemistry-aware (was IsomericSMILES)
    # 'IsomericSMILES'      = legacy stereochemistry-aware
    # 'ConnectivitySMILES'  = modern, no stereo (was CanonicalSMILES)
    # 'CanonicalSMILES'     = legacy, no stereo
    for prop in ("SMILES", "IsomericSMILES",
                 "ConnectivitySMILES", "CanonicalSMILES"):
        try:
            r = session.get(f"{base}/{cid}/property/{prop}/JSON", timeout=30)
            if not r.ok:
                continue
            props = r.json().get("PropertyTable", {}).get("Properties", [])
            if props and props[0].get(prop):
                out["smiles"] = props[0][prop]
                break
        except Exception as e:
            log(f"  ! CID {cid} {prop} fetch failed: {e}")
    if out["smiles"] is None:
        log(f"  ! CID {cid}: no SMILES returned by any property name")

    # Name: try IUPAC, fall back to first synonym
    try:
        r = session.get(f"{base}/{cid}/property/IUPACName/JSON", timeout=30)
        if r.ok:
            out["name"] = r.json()["PropertyTable"]["Properties"][0].get("IUPACName")
    except Exception:
        pass

    if not out["name"]:
        try:
            r = session.get(f"{base}/{cid}/synonyms/JSON", timeout=30)
            if r.ok:
                syns = r.json()["InformationList"]["Information"][0].get("Synonym", [])
                if syns:
                    out["name"] = syns[0]
        except Exception:
            pass

    return out


def run_cmd(cmd, timeout=None):
    """Run subprocess, return CompletedProcess."""
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def receptor_bounding_box(pdb_path):
    """Return (mn, mx, centroid) of a PDB's ATOM records."""
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except ValueError:
                    continue
    if not xs:
        raise RuntimeError(f"No atoms in {pdb_path}")
    n  = len(xs)
    mn = (min(xs), min(ys), min(zs))
    mx = (max(xs), max(ys), max(zs))
    ctr = (sum(xs)/n, sum(ys)/n, sum(zs)/n)
    return mn, mx, ctr, n


def compute_blind_box(pdb_path):
    """Compute blind-docking box centred on receptor centroid, sized to
       cover full extent + padding, capped at MAX_BOX_EDGE per axis."""
    mn, mx, ctr, n_atoms = receptor_bounding_box(pdb_path)
    extents = [mx[i] - mn[i] + 2*BOX_PADDING for i in range(3)]
    extents = [min(e, MAX_BOX_EDGE) for e in extents]
    return {
        "center_x": round(ctr[0], 3),
        "center_y": round(ctr[1], 3),
        "center_z": round(ctr[2], 3),
        "size_x":   round(extents[0], 1),
        "size_y":   round(extents[1], 1),
        "size_z":   round(extents[2], 1),
        "n_atoms":  n_atoms,
    }


# =========================================================================
# RECEPTOR PREPARATION
# =========================================================================
def prepare_receptor(pdb_path, out_pdbqt):
    """Clean waters, convert to PDBQT with Gasteiger charges at pH 7.4."""
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    # Strip HOH and non-cofactor HETATMs
    keep_hetatm = {"MG", "ZN", "CA", "FE", "MN", "CU", "NA", "K"}
    cleaned = out_pdbqt.parent / f"{out_pdbqt.stem}_cleaned.pdb"
    with open(pdb_path) as fi, open(cleaned, "w") as fo:
        for line in fi:
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname == "HOH" or resname not in keep_hetatm:
                    continue
            fo.write(line)

    cmd = [
        "obabel", str(cleaned), "-O", str(out_pdbqt),
        "-xr",                          # rigid receptor
        "-p", str(PH_FOR_PROT),
        "--partialcharge", "gasteiger",
    ]
    r = run_cmd(cmd, timeout=300)
    cleaned.unlink(missing_ok=True)
    if r.returncode != 0 or not out_pdbqt.exists():
        raise RuntimeError(f"obabel receptor failed: {r.stderr[:500]}")
    return out_pdbqt


# =========================================================================
# LIGAND PREPARATION (SMILES -> 3D SDF -> PDBQT)
# =========================================================================
def _find_lig_prep_tool():
    """Locate a working ligand-prep CLI tool. Returns (invocation_list, desc)
    or (None, reason).

    Preference order:
      1. mk_prepare_ligand.py (Meeko CLI) — works on Python 3.8 if Meeko 0.3.3
         is installed; known-good PDBQT for Vina 1.2.x
      2. prepare_ligand4.py (MGLTools/ADFR) — AutoDock reference tool
      3. prepare_ligand (ADFR variant)
    """
    # 1. Meeko CLI — preferred
    mk = shutil.which("mk_prepare_ligand.py")
    if mk:
        return ([mk], f"mk_prepare_ligand.py ({mk})")
    # 2. MGLTools
    direct = shutil.which("prepare_ligand4.py")
    if direct:
        return ([direct], f"prepare_ligand4.py ({direct})")
    direct2 = shutil.which("prepare_ligand")
    if direct2:
        return ([direct2], f"prepare_ligand ({direct2})")
    # 3. Via pythonsh (MGLTools' bundled Python)
    pythonsh = shutil.which("pythonsh")
    if pythonsh:
        for cand in [
            "/Applications/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py",
            os.path.expanduser("~/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"),
            os.path.expanduser("~/ADFRsuite-1.0/CCSBpckgs/AutoDockTools/Utilities24/prepare_ligand4.py"),
        ]:
            if os.path.isfile(cand):
                return ([pythonsh, cand], f"pythonsh + {cand}")
    return (None,
            "No ligand-prep tool found. Recommended fix for Python 3.8 (sbdd-env):\n"
            "    pip install 'meeko==0.3.3'\n"
            "  (Meeko >=0.5 requires Python 3.10+. Meeko 0.3.3 is the last\n"
            "   Python-3.8-compatible release and ships mk_prepare_ligand.py.)\n"
            "  Alternatively install MGLTools to get prepare_ligand4.py.")


def smiles_to_pdbqt(smi, cid, ligands_dir):
    """Embed + MMFF optimise + mk_prepare_ligand.py -> PDBQT.
    Returns (Path, None) on success, (None, err_msg) on failure.

    Uses Meeko's CLI (mk_prepare_ligand.py) rather than its Python API
    because the API changed shape across 0.3 -> 0.4 -> 0.6, and Meeko >= 0.5
    requires Python 3.10+ (uses PEP-585 `dict[str, Any]` syntax at import
    time, which fails on Python 3.8 with 'type object is not subscriptable').
    On Python 3.8 environments like sbdd-env, only Meeko <= 0.3.3 imports
    cleanly, but its CLI front-end is stable and produces Vina-1.2.x-
    compliant PDBQT files that Open Babel does not reliably generate."""
    sdf = ligands_dir / f"CID{cid}.sdf"
    pdbqt = ligands_dir / f"CID{cid}.pdbqt"
    if pdbqt.exists() and pdbqt.stat().st_size > 100:
        return pdbqt, None

    # --- 1. SMILES -> 3D via RDKit ---
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, "RDKit failed to parse SMILES"

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    try:
        params.maxAttempts = 200
    except (AttributeError, TypeError):
        pass
    success = -1
    for _ in range(EMBED_MAX_TRIES):
        success = AllChem.EmbedMolecule(mol, params)
        if success == 0:
            break
    if success != 0:
        return None, "ETKDG embedding failed"

    try:
        AllChem.MMFFOptimizeMolecule(mol, mmffVariant=MMFF_VARIANT, maxIters=500)
    except Exception as e:
        return None, f"MMFF optimisation failed: {e}"

    # --- 2. Write SDF (Meeko CLI accepts SDF directly) ---
    writer = Chem.SDWriter(str(sdf))
    writer.write(mol)
    writer.close()
    if not sdf.exists() or sdf.stat().st_size < 50:
        return None, "RDKit SDF write failed"

    # --- 3. SDF -> PDBQT via Meeko CLI ---
    tool = _LIG_PREP_CMD
    if tool is None:
        return None, _LIG_PREP_ERR

    cmd = list(tool) + ["-i", str(sdf), "-o", str(pdbqt)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    except subprocess.TimeoutExpired:
        return None, "mk_prepare_ligand.py timeout"
    except Exception as e:
        return None, f"mk_prepare_ligand.py exception: {e}"

    if r.returncode != 0 or not pdbqt.exists() or pdbqt.stat().st_size < 100:
        return None, (f"mk_prepare_ligand.py failed (rc={r.returncode}): "
                      f"{(r.stderr or r.stdout)[:300]}")

    # Sanity: valid flexible-ligand PDBQT must have a ROOT block
    with open(pdbqt) as f:
        if not any(ln.startswith("ROOT") for ln in f):
            return None, "PDBQT missing ROOT block"

    return pdbqt, None


# Resolve ligand-prep tool once at import time
_tool_result = _find_lig_prep_tool()
if _tool_result[0] is None:
    _LIG_PREP_CMD, _LIG_PREP_ERR = None, _tool_result[1]
    _LIG_PREP_DESC = "not available"
else:
    _LIG_PREP_CMD, _LIG_PREP_DESC = _tool_result
    _LIG_PREP_ERR = None


# =========================================================================
# DOCKING
# =========================================================================
def parse_vina_log(log_path):
    """Extract best-mode affinity from a Vina log file."""
    if not log_path.exists():
        return None
    best = None
    with open(log_path) as f:
        in_table = False
        for line in f:
            if re.match(r"\s*-+\+-+\+-+", line):
                in_table = True
                continue
            if in_table:
                m = re.match(r"\s*1\s+(-?\d+\.\d+)", line)
                if m:
                    best = float(m.group(1))
                    break
    return best


def dock_one(args):
    """Single-ligand worker. Runs Vina, returns dict with result."""
    (cid, smiles, name, degrees, receptor_pdbqt,
     ligands_dir, poses_dir, box) = args

    result = {
        "cid": cid, "name": name, "smiles": smiles,
        "degrees": ";".join(degrees),
        "vina_score": None, "status": "pending",
        "pose_file": None, "log_file": None, "error": None,
    }

    lig_pdbqt, err = smiles_to_pdbqt(smiles, cid, ligands_dir)
    if lig_pdbqt is None:
        result["status"] = "ligand_prep_failed"
        result["error"] = err
        return result

    out_pose = poses_dir / f"CID{cid}_out.pdbqt"
    log_file = poses_dir / f"CID{cid}_vina.log"

    if out_pose.exists() and log_file.exists():
        # Idempotent: reuse previous result
        score = parse_vina_log(log_file)
        if score is not None:
            result.update(status="ok_cached", vina_score=score,
                          pose_file=str(out_pose), log_file=str(log_file))
            return result

    cmd = [
        "vina",
        "--receptor", str(receptor_pdbqt),
        "--ligand",   str(lig_pdbqt),
        "--out",      str(out_pose),
        "--center_x", str(box["center_x"]),
        "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]),
        "--size_x",   str(box["size_x"]),
        "--size_y",   str(box["size_y"]),
        "--size_z",   str(box["size_z"]),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes",      str(N_POSES),
        "--seed",           str(VINA_SEED),
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=DOCKING_TIMEOUT)
        with open(log_file, "w") as f:
            f.write("=== STDOUT ===\n" + r.stdout +
                    "\n=== STDERR ===\n" + r.stderr)
        if r.returncode != 0:
            result.update(status="vina_failed",
                          error=r.stderr.strip()[-300:])
            return result
        score = parse_vina_log(log_file)
        if score is None:
            result.update(status="parse_failed",
                          error="could not parse best score from Vina log")
            return result
        result.update(status="ok", vina_score=score,
                      pose_file=str(out_pose), log_file=str(log_file))
    except subprocess.TimeoutExpired:
        result.update(status="timeout",
                      error=f"vina exceeded {DOCKING_TIMEOUT}s")
    except Exception as e:
        result.update(status="exception", error=str(e))
    return result


# =========================================================================
# CID SELECTION
# =========================================================================
def select_top_cids(xlsx_path, top_n=TOP_N_PER_DEGREE):
    """Read xlsx, take top_n from each Degree column, dedup, track degrees."""
    df = pd.read_excel(xlsx_path)
    cid_to_degrees = {}     # cid -> list of degree labels
    for col in DEGREE_COLUMNS:
        if col not in df.columns:
            log(f"  ! column '{col}' missing in {xlsx_path.name}; skipping")
            continue
        values = df[col].dropna().head(top_n).tolist()
        for v in values:
            try:
                cid = int(float(v))
            except (ValueError, TypeError):
                log(f"  ! skipping non-numeric value '{v}' in {col}")
                continue
            cid_to_degrees.setdefault(cid, []).append(col)
    return cid_to_degrees


# =========================================================================
# TARGET PIPELINE
# =========================================================================
def process_target(target_name, cfg, session):
    log("")
    log("=" * 70)
    log(f"TARGET: {target_name}")
    log("=" * 70)

    target_dir = OUTPUT_DIR / "runs" / target_name
    receptor_dir = target_dir / "receptor"
    ligands_dir  = target_dir / "ligands"
    poses_dir    = target_dir / "poses"
    for d in (receptor_dir, ligands_dir, poses_dir):
        d.mkdir(parents=True, exist_ok=True)

    # --- Receptor prep ---
    receptor_pdbqt = receptor_dir / f"{target_name}.pdbqt"
    if not receptor_pdbqt.exists():
        log(f"  preparing receptor from {cfg['pdb'].name}")
        prepare_receptor(cfg["pdb"], receptor_pdbqt)
    else:
        log(f"  receptor already prepared: {receptor_pdbqt.name}")

    # --- Blind docking box ---
    box = compute_blind_box(cfg["pdb"])
    log(f"  receptor atoms: {box['n_atoms']}")
    log(f"  blind box: centre=({box['center_x']}, {box['center_y']}, "
        f"{box['center_z']})  size=({box['size_x']}, {box['size_y']}, "
        f"{box['size_z']}) Å")

    # --- CID selection ---
    cid_to_degrees = select_top_cids(cfg["xlsx"])
    log(f"  {len(cid_to_degrees)} unique CIDs after dedup "
        f"(from top {TOP_N_PER_DEGREE} × {len(DEGREE_COLUMNS)} degrees)")

    # --- PubChem lookup ---
    log(f"  fetching SMILES + names from PubChem...")
    cid_records = []
    for i, (cid, degrees) in enumerate(cid_to_degrees.items(), 1):
        info = fetch_cid_info(cid, session)
        if not info["smiles"]:
            log(f"    [{i}/{len(cid_to_degrees)}] CID {cid}: no SMILES - SKIPPED")
            continue
        cid_records.append({
            "cid": cid,
            "name": info["name"],
            "smiles": info["smiles"],
            "degrees": degrees,
        })
        if i % 10 == 0 or i == len(cid_to_degrees):
            log(f"    [{i}/{len(cid_to_degrees)}] done")

    log(f"  {len(cid_records)} CIDs have SMILES and will be docked")

    # --- Parallel docking ---
    log(f"  docking with {N_WORKERS} workers "
        f"(exhaustiveness={EXHAUSTIVENESS}, poses={N_POSES}, seed={VINA_SEED})")

    args_list = [
        (rec["cid"], rec["smiles"], rec["name"], rec["degrees"],
         receptor_pdbqt, ligands_dir, poses_dir, box)
        for rec in cid_records
    ]

    results = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = {ex.submit(dock_one, a): a[0] for a in args_list}
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            results.append(res)
            score_str = (f"{res['vina_score']:.3f}"
                         if res["vina_score"] is not None else "---")
            log(f"    [{i}/{len(futures)}] CID {res['cid']}: "
                f"{res['status']:20s} score={score_str}")
    elapsed = (time.time() - t0) / 60
    log(f"  docking complete in {elapsed:.1f} min")

    # --- Build results DataFrame ---
    df_out = pd.DataFrame(results)
    df_out["target"] = target_name
    df_out["box_center"] = (f"({box['center_x']}, {box['center_y']}, "
                            f"{box['center_z']})")
    df_out["box_size"] = (f"({box['size_x']}, {box['size_y']}, "
                          f"{box['size_z']})")
    cols = ["target", "cid", "name", "degrees", "smiles",
            "vina_score", "status", "box_center", "box_size",
            "pose_file", "log_file", "error"]
    df_out = df_out[cols].sort_values(
        by="vina_score", ascending=True, na_position="last"
    ).reset_index(drop=True)
    df_out.insert(0, "rank", range(1, len(df_out) + 1))

    out_xlsx = target_dir / f"{target_name}_docking_results.xlsx"
    df_out.to_excel(out_xlsx, index=False)
    log(f"  results: {out_xlsx}")
    return df_out


# =========================================================================
# MAIN
# =========================================================================
def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Sanity: required binaries
    for tool in ("vina", "obabel"):
        if shutil.which(tool) is None:
            sys.exit(f"ERROR: '{tool}' not on PATH. "
                     f"Activate your conda env (sbdd-env) first.")

    # Sanity: inputs exist
    for name, cfg in TARGETS.items():
        for key in ("pdb", "xlsx"):
            if not cfg[key].exists():
                sys.exit(f"ERROR: missing input for {name} ({key}): {cfg[key]}")

    # Sanity: ligand-prep self-test so we catch setup issues BEFORE docking
    # 18 ligands. Tests full SMILES->SDF->PDBQT pipeline on ethanol.
    log("#" * 70)
    log("BLIND DOCKING — top CIDs per degree for BMP1 and MMP9")
    log(f"  input:  {INPUT_DIR}")
    log(f"  output: {OUTPUT_DIR}")
    log("#" * 70)
    if _LIG_PREP_CMD is None:
        log(f"  !! ligand-prep tool not found:")
        for line in _LIG_PREP_ERR.splitlines():
            log(f"     {line}")
        sys.exit(1)
    log(f"  ligand-prep tool: {_LIG_PREP_DESC}")
    import tempfile
    test_dir = Path(tempfile.mkdtemp(prefix="prep_test_"))
    p, err = smiles_to_pdbqt("CCO", "TEST_ETHANOL", test_dir)
    if p is None:
        log(f"  !! ligand-prep self-test FAILED: {err}")
        log(f"     aborting before we waste hours on guaranteed failures")
        sys.exit(1)
    log(f"  ligand-prep self-test OK (ethanol PDBQT: {p.stat().st_size} bytes)")
    shutil.rmtree(test_dir, ignore_errors=True)

    session = build_session()

    all_results = {}
    for target_name, cfg in TARGETS.items():
        try:
            all_results[target_name] = process_target(target_name, cfg, session)
        except Exception as e:
            log(f"!! {target_name} FAILED: {e}")
            raise

    # Combined workbook
    combined_xlsx = OUTPUT_DIR / "all_targets_docking_summary.xlsx"
    with pd.ExcelWriter(combined_xlsx, engine="openpyxl") as writer:
        for target_name, df in all_results.items():
            df.to_excel(writer, sheet_name=target_name, index=False)
    log("")
    log("#" * 70)
    log(f"ALL DONE.  combined summary: {combined_xlsx}")
    log("#" * 70)


if __name__ == "__main__":
    main()
