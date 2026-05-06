"""
08_enrich_drugs.py — Enrich top-N candidates per target with external drug info.

For each target, reads the clean shortlist from analysis/<target>/<target>_clean_shortlist.csv
(produced by 07_analyze.py), then for each drug queries:

  - PubChem:   CID, synonyms, PubChem URL (via InChIKey lookup)
  - ChEMBL:    CHEMBL ID, known targets, approved indications (diseases)
  - UniChem:   DrugBank ID → ChEMBL ID crosswalk (more reliable than PubChem xref)

Output (per target):
  analysis/<target>/<target>_enriched_top<N>.csv

Also writes a cross-target master file:
  analysis/all_targets_enriched.csv

Run:
  python 08_enrich_drugs.py                   # process all 9 targets, top 20 each
  python 08_enrich_drugs.py IGFBP7 BGN        # specific targets
  TOP_N=50 python 08_enrich_drugs.py          # override top-N per target
"""

import os
import sys
import time
import pandas as pd
import requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

RDLogger.DisableLog("rdApp.*")

# ==============================================================
# CONFIG
# ==============================================================
BASE_DIR      = Path(os.environ.get("DRUGPIPE_BASE_DIR", Path(__file__).resolve().parent.parent))
ANALYSIS_DIR  = BASE_DIR / "analysis"

TOP_N         = int(os.environ.get("TOP_N", "20"))        # per-target depth
MAX_WORKERS   = 6                                          # parallel drugs
CHEMBL_WORKERS = 3                                         # inner parallel for ChEMBL

# Which targets to process (if no CLI args given)
ALL_TARGETS = [
    "IGFBP7", "BGN", "BMP1",
    "INHBA_with_propeptide", "INHBA_without_propeptide",
    "MMP9", "SPARC_monomer", "TANGL", "TNC",
]


# ==============================================================
# SESSION WITH RETRIES (from your reference script)
# ==============================================================
def build_session():
    s = requests.Session()
    retry = Retry(
        total=5, connect=5, read=5,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=20, pool_maxsize=20)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers.update({"User-Agent": "Mozilla/5.0 DrugPipe-Pipeline"})
    return s

session = build_session()


# ==============================================================
# API HELPERS
# ==============================================================
def safe_get_json(url, timeout=20):
    try:
        r = session.get(url, timeout=timeout)
        if r.status_code == 200:
            return r.json()
        if r.status_code == 404:
            return {}                       # no data for this ID
        return None                         # retryable failure
    except Exception:
        return None


def smiles_to_inchikey(smi):
    """Convert SMILES to InChIKey for reliable PubChem lookup."""
    if not isinstance(smi, str):
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    try:
        return inchi.MolToInchiKey(mol)
    except Exception:
        return None


def get_pubchem_cid_from_inchikey(inchikey):
    """InChIKey → PubChem CID."""
    if not inchikey:
        return None
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
    data = safe_get_json(url)
    if not data:
        return None
    try:
        cids = data.get("IdentifierList", {}).get("CID", [])
        return int(cids[0]) if cids else None
    except Exception:
        return None


def get_pubchem_synonyms(cid):
    if not cid:
        return []
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
    data = safe_get_json(url)
    if not data:
        return []
    try:
        info = data.get("InformationList", {}).get("Information", [])
        if info and isinstance(info, list):
            return info[0].get("Synonym", []) or []
    except Exception:
        pass
    return []


def drugbank_to_chembl(drugbank_id):
    """UniChem crosswalk: DrugBank ID → ChEMBL ID.
    UniChem src IDs: DrugBank=2, ChEMBL=1
    """
    if not drugbank_id:
        return []
    url = f"https://www.ebi.ac.uk/unichem/rest/src_compound_id/{drugbank_id}/2/1"
    data = safe_get_json(url, timeout=15)
    if not data or not isinstance(data, list):
        return []
    return [item.get("src_compound_id") for item in data if item.get("src_compound_id")]


def get_chembl_targets(chembl_id):
    url = (f"https://www.ebi.ac.uk/chembl/api/data/activity.json"
           f"?molecule_chembl_id={chembl_id}&limit=100")
    data = safe_get_json(url, timeout=30)
    targets = set()
    if data:
        for act in data.get("activities", []) or []:
            tname = act.get("target_pref_name")
            if tname:
                targets.add(tname)
    return targets


def get_chembl_diseases(chembl_id):
    url = (f"https://www.ebi.ac.uk/chembl/api/data/drug_indication.json"
           f"?molecule_chembl_id={chembl_id}")
    data = safe_get_json(url, timeout=30)
    diseases = set()
    max_phase = None
    if data:
        for ind in data.get("drug_indications", []) or []:
            efo = ind.get("efo_term")
            if efo:
                diseases.add(efo)
            phase = ind.get("max_phase_for_ind")
            if phase is not None:
                try:
                    p = float(phase)
                    if max_phase is None or p > max_phase:
                        max_phase = p
                except (TypeError, ValueError):
                    pass
    return diseases, max_phase


def get_chembl_molecule_info(chembl_id):
    """Approval status, oral/parenteral, known MoA."""
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    data = safe_get_json(url, timeout=20)
    if not data:
        return {}
    try:
        return {
            "max_phase": data.get("max_phase"),
            "oral": data.get("oral"),
            "parenteral": data.get("parenteral"),
            "first_approval": data.get("first_approval"),
            "indication_class": data.get("indication_class"),
            "mol_type": data.get("molecule_type"),
        }
    except Exception:
        return {}


# ==============================================================
# ENRICH ONE DRUG
# ==============================================================
def enrich_drug(row):
    """Takes a row from the shortlist (dict-like), returns enriched dict."""
    drug_id    = row.get("Drug id") or row.get("drug_id")
    drug_name  = row.get("Drug name") or row.get("drug_name")
    smiles     = row.get("smiles")

    result = {
        "Drug id":       drug_id,
        "Drug name":     drug_name,
        "smiles":        smiles,
        "vina_score":    row.get("vina_score"),
        "pocket_num":    row.get("pocket_num"),
        "tanimoto_score": row.get("tanimoto_score"),
        "mw":            row.get("mw"),
        "lipinski_ok":   row.get("lipinski_ok"),
    }

    # DrugBank → ChEMBL via UniChem
    chembl_ids = drugbank_to_chembl(drug_id) if drug_id else []

    # SMILES → InChIKey → PubChem CID
    ikey = smiles_to_inchikey(smiles)
    cid  = get_pubchem_cid_from_inchikey(ikey) if ikey else None
    synonyms = get_pubchem_synonyms(cid) if cid else []

    # ChEMBL data (targets, diseases, approval)
    all_targets, all_diseases = set(), set()
    max_phase, indication_class, first_approval = None, None, None

    if chembl_ids:
        # Run ChEMBL calls in parallel per molecule
        with ThreadPoolExecutor(max_workers=min(CHEMBL_WORKERS, len(chembl_ids))) as pool:
            tgt_f = [pool.submit(get_chembl_targets, cid_) for cid_ in chembl_ids]
            dis_f = [pool.submit(get_chembl_diseases, cid_) for cid_ in chembl_ids]
            mol_f = [pool.submit(get_chembl_molecule_info, cid_) for cid_ in chembl_ids]
            for f in as_completed(tgt_f):
                all_targets.update(f.result())
            for f in as_completed(dis_f):
                d, phase = f.result()
                all_diseases.update(d)
                if phase is not None and (max_phase is None or phase > max_phase):
                    max_phase = phase
            for f in as_completed(mol_f):
                info = f.result()
                if info.get("indication_class") and not indication_class:
                    indication_class = info["indication_class"]
                if info.get("first_approval") and not first_approval:
                    first_approval = info["first_approval"]
                if info.get("max_phase") is not None:
                    try:
                        mp = float(info["max_phase"])
                        if max_phase is None or mp > max_phase:
                            max_phase = mp
                    except (TypeError, ValueError):
                        pass

    result.update({
        "PubChem_CID":    cid,
        "PubChem_URL":    f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" if cid else None,
        "DrugBank_URL":   f"https://go.drugbank.com/drugs/{drug_id}" if drug_id else None,
        "Synonyms":       "; ".join(synonyms[:5]) if synonyms else None,
        "ChEMBL_IDs":     "; ".join(chembl_ids) if chembl_ids else None,
        "Known_Targets":  "; ".join(sorted(all_targets)) if all_targets else None,
        "N_Known_Targets": len(all_targets),
        "Indications":    "; ".join(sorted(all_diseases)) if all_diseases else None,
        "N_Indications":  len(all_diseases),
        "Max_Clinical_Phase": max_phase,
        "Indication_Class":   indication_class,
        "First_Approval":     first_approval,
    })
    return result


# ==============================================================
# PROCESS ONE TARGET
# ==============================================================
def process_target(target):
    print(f"\n{'='*70}")
    print(f"Enriching {target} — top {TOP_N} candidates")
    print(f"{'='*70}")

    shortlist_csv = ANALYSIS_DIR / target / f"{target}_clean_shortlist.csv"
    if not shortlist_csv.exists():
        # Fallback: read directly from per-pocket final.csv files
        print(f"  clean shortlist missing, falling back to pocket final.csv files")
        parts = []
        for fc in (BASE_DIR / "runs" / target).glob(f"pocket*/{target}_pocket*_final.csv"):
            df = pd.read_csv(fc)
            parts.append(df)
        if not parts:
            print(f"  no data for {target}, skipping")
            return None
        shortlist = pd.concat(parts, ignore_index=True)
        shortlist = shortlist.sort_values("vina_score").drop_duplicates("Drug id").head(TOP_N)
    else:
        shortlist = pd.read_csv(shortlist_csv).head(TOP_N)

    print(f"  {len(shortlist)} drugs to enrich")
    rows = shortlist.to_dict(orient="records")

    results = []
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(enrich_drug, r): r for r in rows}
        for i, f in enumerate(as_completed(futures), 1):
            r = futures[f]
            try:
                enriched = f.result()
                results.append(enriched)
                name = enriched.get("Drug name", "?")
                nt = enriched.get("N_Known_Targets", 0)
                ni = enriched.get("N_Indications", 0)
                phase = enriched.get("Max_Clinical_Phase", "?")
                print(f"  [{i:2d}/{len(rows)}] {name[:40]:<40}  "
                      f"targets={nt:>3}  indications={ni:>3}  phase={phase}")
            except Exception as e:
                print(f"  [{i:2d}/{len(rows)}] FAILED {r.get('Drug name')}: {e}")
                results.append({"Drug id": r.get("Drug id"),
                                "Drug name": r.get("Drug name"),
                                "enrichment_error": str(e)})

    # Preserve shortlist order (by vina_score ascending)
    df = pd.DataFrame(results)
    order_map = {r["Drug id"]: i for i, r in enumerate(rows) if r.get("Drug id")}
    df["_order"] = df["Drug id"].map(order_map).fillna(999)
    df = df.sort_values("_order").drop(columns="_order").reset_index(drop=True)

    out = ANALYSIS_DIR / target / f"{target}_enriched_top{TOP_N}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"  wrote {out.relative_to(BASE_DIR)}  ({time.time()-t0:.0f}s)")

    df["target"] = target
    return df


# ==============================================================
# MAIN
# ==============================================================
def main():
    targets = sys.argv[1:] if len(sys.argv) > 1 else ALL_TARGETS

    all_dfs = []
    for t in targets:
        df = process_target(t)
        if df is not None:
            all_dfs.append(df)

    if all_dfs:
        master = pd.concat(all_dfs, ignore_index=True)
        cols = ["target"] + [c for c in master.columns if c != "target"]
        master = master[cols]
        out = ANALYSIS_DIR / "all_targets_enriched.csv"
        master.to_csv(out, index=False)
        print(f"\n{'='*70}")
        print(f"Master enriched table: {out.relative_to(BASE_DIR)}")
        print(f"  {len(master)} rows across {len(all_dfs)} targets")

        # Quick headline: which targets have approved-drug hits (max_phase >= 4)?
        if "Max_Clinical_Phase" in master.columns:
            approved = master[master["Max_Clinical_Phase"] >= 4]
            if len(approved) > 0:
                print(f"\nApproved drugs appearing as top hits ({len(approved)} rows):")
                cols_show = ["target", "Drug name", "vina_score",
                             "Max_Clinical_Phase", "Indications"]
                cols_show = [c for c in cols_show if c in approved.columns]
                print(approved[cols_show].to_string(index=False))


if __name__ == "__main__":
    main()