# Target2DeNovoDrug — QSAR-Based Virtual Screening Notebook

Locally adapted Jupyter notebook implementation of **Target2DeNovoDrug** (Madaj et al. 2022, *Journal of Biomolecular Structure and Dynamics*) for ligand-based virtual screening. Built and used during the BSc dissertation *"In Silico Drug Discovery in Atherosclerosis"* (Yikai Yang, 2026).

This notebook is the route taken when a target **does have usable bioassay data** in PubChem — i.e. enough known active compounds with measured activity values to train a QSAR model. For the structure-based route used when bioassay data are missing, see `DrugPipe_README.md`.

The workflow has two parts:

```
drug2target_april_2026_optimized_faster.ipynb   ← QSAR ranking (main notebook)
dock_top_drugs.py                               ← blind docking of top hits
build_comparison_table.py                       ← QSAR vs Vina comparison table
```

The notebook handles steps 1–10 (target selection through QSAR ranking). The two helper scripts run after the notebook to validate the QSAR ranking by docking the top-10 predicted compounds and comparing the rankings.

---

## What the notebook does

For one protein target at a time, the notebook:

1. Looks up whether the target has bioassay data in PubChem (manual web check)
2. Runs a BLAST search of the protein sequence against SwissProt to confirm UniProt identity
3. Downloads PubChem bioactivity data for the target's UniProt accession
4. Aggregates compounds across multiple assays, deduplicates by CID, and computes a single activity value per compound
5. Expands the training set by retrieving structurally similar compounds from PubChem (with caching, retries, and parallel fetching)
6. Trains a **QSAR regression model** (linear regression with polynomial features at degrees 1, 2, and 3) on molecular descriptors against measured activity
7. Selects the best-performing model and predicts activity for a wider set of PubChem compounds
8. Exports the top-20 ranked compounds, enriched with names, ChEMBL IDs, known targets, and indication info
9. Optionally runs **AutoDock Vina** docking on top hits against a user-provided processed PDB

Final output: a ranked, annotated CSV of QSAR-predicted candidates per target, plus optional docking poses for the top hits.

---

## Notebook flow at a glance

```
┌────────────────────────────────────────────────┐
│  Cell 0: imports                               │
└────────────────────────────────────────────────┘
                      │
                      ▼
┌────────────────────────────────────────────────┐
│  Cell 1–2: PubChem bioassay availability check │
│  (manual: open the URL in a browser)           │
└────────────────────────────────────────────────┘
                      │
                      ▼
┌────────────────────────────────────────────────┐
│  Cell 3: BLAST query against SwissProt         │
│  → confirms UniProt accession for your protein │
└────────────────────────────────────────────────┘
                      │
                      ▼
┌────────────────────────────────────────────────┐
│  Cell 4–5: download bioassay data from PubChem │
│  → save summary CSV manually first             │
└────────────────────────────────────────────────┘
                      │
                      ▼
┌────────────────────────────────────────────────┐
│  Cell 6–9: build training set                  │
│  (deduplicate, normalise activity column)      │
└────────────────────────────────────────────────┘
                      │
                      ▼
┌────────────────────────────────────────────────┐
│  Cell 7: expand training set with similar CIDs │
│  (parallelised, cached)                        │
└────────────────────────────────────────────────┘
                      │
                      ▼
┌────────────────────────────────────────────────┐
│  Cell 10–13: data prep + sanity checks         │
└────────────────────────────────────────────────┘
                      │
                      ▼
┌────────────────────────────────────────────────┐
│  Cell 14: QSAR model training                  │
│  (degrees 1/2/3, 10-fold cross-validation)     │
└────────────────────────────────────────────────┘
                      │
                      ▼
┌────────────────────────────────────────────────┐
│  Cell 15–18: predict activity for wider set    │
│  → export top-N ranked compounds with metadata │
└────────────────────────────────────────────────┘
                      │
                      ▼
┌────────────────────────────────────────────────┐
│  Cell 20–21: view ranked compound table        │
└────────────────────────────────────────────────┘
                      │
                      ▼
┌────────────────────────────────────────────────┐
│  Cell 22–24: optional AutoDock Vina docking    │
│  on top hits against your local PDB            │
└────────────────────────────────────────────────┘
```

---

## When this notebook is the right tool

Target2DeNovoDrug needs **prior knowledge of bioactive compounds** for your target. If your target doesn't have at least a few hundred compounds with measured activity values in PubChem, the QSAR model will not be trainable in any useful way and you should switch to DrugPipe instead.

In the dissertation:

| Target | Bioassays found | QSAR usable? | Route used |
|---|---|---|---|
| IGFBP7 | 0 | No | DrugPipe |
| TAGLN | 0 | No | DrugPipe |
| SPARC | 0 | No | DrugPipe |
| TNC | 0 | No | DrugPipe |
| BGN | 0 | No | DrugPipe |
| INHBA | 3 (no usable values) | No | DrugPipe |
| **BMP1** | 829 unique compounds | **Yes** | **Target2DeNovoDrug** |
| **MMP9** | 1,060 unique compounds | **Yes** | **Target2DeNovoDrug** |

---

## Folder layout

There is no fixed folder structure for this notebook beyond what you set in the cell variables. A reasonable layout:

```
Target2DeNovoDrug_runs/
├── drug2target_april_2026_optimized_faster.ipynb
│
├── bioassays/                            ← downloaded summary CSVs (one per target)
│   ├── BMP1_pubchem_protacxn_P13497_bioassay.csv
│   └── MMP9_pubchem_protacxn_P14780_bioassay.csv
│
├── PDB/                                  ← processed PDB files for docking
│   ├── BMP1_energy_minimized.pdb
│   └── MMP9_energy_minimized.pdb
│
├── results/                              ← final ranked compound CSVs
│   ├── BMP1_1699_QSAR_CID_Info.csv
│   └── MMP9_4543_QSAR_CID_Info.csv
│
├── docking/                              ← Vina output for top hits
│   └── tmp_<runid>/
│       ├── <TARGET>.pdbqt
│       ├── <CID>.pdbqt
│       └── <CID>_out.pdbqt
│
└── cid_cache.json                        ← auto-created PubChem similarity cache
```

The notebook itself uses absolute paths in several cells (originally `/Users/yikai/Desktop/...`). **Search the notebook for `/Users/yikai/` and replace with your own paths before running.**

---

## Software requirements

The notebook was developed on **macOS** with a CPU-only MacBook Air. AutoDock Vina docking (Cells 23–24) uses MGLTools, which has stricter installation requirements than Open Babel.

### Python environment

```bash
conda create -n qsar-env python=3.9
conda activate qsar-env
pip install jupyter pandas numpy scikit-learn matplotlib seaborn \
            biopython biopandas requests selenium pymol-open-source
```

### Command-line tools

| Tool | Purpose | Install |
|---|---|---|
| AutoDock Vina | Docking | `brew install autodock-vina` |
| Open Babel | Ligand 3D generation | `brew install open-babel` |
| MGLTools 1.5.6 | Receptor/ligand prep | Download from [CCSB](http://mgltools.scripps.edu/downloads) |

The MGLTools paths in Cell 23 hardcode the install directory:

```python
prepare_protein_path = 'mgltools_x86_64Linux2_1.5.6/bin/pythonsh ' \
                       'MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r '
prepare_ligand_path  = 'mgltools_x86_64Linux2_1.5.6/bin/pythonsh ' \
                       'MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -A bonds_hydrogens -U nphs_lps -l'
```

Edit these to match your install location.

### Selenium driver (optional)

Cell 0 imports Selenium with Firefox. This is only used if you run the optional automated PubChem screenshot helpers. If you don't need those, you can ignore Selenium errors at import.

---

## Step-by-step protocol

### Step 1 — Check whether your target has bioassay data

Open this URL in a browser, replacing `IGFBP7` with your gene symbol:

```
https://pubchem.ncbi.nlm.nih.gov/#query=IGFBP7&tab=bioassay
```

If there are no bioassay results, **stop here and use DrugPipe instead**. If there are several bioassays with measured `Standard Value` or `PUBCHEM_ACTIVITY_SCORE` columns, continue.

Download the **summary CSV** of bioassays for your target by clicking *Download* on the bioassay tab. Save it somewhere you can point the notebook at, e.g.:

```
bioassays/BMP1_pubchem_protacxn_P13497_bioassay.csv
```

### Step 2 — Run Cell 0 (imports)

Make sure all imports succeed. If `selenium` or `pymol` fail and you don't plan to use them, comment out those lines.

### Step 3 — Run Cell 3 (BLAST)

Paste your protein sequence in FASTA format when prompted (no header line). The cell:

- Sends the sequence to NCBI BLAST against SwissProt
- Saves the result to `my_blast_result.xml`
- Extracts the top 3 hits and reports them

This step confirms that you have the right UniProt accession before downloading bioassay data. Note the accession code (e.g. `P13497` for BMP1).

### Step 4 — Run Cell 4 (download bioactivity by accession)

Edit the `check_id` variable to your UniProt accession:

```python
check_id = 'P13497'   # BMP1
```

The cell builds a PubChem SDQ download URL for this accession and saves the result to `downloaded1_example.csv`.

### Step 5 — Run Cell 5 (load + aggregate bioassay data)

Edit the `summary_file` path to point at the summary CSV you downloaded in Step 1:

```python
summary_file = "/path/to/BMP1_pubchem_protacxn_P13497_bioassay.csv"
```

The cell:

- Reads the summary CSV
- Limits to the first 30 assays for testing (edit the slice if you want all)
- Downloads each assay's compound list via the PubChem REST API
- Picks the activity column (`Standard Value` if available, else `PUBCHEM_ACTIVITY_SCORE`)
- Deduplicates by `PUBCHEM_CID`, taking the mean activity per compound
- Stores the result in `data`

If you see "Total rows before deduplication: 0", your assays don't have usable activity columns — switch to DrugPipe.

### Step 6 — Run Cell 7 (expand training set with similar compounds)

This is the longest step. For each compound in `data`, it fetches **up to 30 structurally similar CIDs** from PubChem, gets their molecular properties, and assigns the parent compound's activity value to them. This expands a training set of ~1,000 compounds to ~4,000–10,000.

Tunable parameters at the top of the cell:

```python
original_cid_limit   = 10_000_000   # max original CIDs to use
target_ok            = 1_000_000    # stop after this many accepted rows
cid_keep             = 30           # max similar CIDs per original CID
BATCH_SIZE           = 100          # CIDs per batch API call (PubChem max)
MAX_WORKERS_SIMILAR  = 10           # parallel threads for similar-CID fetching
MAX_WORKERS_PROPS    = 8            # parallel threads for property batch fetching
MAX_WORKERS_ASSAY    = 25           # parallel threads for assay fetching
```

Results are cached in `cid_cache.json` next to the notebook. Re-running this cell after a crash resumes from the cache (this is essential — a cold run takes 30+ min for 1,000 compounds).

The cell ends with `x_data` (descriptors) and `y_data` (activity values) ready for QSAR.

### Step 7 — Run Cells 10, 12, 13 (prep + sanity checks)

Quick prints to confirm shapes and that the activity column has reasonable variation. If `train_file["y"]` has only one unique value, the QSAR cannot learn anything.

### Step 8 — Run Cell 14 (QSAR model training)

The main modelling cell. It:

- Splits `(x_, y_)` into 90% train / 10% test
- For each polynomial degree (1, 2, 3):
  - Builds combinations of features
  - Cross-validates with 10-fold KFold
  - Fits `Pipeline(StandardScaler → PolynomialFeatures(degree) → LinearRegression)`
  - Records best test R²
- Uses the best degree to generate predictions for `final_data_frame` (the wider compound set)
- Stores results with predicted activity per CID

In the dissertation:

| Target | Best degree | Test R² | Training size after expansion |
|---|---|---|---|
| BMP1 | 2 | 0.1530 | 1,699 |
| MMP9 | 3 | 0.2631 | 4,543 |

These R² values are low. Treat predictions as exploratory ranking, not as robust binding affinities.

### Step 9 — Run Cell 15 / 17 (export top hits with metadata)

Edit `output_csv` to your desired path:

```python
output_csv = "/path/to/results/BMP1_1699_QSAR_CID_Info.csv"
top_n      = 20
degree_col = "Degree 2"   # or "Degree 3" depending on best model
```

The cell takes the top-N predicted compounds and enriches them in parallel with:

- PubChem synonyms / common names
- ChEMBL IDs (via PubChem cross-reference)
- Known ChEMBL targets
- ChEMBL indications + max clinical phase

Output is a CSV ready to inspect.

Cell 18 is an alternative version that re-enriches an already-exported Excel file without re-running the QSAR step.

### Step 10 — View results (Cells 20–21)

```python
QSAR_df.head(20)
```

The top-ranked CIDs are your candidate compounds for this target.

### Step 11 — Batch-dock the top hits with `dock_top_drugs.py`

Cells 23–24 in the notebook can dock individual CIDs but they're awkward for batch work. `dock_top_drugs.py` is the production version: it reads the QSAR Excel output directly, fetches SMILES from PubChem, prepares ligands and receptor automatically, and runs AutoDock Vina in parallel.

**What it does**

For each target (BMP1 and MMP9 by default):

1. Reads `<TARGET>_drug_screening_results.xlsx` and pulls the top-10 CIDs from each of the three Degree columns
2. Deduplicates across degrees (a CID predicted by all three only gets docked once)
3. Fetches SMILES and IUPAC names for each CID from PubChem (with retries and fallback property names — handles both modern `SMILES` and legacy `IsomericSMILES` PubChem schemas)
4. Cleans the receptor PDB (strips waters, keeps common metal cofactors), converts to PDBQT via Open Babel
5. Computes a **blind-docking box** sized to cover the full protein bounding box plus 5 Å padding, capped at 80 Å per axis
6. Converts each ligand SMILES → 3D SDF → PDBQT (RDKit + MMFF94s minimisation, Open Babel for charges)
7. Runs AutoDock Vina in parallel (8 workers by default) with exhaustiveness 8, 5 poses per ligand, fixed seed 42 for reproducibility
8. Writes per-target results and a combined summary

Edit the paths at the top:

```python
INPUT_DIR  = Path("/path/to/your/Drug_screening results")
OUTPUT_DIR = Path("/path/to/your/Drug_screening results/docking_out")
```

Run it:

```bash
conda activate sbdd-env
python dock_top_drugs.py
```

Tunable parameters at the top:

| Parameter | Default | What it controls |
|---|---|---|
| `TOP_N_PER_DEGREE` | 10 | How many CIDs to take from each Degree column |
| `BOX_PADDING` | 5.0 Å | Padding around receptor bounding box |
| `MAX_BOX_EDGE` | 80.0 Å | Vina's practical box-size cap (larger boxes fail) |
| `EXHAUSTIVENESS` | 8 | Vina default; lower = faster, higher = more thorough |
| `N_POSES` | 5 | Poses Vina reports per ligand |
| `DOCKING_TIMEOUT` | 600 s | Kills runaway dockings |
| `VINA_SEED` | 42 | Fixed seed for reproducibility |
| `N_WORKERS` | 8 | Parallel docking jobs |

**Outputs (per target)**:

```
docking_out/runs/<TARGET>/
├── receptor/<TARGET>.pdbqt          ← prepared receptor
├── ligands/
│   ├── CID<XXXX>.sdf                ← 3D structures from RDKit
│   └── CID<XXXX>.pdbqt              ← Vina-ready ligands
├── poses/
│   ├── CID<XXXX>_out.pdbqt          ← docked poses (open in PyMOL)
│   └── CID<XXXX>_vina.log           ← per-ligand Vina output
└── <TARGET>_docking_results.xlsx    ← per-target ranked table
```

Plus a combined cross-target file:

```
docking_out/all_targets_docking_summary.xlsx     ← one sheet per target
```

Each result row contains: `rank, target, cid, name, degrees, smiles, vina_score, status, box_center, box_size, pose_file, log_file, error`.

The script runs a **ligand-prep self-test on ethanol** before starting the real docking. If the self-test fails, the script aborts immediately rather than waste hours on guaranteed failures — useful when the conda environment isn't set up correctly.

### Step 12 — Compare QSAR ranking vs Vina ranking with `build_comparison_table.py`

This script answers the question "did docking confirm the QSAR ranking?" — it cross-references the QSAR top-10 against the Vina docking scores from Step 11 and produces a side-by-side comparison.

**What it does**

For each target and each Degree column (1, 2, 3):

1. Reads the QSAR top-10 CIDs in their original input-file order (rank 1 = best QSAR prediction)
2. Reads the Vina docking results from Step 11
3. Re-sorts the same CIDs by Vina score (most negative first); CIDs that failed to dock sink to the bottom
4. Builds a wide table with three columns per degree: `DegN_QSAR`, `DegN_Vina`, `DegN_Vina_Score`
5. Builds a long-format `rank_changes` sheet showing how each CID's rank shifted between the two methods

Edit the paths at the top to match your `dock_top_drugs.py` output:

```python
INPUT_DIR = Path("/path/to/your/Drug_screening results")
DOCK_DIR  = Path("/path/to/your/Drug_screening results/docking_out")
```

Run:

```bash
python build_comparison_table.py
```

**Outputs**:

```
docking_out/comparison/
├── BMP1_qsar_vs_vina.xlsx
└── MMP9_qsar_vs_vina.xlsx
```

Each workbook has two sheets:

- **qsar_vs_vina** — 10 rows × 9 columns. For each Degree column, you see the QSAR-ranked CID, the same CIDs re-sorted by Vina score, and the actual score. CIDs without a usable Vina score show one of three labels in the Score column:
  - `NO_SMILES` — couldn't generate a 3D structure
  - `FAILED` — Vina ran but didn't produce a parseable score
  - `NOT_DOCKED` — CID wasn't found in the docking results
- **rank_changes** — long-format table with columns `degree, cid, qsar_rank, vina_rank, rank_diff`. A positive `rank_diff` means the compound climbed when docking re-scored it; negative means it fell.

This is the table that goes into the dissertation comparison figures. Use it to identify which top QSAR predictions were also strong Vina hits (rows where `DegN_QSAR == DegN_Vina`) and which weren't.

### Step 11 (alternative) — Quick-and-dirty single-CID docking from inside the notebook (Cells 23–24)

Cell 23 defines three classes: `ProteinPreparer`, `LigandPreparer`, `VinaDocker`. They handle:

- Cleaning the receptor PDB and converting to PDBQT (via MGLTools)
- Downloading the ligand SDF from PubChem, generating 3D coords (Open Babel), force-field minimisation (GAFF), conversion to PDBQT
- Running Vina, merging poses, generating a PLIP interaction report

Cell 24 is the runner. Edit:

```python
top_cid = [146478919]   # CIDs you want to dock (from your QSAR results)
local_pdb_path = "/path/to/your/processed/protein.pdb"
```

The cell creates a temp directory, prepares both files, and docks every CID in `top_cid` against the receptor. Outputs land in the temp directory (not deleted automatically — copy them out if you want to keep them).

---

## Where to find your final results

**QSAR-ranked compounds with metadata** (the main deliverable from the notebook):

```
<your output_csv path>     e.g. results/BMP1_1699_QSAR_CID_Info.csv
```

This CSV has one row per top compound, with predicted activity, PubChem CID, common name, ChEMBL ID, known targets, indications, and max clinical phase.

**Per-target docking results** (from `dock_top_drugs.py`):

```
docking_out/runs/<TARGET>/<TARGET>_docking_results.xlsx
```

Ranked table of every docked CID with Vina score, status, pose file paths.

**Cross-target docking summary**:

```
docking_out/all_targets_docking_summary.xlsx
```

One sheet per target — the headline cross-protein comparison.

**QSAR vs Vina ranking comparison** (from `build_comparison_table.py`):

```
docking_out/comparison/<TARGET>_qsar_vs_vina.xlsx
```

The side-by-side comparison used in the dissertation comparison figures.

**Docking poses for individual hits** (open in PyMOL alongside the receptor):

```
docking_out/runs/<TARGET>/poses/CID<XXXX>_out.pdbqt
```

Or if you used the notebook's Cell 24 instead, the temp directory printed at the end of that cell.

---

## Common pitfalls

**`Found 0 target assay(s)`** in Cell 5 — your summary CSV's `aidname` column doesn't match. The original cell uses `str.contains("")` (matches everything), so this is rare. If you edited it to filter by gene name and your assays use synonyms, broaden the filter.

**Cell 7 hangs / times out on PubChem** — PubChem rate-limits at ~5 req/s. Lower `MAX_WORKERS_SIMILAR` from 10 to 4 if you see frequent 503s. The session has retry/backoff but extreme rate limiting still kills throughput.

**`cid_cache.json` corrupted** — the cell handles this and starts fresh, but if you killed the kernel mid-write, manually delete the file.

**QSAR R² is below zero** — the model is worse than predicting the mean. This means either: (a) your descriptors don't carry signal for this target, (b) the activity values you scraped mix incompatible assays (e.g. IC50 in nM and % inhibition), or (c) your training set is too small. Look at `train_file.corr(numeric_only=True)["y"]` from Cell 13 — if no descriptor correlates above ~0.1, abandon QSAR for this target.

**MGLTools `prepare_receptor4.py` fails** — MGLTools is sensitive to PDB formatting. Strip waters, ligands, and alternate conformers from your PDB first (PyMOL: `remove resn HOH; remove not alt ""+A`). Use `pdb4amber` if available.

**Vina segfaults on Cell 24** — usually means the receptor PDBQT or the ligand PDBQT is malformed. Open both in PyMOL and visually check; the receptor should be a complete protein, the ligand a single small molecule.

**Selenium / PyMOL import errors at Cell 0** — neither is required for the QSAR pipeline. Comment out those imports if they fail.

---

## How this differs from the published Target2DeNovoDrug

The notebook is an **adapted reimplementation** rather than a direct port. Differences from Madaj et al. (2022):

- The original uses a deep neural network for *de novo* generation; this implementation uses linear regression with polynomial features (more interpretable, much smaller training requirement, but lower predictive power).
- The original retrieves training data from a curated activity DB; this implementation pulls live from PubChem with on-the-fly deduplication.
- The optional docking step uses MGLTools + AutoDock Vina rather than the original's docking backend.

If you have a GPU and want the published deep-learning version, use the original codebase: https://github.com/Bibhash123/Target2DeNovoDrug

---

## References

Madaj, R., Geoffrey, B., Sanker, A. and Valluri, P. P. (2022). *Target2DeNovoDrug: a novel programmatic tool for in silico deep learning based de novo drug design for any target of interest.* Journal of Biomolecular Structure and Dynamics, 40, 14173–14184. https://doi.org/10.1080/07391102.2021.2009566

Cherkasov, A. et al. (2014). *QSAR modeling: where have you been? Where are you going to?* Journal of Medicinal Chemistry, 57, 4977–5010. https://doi.org/10.1021/jm4004285

---

## Citation

If you use this notebook, please cite the original Target2DeNovoDrug paper above and:

> Yang, Y. (2026). *In Silico Drug Discovery in Atherosclerosis.* BSc Dissertation, Faculty of Life Sciences & Medicine.

Repository: https://github.com/yikai-yang/DrugPipe-Target2DeNovoDrug-workflow-by-YIKAI-YANG
