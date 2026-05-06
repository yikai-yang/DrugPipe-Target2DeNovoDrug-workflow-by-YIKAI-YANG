# DrugPipe — Structure-Based Virtual Screening Pipeline

Locally adapted reimplementation of **DrugPipe** (Pham et al. 2025, *Biology Methods and Protocols*) for structure-based drug repurposing against poorly-characterised protein targets. Built and used during the BSc dissertation *"In Silico Drug Discovery in Atherosclerosis"* (Yikai Yang, 2026).

This pipeline is the route taken when a target has **no usable bioassay data** for QSAR modelling — i.e. when you need to start from the protein structure rather than from known active compounds. For the QSAR/ligand-based route, see `Target2DeNovoDrug_README.md`.

---

## What the pipeline does

For each protein target, the pipeline:

1. Detects druggable pockets on the protein surface
2. Generates pocket-specific *de novo* ligands using a diffusion model (DiffSBDD)
3. Docks those generated ligands with AutoDock Vina to establish a binding-energy baseline per pocket
4. Searches the entire DrugBank library (9,716 approved drugs) for compounds chemically similar to the generated ligands using Morgan fingerprints
5. Re-docks the top-100 most similar DrugBank drugs into each pocket to obtain real Vina scores
6. Aggregates results, applies drug-likeness filters, and enriches the top hits with external annotations (PubChem, ChEMBL, indications, clinical phase)

Final output: a ranked, annotated CSV of repurposing candidates per target, plus a cross-target comparison.

---

## Pipeline overview

```
                     ┌──────────────────────────┐
   Input PDB  ───▶   │  01_prepare.py           │   fpocket + receptor PDBQT
                     └──────────────┬───────────┘
                                    │
                                    ▼
                     ┌──────────────────────────┐
                     │  02_generate.py          │   DiffSBDD → ~50 ligands per pocket
                     └──────────────┬───────────┘
                                    │
                                    ▼
                     ┌──────────────────────────┐
                     │  03_dock_generated.py    │   AutoDock Vina baseline scores
                     └──────────────┬───────────┘
                                    │
                                    ▼
                     ┌──────────────────────────┐
                     │  04_morgan_rank.py       │   DrugBank similarity ranking
                     └──────────────┬───────────┘
                                    │
                                    ▼
                     ┌──────────────────────────┐
                     │  05_dock_drugbank.py     │   Re-dock top-100 DrugBank drugs
                     └──────────────┬───────────┘
                                    │
                                    ▼
                     ┌──────────────────────────┐
                     │  06_aggregate.py         │   Per-target summary
                     └──────────────┬───────────┘
                                    │
                                    ▼
                     ┌──────────────────────────┐
                     │  07_analyze.py           │   Drug-likeness + chemotype clusters
                     └──────────────┬───────────┘
                                    │
                                    ▼
                     ┌──────────────────────────┐
                     │  08_enrich_drugs.py      │   PubChem / ChEMBL annotations
                     └──────────────────────────┘
```

`batch_run.py` orchestrates Steps 3–5 across all targets in one go (recommended for a multi-target screen).

---

## Folder layout

The pipeline expects a fixed folder structure. Set `BASE_DIR` in `config.py` to wherever you put it.

```
Whole_DrugPipe/
├── PDB/                                 ← put your input .pdb files here
│   ├── IGFBP7_energy_minimized.pdb
│   ├── BGN_energy_minimized.pdb
│   └── ...
│
├── shared/                              ← shared assets, set up once per machine
│   ├── drugbank.csv                     ← 9,716 approved drugs (name, id, SMILES)
│   └── DiffSBDD/                        ← cloned DiffSBDD repo
│       └── checkpoints/
│           └── crossdocked_fullatom_cond.ckpt
│
├── scripts/                             ← all pipeline scripts live here
│   ├── config.py
│   ├── 01_prepare.py
│   ├── 02_generate.py
│   ├── 03_dock_generated.py
│   ├── 04_morgan_rank.py
│   ├── 05_dock_drugbank.py
│   ├── 06_aggregate.py
│   ├── 07_analyze.py
│   ├── 08_enrich_drugs.py
│   └── batch_run.py
│
├── runs/                                ← auto-created per-target output
│   └── <TARGET_NAME>/
│       ├── input/<TARGET>.pdb           ← canonical copy of input
│       ├── receptor/<TARGET>.pdbqt      ← prepared receptor
│       ├── fpocket/<TARGET>_out/        ← raw fpocket output
│       ├── <TARGET>_pockets_summary.csv ← ranked pockets (Step 1)
│       ├── pocket1/                     ← one folder per selected pocket
│       │   ├── <TARGET>_pocket1_marker.sdf
│       │   ├── <TARGET>_pocket1_ligands.sdf
│       │   ├── <TARGET>_pocket1_ligands_clean.sdf
│       │   ├── <TARGET>_pocket1_generated_docking.csv
│       │   ├── <TARGET>_pocket1_morgan.csv
│       │   ├── <TARGET>_pocket1_final.csv
│       │   ├── generated_ligands_pdbqt/
│       │   ├── generated_poses/
│       │   └── drugbank_poses/
│       ├── pocket2/                     ← (same structure)
│       └── drugbank_ligands_pdbqt/      ← shared PDBQT cache (reused across pockets)
│
└── analysis/                            ← cross-target analysis output
    ├── <TARGET>/
    │   ├── <TARGET>_all_ranked.csv
    │   ├── <TARGET>_clean_shortlist.csv
    │   ├── <TARGET>_top30_with_chemotypes.csv
    │   ├── <TARGET>_pocket_stats.csv
    │   └── <TARGET>_enriched_top20.csv
    ├── cross_target_summary.csv
    └── all_targets_enriched.csv
```

---

## Software requirements

The pipeline was developed and tested on **macOS** with a CPU-only MacBook Air (10 cores, no GPU). It should run on Linux with minor path adjustments.

### Conda environment

Create a fresh conda env (recommended name: `sbdd-env`):

```bash
conda create -n sbdd-env python=3.9.6
conda activate sbdd-env
```

### Python packages

| Package | Version | Purpose |
|---|---|---|
| `rdkit` | 2026.03.1 | Fingerprints, SMILES, descriptors |
| `pandas` | 3.0.2 | CSV processing |
| `numpy` | 2.4.4 | Numerical processing |
| `biopython` | 1.79 | PDB handling |
| `pytorch-lightning` | 2.4.0 | DiffSBDD runtime |
| `pytorch-scatter` | 2.1.2 | DiffSBDD dependency |
| `meeko` | 0.7.1 | Docking file prep |
| `requests` | 2.33.1 | External database queries |

### Command-line tools (must be in `$PATH`)

| Tool | Version | Install |
|---|---|---|
| AutoDock Vina | 1.2.5 | `brew install autodock-vina` (mac) or download from CCSB |
| Open Babel | 3.1.1 | `brew install open-babel` |
| fpocket | 4.2.3 | `brew install fpocket` |

Verify installation:

```bash
which vina obabel fpocket
vina --version
```

### DiffSBDD setup

```bash
cd shared/
git clone https://github.com/arneschneuing/DiffSBDD.git
cd DiffSBDD
# place crossdocked_fullatom_cond.ckpt in DiffSBDD/checkpoints/
```

The checkpoint can be downloaded from the [DiffSBDD release page](https://github.com/arneschneuing/DiffSBDD/releases). Expected file: `crossdocked_fullatom_cond.ckpt` (~17.9 MB).

### DrugBank library

`shared/drugbank.csv` should be a CSV with at minimum these columns:

```
Drug name, Drug id, smiles
```

The original DrugPipe paper uses the 9,716-compound library distributed with their repository.

---

## Configuration

Everything pipeline-wide is in `scripts/config.py`. Open it once, edit the **BASE_DIR** and per-machine paths, then leave it alone.

The two values you change *per target*:

```python
TARGET_NAME = "IGFBP7"
PDB_FILE    = "/Users/yikai/Desktop/Whole_DrugPipe/PDB/IGFBP7_energy_minimized.pdb"
```

You can also set them via environment variables (this is what `batch_run.py` uses):

```bash
export DRUGPIPE_TARGET="IGFBP7"
export DRUGPIPE_PDB="/path/to/IGFBP7.pdb"
```

### Key parameters and defaults

| Parameter | Default | What it controls |
|---|---|---|
| `N_POCKETS` | 2 | Top-N druggable pockets to process per target (use 3 for the primary target) |
| `MIN_DRUGGABILITY` | 0.0 | Skip pockets below this fpocket score; 0 = take top-N regardless |
| `N_LIGANDS` | 30 | DiffSBDD samples per pocket (paper uses 50–100; expect 30–70% to survive sanitisation) |
| `DIFFSBDD_TIMESTEPS` | 500 | More = slower but better quality; 500 is the DiffSBDD default |
| `BOX_SIZE` | 22.0 Å | Cubic Vina docking box side length |
| `EXHAUSTIVENESS` | 8 | Vina default; lower = faster, higher = more thorough |
| `N_POSES` | 5 | Poses Vina reports per ligand |
| `TOP_N_DOCK` | 100 | Top Morgan-FP hits per pocket sent to docking |
| `MORGAN_RADIUS` | 2 | ECFP4 — standard for drug similarity |
| `MORGAN_NBITS` | 2048 | Fingerprint bit length |
| `N_WORKERS` | 8 | Parallel docking jobs |
| `DOCKING_TIMEOUT` | 600 s | Kills runaway dockings |

To sanity-check the configuration:

```bash
cd scripts/
python config.py
```

This prints all current settings and warns about any missing files.

---

## Running the pipeline — single target

For one target, run the steps in order from inside `scripts/`:

```bash
cd scripts/

# Step 1: pocket detection + receptor preparation (~1 min)
python 01_prepare.py

# Step 2: DiffSBDD ligand generation (~10–15 min per pocket on CPU)
python 02_generate.py

# Step 3: dock generated ligands (~5–10 min for 50 ligands × 3 pockets, 8 workers)
python 03_dock_generated.py

# Step 4: Morgan FP ranking against DrugBank (~2–3 min per pocket)
python 04_morgan_rank.py

# Step 5: dock top-100 DrugBank drugs per pocket (~20–30 min for 3 pockets, 8 workers)
python 05_dock_drugbank.py

# Step 6: aggregate per-target results
python 06_aggregate.py

# Step 7: drug-likeness analysis and chemotype clustering
python 07_analyze.py

# Step 8: external enrichment (PubChem + ChEMBL)
python 08_enrich_drugs.py
```

Every step is **idempotent** — re-running a step with completed outputs already present will skip the work. To force re-run, delete the relevant output files first.

After Step 1, **review `runs/<TARGET>/<TARGET>_pockets_summary.csv`** before continuing. If no pockets reach a useful druggability score (>0.3), the target may not be tractable for structure-based screening.

---

## Running the pipeline — many targets at once

For batch screens across many targets, prepare each target's PDB and run Step 1 + Step 2 manually for each (these are too sensitive to fold into the batch script). Then:

```bash
# Edit the TARGETS list at the top of batch_run.py, then:
python batch_run.py
```

`batch_run.py` runs Steps 3–5 across the listed targets, with:

- 2-hour hard timeout per step (kills stuck runs, continues to the next target)
- Per-target exception handling (one failure doesn't abort the others)
- Full stdout/stderr captured per target in `runs/<TARGET>/batch_<step>.log`
- Truncated output detection (deletes corrupt CSVs <200 bytes so re-runs are clean)
- Final `batch_summary.csv` with the status of every target × step

Watch progress live:

```bash
tail -f batch_run.log
```

Once `batch_run.py` finishes, run Steps 6–8 once across everything:

```bash
python 07_analyze.py                 # all targets in runs/
python 08_enrich_drugs.py            # all targets, top-20 each
# or for specific targets:
python 07_analyze.py IGFBP7 BGN
TOP_N=50 python 08_enrich_drugs.py IGFBP7
```

---

## What each step produces

### Step 1 — `01_prepare.py`

Cleans the PDB (strips waters, keeps common metal cofactors), prepares the receptor PDBQT via Open Babel, runs fpocket, and ranks pockets by druggability.

**Outputs**:
- `runs/<TARGET>/input/<TARGET>.pdb` — canonical copy of input
- `runs/<TARGET>/receptor/<TARGET>.pdbqt` — prepared receptor
- `runs/<TARGET>/fpocket/<TARGET>_out/` — raw fpocket output (pocket atm.pdb files etc.)
- `runs/<TARGET>/<TARGET>_pockets_summary.csv` — ranked pocket table with columns: `druggability_rank, pocket_num, selected, druggability, score, volume, num_alpha_spheres, hydrophobicity, polarity, center_x, center_y, center_z`

The `selected` column flags which pockets get processed downstream.

### Step 2 — `02_generate.py`

Writes a fake-ligand SDF marker at each selected pocket's centre, runs DiffSBDD to generate ~30 candidate ligands, and sanitises them through RDKit.

**Outputs (per selected pocket N)**:
- `<TARGET>_pocket<N>_marker.sdf` — pocket centre marker (single carbon atom)
- `<TARGET>_pocket<N>_ligands.sdf` — raw DiffSBDD output
- `<TARGET>_pocket<N>_ligands_clean.sdf` — RDKit-sanitised valid molecules

### Step 3 — `03_dock_generated.py`

Converts each generated ligand to PDBQT, docks against the receptor with a pocket-centred 22 Å box, extracts the best Vina score per ligand.

**Outputs (per pocket)**:
- `<TARGET>_pocket<N>_generated_docking.csv` — columns: `pocket_num, ligand_id, smiles, ba, x, y, z, status`
- `pocket<N>/generated_ligands_pdbqt/` — PDBQT files
- `pocket<N>/generated_poses/` — Vina pose `.pdbqt` outputs

### Step 4 — `04_morgan_rank.py`

For every DrugBank drug, finds the most-similar generated ligand by Tanimoto on Morgan fingerprints, inherits that ligand's binding energy, and ranks by a 50/50 score combining similarity and binding energy.

**Outputs (per pocket)**:
- `<TARGET>_pocket<N>_morgan.csv` — columns: `rank_in_pocket, pocket_num, Drug name, Drug id, smiles, tanimoto_score, inherited_ba, combined_score, matched_ligand_idx, matched_ligand_smi`

### Step 5 — `05_dock_drugbank.py`

Takes the top-100 DrugBank candidates per pocket, prepares 3D PDBQTs (cached and reused across pockets), docks them with Vina.

**Outputs (per pocket)**:
- `<TARGET>_pocket<N>_final.csv` — columns: `final_rank, pocket_num, Drug name, Drug id, smiles, tanimoto_score, inherited_ba, vina_score, combined_score, dock_status`
- `pocket<N>/drugbank_poses/<DBID>_out.pdbqt` — Vina poses for each docked drug

### Step 6 — `06_aggregate.py`

Combines per-pocket results into per-target tables.

### Step 7 — `07_analyze.py`

Annotates drug-likeness (MW, logP, HBD/HBA, rotatable bonds, Lipinski, common PAINS scaffolds), clusters the top-30 by chemotype, and builds a "clean shortlist" of drug-like non-PAINS hits.

**Outputs (per target)**:
- `analysis/<TARGET>/<TARGET>_all_ranked.csv` — best pocket per drug, ranked
- `analysis/<TARGET>/<TARGET>_clean_shortlist.csv` — drug-like, non-PAINS, vina ≤ −6.0
- `analysis/<TARGET>/<TARGET>_top30_with_chemotypes.csv` — top-30 with cluster IDs
- `analysis/<TARGET>/<TARGET>_pocket_stats.csv` — per-pocket statistics

**Cross-target output**:
- `analysis/cross_target_summary.csv` — best hit per target across the whole screen

### Step 8 — `08_enrich_drugs.py`

For each of the top-N candidates per target (default N=20), queries:

- **PubChem** — CID, synonyms, PubChem URL (via InChIKey lookup)
- **UniChem** — DrugBank ID → ChEMBL ID crosswalk
- **ChEMBL** — known targets, approved indications, max clinical phase, indication class

**Outputs**:
- `analysis/<TARGET>/<TARGET>_enriched_top20.csv`
- `analysis/all_targets_enriched.csv` — master cross-target table

The end of Step 8 prints any compounds that are **already approved** (max clinical phase ≥ 4) and appear among your top hits — these are your most actionable repurposing candidates.

---

## Where to find your final results

For a quick look at the best candidate per target:

```
analysis/cross_target_summary.csv
```

For the fully annotated top hits across all targets, with PubChem/ChEMBL info:

```
analysis/all_targets_enriched.csv
```

For a single target's full ranked list:

```
analysis/<TARGET>/<TARGET>_clean_shortlist.csv
```

For visualising poses, open the receptor and a pose file together in PyMOL:

```bash
pymol runs/<TARGET>/receptor/<TARGET>.pdbqt \
      runs/<TARGET>/pocket2/drugbank_poses/<DRUGBANK_ID>_out.pdbqt
```

---

## Common pitfalls

**`receptor PDBQT is empty or tiny`** during Step 1 — Open Babel failed. Check `obabel --version`; on some macOS installs the partial-charge flag fails silently. Try running `obabel <file>.pdb -O <file>.pdbqt -xr -p 7.4 --partialcharge gasteiger` manually to see the error.

**fpocket finds zero pockets** — usually the structure has chain breaks or non-standard residues confusing fpocket. Open the cleaned PDB in PyMOL and check.

**DiffSBDD hangs with no output** — the unbuffered subprocess monitor in `02_generate.py` will print `...still running, no new output in 30s...`. CPU-only DiffSBDD genuinely takes ~10s/ligand. If it's been silent for >5 min, kill it and check `DIFFSBDD_TIMESTEPS` (try lowering to 200 to test).

**Vina scores all NaN in Step 3** — your Vina box is probably wrong. Check pocket centres in `<TARGET>_pockets_summary.csv` are inside the protein bounding box.

**Step 8 enrichment is slow / fails on some drugs** — PubChem/ChEMBL rate-limit at ~5 req/s. The script uses retries with backoff. If many fail, lower `MAX_WORKERS` from 6 to 3.

**Truncated output CSVs after a crash** — `batch_run.py` deletes anything <200 bytes automatically. If you ran a step manually and it crashed, check the file size with `ls -l` and delete by hand if needed; the next run will skip files that exist and have data.

---

## Reference

Pham, P., Cheng, J., Wu, S., Zhang, Y. and Chen, J. (eds.) 2025. *DrugPipe: generative artificial intelligence-assisted virtual screening pipeline for generalizable and efficient drug repurposing.* Biology Methods and Protocols, 10, bpaf038. https://doi.org/10.1093/biomethods/bpaf038

This implementation differs from the original DrugPipe in that it uses **Morgan fingerprints** for the similarity step rather than the original graph neural network. The published `encoder_best.pt` GNN weights were not available, and pre-training a replacement requires GPU time beyond a typical undergraduate setup.

---

## Citation

If you use this pipeline, please cite the original DrugPipe paper above and:

> Yang, Y. (2026). *In Silico Drug Discovery in Atherosclerosis.* BSc Dissertation, Faculty of Life Sciences & Medicine.

Repository: https://github.com/yikai-yang/DrugPipe-Target2DeNovoDrug-workflow-by-YIKAI-YANG
