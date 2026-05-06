# DrugPipe + Target2DeNovoDrug — In Silico Drug Discovery Workflow

This repository contains two complementary in silico drug discovery pipelines developed for the dissertation *In Silico Drug Discovery in Atherosclerosis* (Yikai Yang, 2026).

**DrugPipe** starts from a protein structure and generates, docks, and ranks novel and known drug candidates entirely from structure. **Target2DeNovoDrug** starts from public bioassay data, builds a QSAR model to rank existing compounds by predicted activity, then validates top hits with blind docking.

| I have… | Use this pipeline |
|---------|-----------------|
| A protein PDB file, no bioassay data | `drugpipe/` |
| PubChem bioassay data for my target | `target2denovodrug/` |

---

## Repository layout

```
drugpipe/               Structure-based pipeline (PDB → DiffSBDD → Vina → DrugBank)
  config.py             Shared config; set DRUGPIPE_TARGET and DRUGPIPE_PDB
  01_prepare.py … 08_enrich_drugs.py
  batch_run.py          Multi-target batch runner
  README.md

target2denovodrug/      QSAR + docking pipeline (bioassay → model → Vina)
  drug2target_april_2026.ipynb
  dock_top_drugs.py
  build_comparison_table.py
  README.md

shared/                 Gitignored; place drugbank.csv and DiffSBDD here
  README.md             ← instructions for what to download

LICENSE
README.md               ← this file
```

---

## Required external downloads

Both pipelines share external assets that cannot be distributed here:

- **DrugBank CSV** — download from https://go.drugbank.com/releases/latest under your own academic licence; save to `shared/drugbank.csv`.
- **DiffSBDD** — clone to `shared/DiffSBDD/` and place `crossdocked_fullatom_cond.ckpt` in `shared/DiffSBDD/checkpoints/`. See `shared/README.md`.

Set the base-directory variable before running DrugPipe scripts:
```bash
export DRUGPIPE_BASE_DIR=/path/to/your/clone
```

See each subfolder's `README.md` for full setup and usage instructions.

---

## Citation

Yang, Y. (2026). *In Silico Drug Discovery in Atherosclerosis* [Dissertation].

## License

MIT — see `LICENSE`.
