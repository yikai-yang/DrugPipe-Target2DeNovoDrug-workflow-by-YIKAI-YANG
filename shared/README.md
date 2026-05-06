# shared/

This directory is excluded from version control (`shared/*` in `.gitignore`) except for this README.
You must download or clone the following assets yourself before running the pipeline.

## Required downloads

### 1. DrugBank CSV

Download the full DrugBank database CSV from:
https://go.drugbank.com/releases/latest

You will need a free academic account. Save the file to:

```
shared/drugbank.csv
```

### 2. DiffSBDD model and checkpoint

Clone the DiffSBDD repository into this folder:

```bash
git clone https://github.com/arneschneuing/DiffSBDD shared/DiffSBDD
```

Then download the pre-trained checkpoint `crossdocked_fullatom_cond.ckpt` from the DiffSBDD release page and place it at:

```
shared/DiffSBDD/checkpoints/crossdocked_fullatom_cond.ckpt
```

## Environment variable

Set `DRUGPIPE_BASE_DIR` to the absolute path of your clone before running any pipeline script:

```bash
export DRUGPIPE_BASE_DIR=/path/to/your/clone
```

The scripts will fall back to the repository root if this variable is not set.
