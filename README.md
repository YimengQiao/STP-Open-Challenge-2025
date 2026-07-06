# STP Open Challenge 2025

> Spatial Transcriptomics → Proteomics (STP): predicting spatial protein abundance from spatial transcriptomics

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](./LICENSE)


## Overview

This is the official Github repo for the [STP Open Challenge](https://www.codabench.org/competitions/10696/). Given spatially resolved transcriptomics (bin-level spatial RNA) together with the matched H&E image, the task is to **predict the spatial abundance of 44 protein markers** at each spatial location (bin).

- **Inputs:** per-bin RNA expression + H&E image patches (optional), with full-resolution pixel coordinates (`pxl_row_in_fullres`, `pxl_col_in_fullres`).
- **Target:** 44 protein markers per bin (e.g. CD31, CD68, GFAP, Ki67, PD-L1, …).
- **Metric:** per-protein Spearman correlation coefficient (SCC); the ranking score is the **mean SCC of the top-10 best-predicted proteins** per submission (see [Evaluation](#evaluation)).
- **Samples:** the training sample and the held-out test sample come from **different tissue sections**, so the challenge probes cross-section generalization.

---

## Repository structure

```
STP-Open-Challenge-2025/
├── data/                         # Data access & preprocessing (no large files in git)
│   ├── README.md                 # Google Drive links + download & layout instructions
│   ├── data_dictionary.md        # Field definitions, protein list, coordinate system
│   └── preprocessing/            # ST-SP alignment
│       ├── registration/         # Image registration between H&E image from Visium HD and DAPI channel image from CODEX
│       └── build_anndata/        # CSV → AnnData (.h5ad); train/valid/test split
│
├── models/
│   ├── README.md                 # Index of team solutions + reproduction status
│   ├── ensemble/                 # Our consensus-based ensemble model over team predictions
│   │   └── README.md
│   └── team_solutions/           # Individual team code (will release with consent)
│
├── evaluation/
│   ├── README.md                 # How to score predictions locally / on Codabench
│   └── scoring.py                # Reference SCC scoring program
│
├── requirment.txt                # Environmental config
├── LICENSE                       # Apache-2.0 (repository code)
└── README.md                     # You are here

```
---
<!-- ## Getting started

```bash
# 1. Clone
git clone <REPO_URL>
cd STP-Open-Challenge-2025

# 2. Create the environment (Python 3.10+ recommended)
#    Core deps: numpy, pandas, scipy, scikit-learn, scanpy, anndata, matplotlib
pip install -r requirements.txt   # or: conda env create -f environment.yml

# 3. Download the data (see data/README.md)
python data/download_data.py      # pulls raw data from Google Drive into data/raw/
```

> Large artifacts (`*.h5ad`, `*.tif`, `*.zip`, model checkpoints, prediction outputs) are **not** tracked in git. They are distributed via Google Drive — see [`data/README.md`](./data/README.md).

---
-->
## Data

Datasets are hosted on Google Drive and downloaded on demand.

- **Access:** [Google Drive](https://drive.google.com/drive/folders/1eq6sbTUaWCCOKcnkei6B65rozx-VX70K?usp=sharing)
- **Layout & download:** [`data/README.md`](./data/README.md)
- **Field / marker reference:** [`data/data_dictionary.md`](./data/data_dictionary.md)

### Preprocessing

The `data/preprocessing/` pipeline turns raw acquisitions into model-ready inputs:

| Step | Location | Purpose |
|---|---|---|
| Manually annotation | / | Manually annotate the corresponding captured regions on the H&E image from Visium HD and DAPI channel from CODEX with QuPath |
| Image registration | `data/preprocessing/registration/` | ST-SP image registration |
| Spatial coordinate mapping | / | Map the spatial coordinates of the Visium HD data onto CODEX data, and aggregate spatial proteomic abundance values with QuPath |
| Build AnnData | / | Assemble `.h5ad` objects with `Scanpy`, generating ST expression matrix & SP expression matrix |

---

## Models

### Team solutions — `models/team_solutions/` (under preparation)

Code from participating teams is released **only with each team's consent**, and each team's subfolder carries its **own license**. Teams that prefer not to vendor their code will be listed with a link to their own repository instead. See [`models/README.md`](./models/README.md) for the per-team index, licenses, and reproduction status.

> If you are a participating team and want to add, update, or remove your solution, please open an issue or PR.

### Ensemble model — `models/ensemble/`

Our in-house ensemble combines the per-protein predictions of the participating teams into a single, stronger predictor. Details and usage in [`models/ensemble/README.md`](./models/ensemble/README.md).

You can run the ensemble model by:
```
python ensemble.py \
    --prediction-dir /path/to/team_predictions \
    --truth-path /path/to/ground_truth.csv \
    --scoring-script /path/to/scoring.py \
    --output-root /path/to/outputs \
    --dataset-name my_dataset \
    --consensus-temperature 1.0 \
    --spatial-k 5
```

---

## Evaluation

Scoring mirrors the Codabench setup and is implemented in [`evaluation/scoring.py`](./evaluation/scoring.py).

**Prediction format** — a single CSV with these ID columns followed by the 44 protein columns:

```
barcode, pxl_row_in_fullres, pxl_col_in_fullres, <protein_1>, ..., <protein_44>
```

**Scoring procedure**

1. Align prediction and ground truth on `barcode` (duplicate barcodes dropped, keeping first).
2. Compute the Spearman correlation (SCC) between predicted and true values **per protein**.
3. Rank proteins by SCC and take the **mean of the top-10** as the final ranking score.

**Run locally**

```bash
python evaluation/scoring.py --name my_submission --gt_path ground_truth_path --pred_path prediction_path
# see evaluation/README.md for input/output conventions
```

More details (Codabench bundle, reference data handling) in [`evaluation/README.md`](./evaluation/README.md).

---

## Reproducing the challenge

1. Download data → [`data/README.md`](./data/README.md)
2. Run preprocessing → `data/preprocessing/`
3. Run a team solution or the ensemble → `models/`
4. Score predictions → `evaluation/`

---

## License

- **Code** in this repository is released under the [Apache License 2.0](./LICENSE).
- **Team solutions** under `models/team_solutions/` are governed by their respective per-folder licenses.
- **Data** is distributed under <DATA_LICENSE> (see [`data/README.md`](./data/README.md)); please cite the original data sources.

---

<!-- ## Citation

If you use this repository, the data, or the challenge in your work, please cite:

```bibtex
@misc{stp_open_challenge_2025,
  title        = {STP Open Challenge 2025: Predicting Spatial Protein Abundance from Spatial Transcriptomics},
  author       = {},
  year         = {2025},
  howpublished = {\url{<REPO_URL>}}
}
```

---
-->
## Acknowledgements

We thank all participating teams for their contributions, and the organizers and data providers of the STP Open Challenge 2025.
