# Single-cell transcriptomic analysis code

This repository contains the R scripts used for the single-cell transcriptomic analyses and figure generation in this study.

## Repository structure

- `scripts/01_Figure_1-2.R`: Data preprocessing, quality control, integration, clustering and cell-type annotation.
- `scripts/02_Figure_3.R`: Slingshot, Monocle2 and tradeSeq trajectory analyses.
- `scripts/03_Figure_4.R`: CellChat-based cell-cell communication analysis.
- `scripts/04_Figure_5.R`: copyKAT, GSVA, hdWGCNA and survival analyses.
- `config/README.md`: Description of configuration requirements.
- `data/README.md`: Description of the required input data.

## Requirements

The analyses were performed in R. Major packages include:

- Seurat
- Harmony
- DoubletFinder
- SingleR
- Slingshot
- Monocle
- tradeSeq
- CellChat
- copyKAT
- GSVA
- hdWGCNA

## Data availability

The raw and processed sequencing data are not stored in this repository.

Users should download the corresponding datasets from the data repository described in the manuscript.

Large intermediate objects, including RDS files, are not included.

## Usage

Place the required input files in the corresponding `data` directory and run the scripts in numerical order:

```r
source("scripts/01_Figure_1-2.R")
source("scripts/02_Figure_3.R")
source("scripts/03_Figure_4.R")
source("scripts/04_Figure_5.R")
```

Input file names, metadata columns and analysis parameters should be adjusted according to the user's dataset.
