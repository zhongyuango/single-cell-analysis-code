\# Single-cell transcriptomic analysis code



This repository contains the R scripts used for the single-cell

transcriptomic analyses and figure generation in this study.



\## Repository structure



\- `scripts/Figure1\_2.R`: Data preprocessing, quality control,

&#x20; integration, clustering and cell-type annotation.

\- `scripts/Figure3.R`: Slingshot, Monocle2 and tradeSeq trajectory analyses.

\- `scripts/Figure4.R`: CellChat-based cell-cell communication analysis.

\- `scripts/Figure5.R`: copyKAT, GSVA, hdWGCNA and survival analyses.

\- `config/sample\_info\_template.csv`: Template for sample information.

\- `data/README.md`: Description of the required input data.



\## Requirements



The analyses were performed in R. Major packages include:



\- Seurat

\- Harmony

\- DoubletFinder

\- SingleR

\- Slingshot

\- Monocle

\- tradeSeq

\- CellChat

\- copyKAT

\- GSVA

\- hdWGCNA



\## Data availability



The raw and processed sequencing data are not stored in this repository.

Users should download the corresponding datasets from the data repository

described in the manuscript.



Large intermediate objects, including RDS files, are not included.



\## Usage



Place the required input files in the corresponding `data` directory

and run the scripts in numerical order:



```r

source("scripts/Figure1\_2.R")

source("scripts/Figure3.R")

source("scripts/Figure4.R")

source("scripts/Figure5.R")

