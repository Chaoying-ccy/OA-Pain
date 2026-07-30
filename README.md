# Monocyte miRNA Analysis in Osteoarthritis Pain

## Overview

This repository contains R scripts used to investigate pain-associated miRNA expression patterns across different monocyte subsets.

The analysis compares High Pain and Low Pain groups within three monocyte subsets:

- Classical monocytes (CLAS)
- Intermediate monocytes (INT)
- Non-classical monocytes (NONCLAS)

The analysis workflow includes principal component analysis (PCA), differential expression (DE) analysis, volcano plot visualisation, and overlap analysis of candidate miRNAs across the three monocyte subsets.

## Analysis Workflow

The scripts are designed to be run from the project root directory in the following order:

### 1. Principal Component Analysis (PCA)

Script:

`scripts/01_PCA.R`

This script performs PCA separately for CLAS, INT, and NONCLAS monocytes to examine overall miRNA expression patterns across pain groups.

The PCA plots are combined into a single figure for comparison across the three monocyte subsets.

### 2. Differential Expression Analysis

Script:

`scripts/02_DE_Volcano.R`

This script performs differential expression analysis using the `limma` package.

The comparison is:

**High Pain vs Low Pain**

Therefore:

- Positive logFC indicates relatively higher miRNA expression in the High Pain group.
- Negative logFC indicates relatively higher miRNA expression in the Low Pain group.

Candidate miRNAs are identified using the following exploratory thresholds:

- Nominal P < 0.05
- |logFC| > 0.5

Volcano plots are generated separately for CLAS, INT, and NONCLAS monocytes.

### 3. Overlap Analysis

Script:

`scripts/03_Overlap.R`

This script compares candidate miRNAs identified across the three monocyte subsets.

The analysis identifies:

- Cell-type-specific candidate miRNAs
- Candidate miRNAs shared between two or more monocyte subsets
- Direction of expression changes across monocyte subsets

A Venn diagram and a shared-miRNA table are generated to summarise the overlap.

## Repository Structure

Only analysis scripts and documentation are included in the GitHub repository.

```text
OA-Pain/
├── scripts/
│   ├── 01_PCA.R
│   ├── 02_DE_Volcano.R
│   └── 03_Overlap.R
├── README.md
└── .gitignore
```

## Data and Results Availability

The original miRNA expression data, sample metadata, and analysis outputs are not included in this repository due to data access restrictions.

This repository therefore contains analysis code only.

The scripts assume that the required input data are stored locally in a `data/` directory. Generated analysis outputs are saved locally in a `results/` directory.

Both `data/` and `results/` are excluded from version control using `.gitignore`.

## Expected Local Project Structure

To run the scripts locally, the project is organised as follows:

```text
OA-Pain/
├── data/
│   ├── miRNA Data Sample Key.csv
│   └── miRNA_normalized_counts.csv
├── scripts/
│   ├── 01_PCA.R
│   ├── 02_DE_Volcano.R
│   └── 03_Overlap.R
├── results/
│   ├── PCA/
│   ├── DE/
│   └── Overlap/
├── README.md
└── .gitignore
```

## Running the Analysis

All scripts should be run from the `OA-Pain` project root directory.

The recommended order is:

```text
01_PCA.R
    ↓
02_DE_Volcano.R
    ↓
03_Overlap.R
```

For example, the scripts can be run from the command line using:

```bash
Rscript scripts/01_PCA.R
Rscript scripts/02_DE_Volcano.R
Rscript scripts/03_Overlap.R
```

The scripts use relative paths so that no user-specific absolute file paths are required.

## Software and R Packages

The analyses were performed in R.

The main R packages used include:

- `limma (version 3.66.0)`
- `ggplot2 (version 4.0.3)`
- `dplyr (version 1.2.1)`
- `ggrepel (version 0.9.8)`
- `patchwork (version 1.3.2)`
- `tidyr (version 1.3.2)`
- `ggforce (version 0.5.0)`
- `ragg (version 1.5.2)`
- `R version 4.5.1 (2025-06-13)`
## Statistical Note

The differential expression analysis was exploratory.

Candidate miRNAs were identified using nominal P values and effect-size thresholds (P < 0.05 and |logFC| > 0.5). These candidate signals should therefore not be interpreted as independently validated biomarkers.

## Reproducibility

The analysis scripts are provided to document the computational workflow used in this project. Because the underlying expression data and sample metadata are subject to data access restrictions, the complete analysis cannot be reproduced from this repository alone.
