# Monocyte miRNA Analysis in Osteoarthritis Pain

## Overview

This repository contains R scripts used to investigate pain-associated miRNA expression patterns across different monocyte subsets.

The analysis compares High Pain and Low Pain groups within three
monocyte subsets:

- Classical monocytes (CLAS)
- Intermediate monocytes (INT)
- Non-classical monocytes (NONCLAS)

The workflow includes principal component analysis (PCA),
differential expression analysis, volcano plot visualisation,
and overlap analysis of candidate miRNAs across monocyte subsets.

## Analysis Workflow

The analysis scripts should be run in the following order:

### 1. PCA

`scripts/01_PCA.R`

Performs principal component analysis separately for CLAS, INT, and NONCLAS monocytes to examine overall miRNA expression patterns between pain groups.

### 2. Differential Expression Analysis

`scripts/02_DE_Volcano.R`

Performs differential expression analysis using the `limma` package.

Comparison:

High Pain vs Low Pain

Interpretation:

- Positive logFC: higher expression in the High Pain group
- Negative logFC: higher expression in the Low Pain group

Candidate miRNAs were identified using the exploratory thresholds:

- P < 0.05
- |logFC| > 0.5

Volcano plots are generated for each monocyte subset.

### 3. Overlap Analysis

`scripts/03_Overlap.R`

Compares candidate miRNAs identified in CLAS, INT, and NONCLAS monocytes.

The analysis identifies:

- Cell-type-specific candidate miRNAs
- miRNAs shared between two or more monocyte subsets
- Direction of expression changes across monocyte subsets

A Venn diagram and shared-miRNA table are generated.

## Repository Structure

```text
OA-Pain/
├── scripts/
│   ├── 01_PCA.R
│   ├── 02_DE_Volcano.R
│   └── 03_Overlap.R
│
├── results/
│   ├── PCA/
│   ├── DE/
│   └── Overlap/
│
├── README.md
└── .gitignore

## Data and Results Availability

The original miRNA expression data, sample metadata, and analysis
outputs are not included in this repository due to data access
restrictions.

This repository contains analysis scripts only. The scripts assume that the required input data are stored locally in a `data/` directory, while generated outputs are saved locally in a `results/` directory.

Both `data/` and `results/` are excluded from version control.

The scripts assume that the required input files are stored locally in a data/ directory.

Expected local structure:
data/
├── miRNA Data Sample Key.csv
└── miRNA_normalized_counts.csv

The data/ directory is excluded from version control using .gitignore.

Software

The analyses were performed in R.

Main R packages used include:

limma
ggplot2
dplyr
ggrepel
patchwork
tidyr
ggforce
ragg

Notes

The differential expression analysis was exploratory. Candidate miRNAs were identified using nominal P values and effect-size thresholds and should not be interpreted as independently validated
biomarkers.