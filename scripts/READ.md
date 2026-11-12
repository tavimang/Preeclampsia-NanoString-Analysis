# Early-Onset vs. Late-Onset Preeclampsia: NanoString Analysis

## Project Overview

This project compares immune-related gene expression in early-onset preeclampsia, late-onset preeclampsia, and healthy placental samples.

## Dataset

- GEO accession: GSE190639
- Data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190639
- Samples: 32 (13 early-onset PE, 6 late-onset PE, 13 healthy controls)
- Genes analyzed: 730

Original NanoString counts and normalized expression data were obtained from GEO.

## Software

- Python 3.11
- R 4.5.3
- limma 3.66.0
- GSEApy 1.3.2

Full dependencies are available in `requirements.txt`, `environment.yml`, and `r_versions.txt`.

## Analysis Workflow

Run the analyses in the following order:

1. `01_data_loading_and_qc.ipynb` — Data loading, normalization assessment, and quality control
2. `02_exploratory_analysis.ipynb` — Exploratory analysis and PCA
3. `scripts/01_limma_analysis.R` — Differential gene expression and gestational-age sensitivity analysis
4. `03_differential_expression.ipynb` — Differential expression results and figures
5. `04_pathway_enrichment.ipynb` — Gene Set Enrichment Analysis

Run the R script from the inner `Preeclampsia_NanoString` folder using `Rscript ../scripts/01_limma_analysis.R`.

The `data` folder contains analysis inputs, `results` contains output tables, and `figures` contains the generated visualizations. Original data can be accessed from the GEO link above; download commands are documented in the data-loading notebook.

## Main Findings

Six genes were significant in the original early-onset versus late-onset comparison. After including gestational age, no genes remained significant at FDR below 0.05.

## AI-Use Statement

ChatGPT was used as a supplemental resource for coding assistance, troubleshooting, and improving clarity. I reviewed and interpreted the findings and made the final research decisions.