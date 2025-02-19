# Analysis of Age and Sex Signature in Human Muscle Tissue and Cell

This repository contains an RMarkdown workflow for analysis using RNA-seq data in humans. The analysis explores age and sex phenotypes across tissue and cell samples. Fastq files from this analysis is available at NCBI GEO (GSE287342). Data for time-course and  differentiation score analysis is available at NCBI GEO (GSE168897).

## Features

- **Differential Expression Analysis:** Conducted using `DESeq2`.
- **Principal Component Analysis (PCA):** Visualization of sample variance using `PCAtools`.
- **Volcano Plots:** Created with `EnhancedVolcano` to visualize significant genes.
- **Gene Set Enrichment Analysis (GSEA):** Using `mitch` and `msigdbr`.

## Repository Contents

- `1_DE_analysis.Rmd`: Differential expression analysis script.
  - Figures 1 to 3
  - Supplementary Figure 1
- `2_timeCourse.Rmd`: Time-course analysis to account for cell differentiation.
  - Figure 4
  - Supplementary Figures 2 and 3
- `3_signatures.Rmd`: Multi-comparison analysis to get age and sex signature.
  - Figures 5 and 6
  - Supplementary Figure 5 and 6
- `4_differentiationScore.Rmd`: Differentiation score to assign to samples.
  - Supplementary Figure 4
  - Table 2
- `functions.R`: Custom utility functions.

## Prerequisites

This project requires R (>= 4.0) and the following R packages:

**Core Packages:**

*   `prettydoc`
*   `tidyverse` (includes `dplyr` and `ggplot2`)
*   `DESeq2`
*   `RUVSeq`
*   `PCAtools`
*   `EnhancedVolcano`
*   `ComplexHeatmap`
*   `knitr`
*   `mitch`
*   `msigdbr`

**Visualization & Plotting:**

*   `ggrepel`
*   `gplots`
*   `ggpubr`
*   `enrichplot`
*   `cowplot`
*   `RColorBrewer`
*   `circlize`
*   `viridis`

**Bioinformatics & Annotation:**

*   `RNAseqQC`
*   `clusterProfiler`
*   `DOSE`
*   `org.Hs.eg.db`
*   `DEGreport`


## Authors

Developed by MSoria

