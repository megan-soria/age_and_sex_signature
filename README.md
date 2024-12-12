# Differential Expression Analysis: Age and Sex Signature

This repository contains an RMarkdown workflow for performing differential expression analysis using RNA-seq data. The analysis explores age and sex phenotypes across tissue and cell samples.

## Features

- **Differential Expression Analysis:** Conducted using `DESeq2`.
- **Principal Component Analysis (PCA):** Visualization of sample variance using `PCAtools`.
- **Volcano Plots:** Created with `EnhancedVolcano` to visualize significant genes.
- **Gene Set Enrichment Analysis (GSEA):** Using `mitch` and `msigdbr`.
- **Heatmaps:** High-quality visualizations via `ComplexHeatmap`.

## Repository Contents

- `analysis.Rmd`: The main analysis script.
- `functions.R`: Custom utility functions.
- `rds/`: Folder containing preprocessed RDS files.
- `figures/`: Folder for saving generated plots.
- `data/`: Output CSV files for differential expression results.
- `renv.lock`: Dependency management file.

## Prerequisites

This project requires R (>= 4.0) and the following R packages:

- `prettydoc`
- `readxl`
- `reshape2`
- `tidyverse`
- `DESeq2`
- `edgeR`
- `RUVSeq`
- `PCAtools`
- `EnhancedVolcano`
- `ComplexHeatmap`
- `RColorBrewer`
- `circlize`
- `kableExtra`
- `knitr`
- `mitch`
- `msigdbr`

### Installing Dependencies

You can install all required dependencies using the command below:

```r
install.packages(c("prettydoc", "readxl", "reshape2", "tidyverse", "DESeq2",
                   "edgeR", "RUVSeq", "PCAtools", "EnhancedVolcano",
                   "ComplexHeatmap", "RColorBrewer", "circlize", "kableExtra",
                   "knitr", "mitch", "msigdbr"))
```

Alternatively, restore the environment using `renv`:

```r
renv::restore()
```

## Usage

1. Clone this repository:
   ```bash
   git clone https://github.com/username/repository-name.git
   cd repository-name
   ```
2. Restore the R environment:
   ```r
   renv::restore()
   ```
3. Open `analysis.Rmd` in RStudio and run the workflow.

## Outputs

- **PCA Plots:** Visualize sample variance.
- **Volcano Plots:** Highlight significant genes.
- **Differential Expression Results:** Saved as CSV files in the `data/` folder.
- **GSEA Results:** Pathway enrichment results visualized and exported.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Authors

Developed by MSoria

