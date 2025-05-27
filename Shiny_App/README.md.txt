# Shiny Interface for R-ODAF DEG Pipeline

This Shiny app provides a graphical interface for the R-ODAF DESeq2 pipeline.  
It enables fixed-parameter differential gene expression analysis with built-in visualization and reporting.

## Features

- Input support for RNA-seq count matrices and sample metadata
- Fixed filters for regulatory consistency:
  - Minimum coverage: 5 million reads/sample
  - CPM threshold: 1
  - FDR cutoff: 0.01
  - PCA outlier variance threshold: 20%
- Interactive DEG exploration:
  - PCA before/after outlier removal
  - Final DEG tables with export
  - Heatmaps, top gene, and average expression plots
- Exportable results and summary reports

## Features for TempO-Seq will be added soon

## How to Run

1. Open `app.R` in RStudio
2. Click **“Run App”**
3. Upload:
   - A count matrix (genes × samples)
   - A metadata file (samples × attributes)
4. Set design column and control/case labels
5. Click **“Run analysis”**

## Notes

- All thresholds are hardcoded to ensure reproducibility and reduce user-induced variability.
- Designed for RNA-seq input. (Tempo-Seq support under consideration)
- No changes were made to the core R-ODAF pipeline logic.

## Citation

Please cite the original R-ODAF publication:

Verheijen, M. C., Meier, M. J., Asensio, J. O., Gant, T. W., Tong, W., Yauk, C. L., & Caiment, F. (2022). R-ODAF: Omics data analysis framework for regulatory application. Regulatory Toxicology and Pharmacology, 131, 105143.