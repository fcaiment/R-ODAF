# Transcriptomic data analysis with R-ODAF

**Author: Saad Lodhi**

## 1. Introduction

The R-ODAF app is a Shiny-based tool for standardised RNA-seq differential expression analysis. It implements the R-ODAF framework in a point-and-click interface so users can go from counts + metadata to quality-controlled differentially expressed genes (DEGs), plots and exportable tables.

You use this tool when you have a two-group RNA-seq experiment (for example, control vs treated) and need a reproducible, harmonised analysis suitable for risk assessment, pathway analysis or adverse outcome pathway (AOP) workflows within VHP4Safety or in general. It addresses the problem of inconsistent in-house RNA-seq pipelines by enforcing fixed thresholds (minimum coverage, CPM, FDR, outlier rules) and a transparent, DESeq2-based workflow.

Key assumptions and limitations:

- Requires a counts matrix (genes × samples) and a metadata table (samples × attributes) with matching sample IDs.

- Designed for two-group comparisons (control vs case).

- Automatic gene symbol mapping is implemented for human, mouse and rat Ensembl gene IDs; other species are analysed at the ID level only.

- QC and filter settings (coverage, CPM, FDR, outlier variance) are fixed and not user-configurable.

## 2. Accessing the Tool

The app is available at:

https://r-odaf.nl/

It runs in a modern browser (Chrome, Firefox, Edge). No installation is needed.

You need:

- a counts file (CSV or tab-delimited): rows = gene IDs, columns = sample IDs;

- a metadata file (CSV or tab-delimited): rows = sample IDs, columns = variables (for example, Treatment, Dose, Time);

- a metadata column that defines the comparison (for example, treatment with levels “control” and “case”).

To load data:

- Open https://r-odaf.nl/
 and go to Analysis.

- Upload the counts file under Counts (genes × samples).

- Upload the metadata file under Metadata (samples × attributes).

- Open Data preview and check that the sample IDs match between counts and metadata.

## 3. Tool Functionalities
### 3.1 Data specification (upload + design)

Purpose:

Connect your counts and metadata, and define which samples are control and which are case.

How it works

- The app reads both files, keeps only overlapping sample IDs, and builds a DESeq2 model using the selected design column. The control level is used as the reference, and the case level is compared against it.

How to use

- Upload both files in Analysis, then set Design column, Control label and Case label in the sidebar (exactly as they appear in your metadata). - Use Data preview to sanity-check the structure.
- Example: Design column = treatment, Control = control, Case = doxorubicin.

### 3.2 Running the R-ODAF pipeline

Purpose
Execute the full, standardised DE pipeline in one step.

How it works
When you click Run analysis, the app:

- filters out low-coverage samples,

- runs an initial DESeq2 model,

- performs PCA-based outlier detection and removes problematic samples,

- applies a 75% CPM-based gene filter,

- re-runs DESeq2 with fixed FDR = 0.01 and R-ODAF spike filters,

- generates DEGs, normalised counts and plots.

How to use:

After upload and design selection, click Run analysis. When the progress bar reaches 100% and shows “Done.”, all result tabs are filled.

### 3.3 Quality control and outliers (PCA & heatmaps)

Purpose:
Check whether your samples behave as expected and whether any were removed as outliers.

How it works:

The PCA & outliers tab shows PCA plots before and after outlier removal plus a simple outlier log. The DE results tab contains:

- a correlation heatmap for all samples (global similarity check),

- a DEG heatmap showing z-scaled expression of significant DEGs.

- All plots can be downloaded as PNG, TIFF, SVG or JPEG.

### 3.4 Differential expression and gene-level views

Purpose:

Inspect and export DE results, and visualise top genes per group.

DE results shows:

- the total number of genes tested and final DEGs (FDR < 0.01),

- an interactive table of DEGs with log₂ fold change, p-values, FDR and (for human) gene symbols,

- download buttons for DEGs and normalised counts (CSV or TXT).

- Top gene expression shows boxplots of normalised counts per group for windows of top-ranked DEGs.

- Avg expression shows barplots of mean normalised counts per group for the same genes.

- In both, gene labels can be set to Ensembl IDs or gene symbols (e.g., for human data), and plots can be downloaded.


- In DE results, review the summary and table, then export DEGs and normalised counts if needed. 
- In Top gene expression and Avg expression, set the number of DEGs and Gene range slider to browse through top genes and check that differences are consistent and not driven by single samples.

### 3.5 Summary, documentation, citation and support

Purpose
Provide a concise record of the analysis and all supporting information.

How it works

Summary: key settings (design, thresholds), number of tested genes and DEGs, recommended citations, plus download buttons for a text summary and sessionInfo() for reproducibility.

Tutorials: renders markdown tutorials (for example, the human doxorubicin organoid example). 
A step-by-step guided example is available directly in the app under the Tutorials tab.

Documentation: short in-app description of input format and fixed filters.

Citation: recommended references for R-ODAF-Shiny and the original R-ODAF paper.

License: MIT licence text.

Contact / Help: GitHub issues link and email address.

How to use
After running an analysis, use Summary to download the report and sessionInfo for archiving or supplementary material. Use Citation when writing manuscripts and Contact / Help to report issues or ask questions.

## 4. Interpreting the Output

The app produces:

- PCA plots (before/after outlier removal) and an outlier log,

- a DE summary + DEG table (with log₂ fold changes and FDR),

- sample correlation and DEG heatmaps,

- per-gene box and bar plots for top DEGs,

- exportable DEGs, normalised counts, and a text summary + sessionInfo.

In a VHP4Safety context:

- Use PCA & correlation to confirm that the data are coherent and that obvious technical failures have been removed.

- Use the DEG table to identify statistically robust treatment-responsive genes (FDR < 0.01) for input into pathway, AOP or modelling tools.


Caveats:

- The current app is designed for two-group designs; complex designs (time, dose, interactions) should be analysed in R outside the app.

- Fixed thresholds may not be ideal for every dataset; always interpret results alongside experimental QC and biological knowledge.

## 5. Summary / Conclusion

The R-ODAF app delivers a standardised, transparent and easy-to-use pipeline for RNA-seq differential expression:

same QC and filtering rules across studies,

a clear DESeq2-based workflow with fixed FDR and outlier handling,

ready-to-use plots, tables and reports for downstream work.

Within VHP4Safety, it should be used whenever a two-group RNA-seq comparison needs to be analysed in a consistent, reproducible way and linked to further tools (pathway enrichment, qAOPs, cross-model comparisons). Users can start at https://r-odaf.nl/
, follow the built-in tutorial, and use the Summary, Citation and Contact / Help tabs to document and support their analyses.

In summary, R-ODAF provides the standardised transcriptomic step in the VHP4Safety process flow, delivering inputs that connect directly to Toxicodynamics/AOPs.