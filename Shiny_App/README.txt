R-ODAF-Shiny: Graphical Interface for the R-ODAF DEG Pipeline

This Shiny app provides a graphical interface for the R-ODAF DESeq2 pipeline.
It enables fixed-parameter differential gene expression (DEG) analysis with built-in visualization and reporting.
All analyses follow R-ODAF’s regulatory criteria for reproducibility and traceability.

Live Deployment

The containerized app is available at:
https://r-odaf.nl

Uploads are processed in-session and discarded automatically after session end.

Key Features

RNA-seq counts + metadata input

Filters for regulatory consistency:

Minimum coverage: 5M reads/sample

CPM threshold: 1

FDR cutoff: 0.01

PCA outlier variance threshold: 20%

Interactive exploration:

PCA (before/after outlier removal)

Final DEG tables (exportable)

Heatmaps, top-gene and average-expression plots

Exportable results and summary reports

TempO-Seq module is planned.

Run Locally (short version)
install.packages("renv", repos = "https://cloud.r-project.org")
renv::restore()        # restore exact package versions
shiny::runApp(".")     # launch the app


Inputs expected

counts: genes × samples (CSV/TSV)

metadata: samples × attributes (must include a condition column)

Reproducibility

This repository includes:

renv.lock – snapshot of all package versions

renv/activate.R – automatic environment activation

docs/sessionInfo.txt – R/session details

Fixed analysis thresholds to ensure identical results across runs

See also:

docs/PRIVACY.md – data handling

docs/REPRODUCIBILITY.md – fixed parameters and validation notes

Citation

Please cite the R-ODAF framework and this Shiny implementation:


R-ODAF-Shiny (this work)
Saad Lodhi, Marcha Verheijen, Theo M. de Kok, Florian Caiment, Danyel Jennen (2025).
R-ODAF-Shiny: A Graphical Interface for Reproducible Transcriptomics Analysis.
(Manuscript in preparation, Maastricht University).

R-ODAF (framework)
Verheijen M.C., Meier M.J., Asensio J.O., Gant T.W., Tong W., Yauk C.L., Caiment F. (2022).
R-ODAF: Omics Data Analysis Framework for Regulatory Application.
Regulatory Toxicology and Pharmacology, 131, 105143. https://doi.org/10.1016/j.yrtph.2022.105143

Credits

Developed by Saad Lodhi under supervision of Marcha Verheijen, Theo M. de Kok, Florian Caiment, and Danyel Jennen.
Part of the R-ODAF at Maastricht University.
