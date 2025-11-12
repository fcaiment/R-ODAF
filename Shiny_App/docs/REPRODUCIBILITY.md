# Reproducibility and Analytical Criteria

The **R-ODAF-Shiny** interface strictly follows the regulatory analytical rules established in the R-ODAF framework.  
All thresholds and filters are fixed to prevent user-induced variability and ensure repeatable results across sessions.

### Parameters
| Criterion | Value / Rule | Description |
|------------|---------------|-------------|
| Minimum coverage | 5 million reads per sample | Ensures adequate sequencing depth |
| CPM threshold | 1 | Filters out very lowly expressed genes |
| FDR cutoff | 0.01 | Controls false discovery rate in DESeq2 results |
| Presence filter | ≥ 75 % of samples | Removes sparsely detected genes |
| PCA outlier threshold | 20 % variance | Excludes outlier samples prior to DEG analysis |

### Provenance and environment control
- `renv.lock` tracks exact R and package versions.  
- `docs/sessionInfo.txt` captures the runtime environment.  
- Each Shiny session generates internal provenance logs for timestamp, thresholds, and result checksums.

### Validation
This interface reproduces identical results to the R-ODAF command-line workflow for identical input datasets.  
Consistency verified using test datasets from TransQST and internal validation runs.
