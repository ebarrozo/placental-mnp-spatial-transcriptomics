# Reproducibility notes

The primary workflow starts from `spatial_mnp_publication_processed_seurat_v1.rds`, the public processed Seurat object. This is intentional: upstream RCTD annotation, SCT conversion, and cross-platform integration required high-memory batch processing and are documented under `provenance/` rather than required for routine figure/table regeneration.

Recommended reviewer/public command:

```r
install.packages("renv")
renv::restore()
source("scripts/00_run_all.R")
```

Generated outputs are written to `outputs/figures`, `outputs/tables`, and `outputs/logs`.
