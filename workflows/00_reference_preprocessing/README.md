# Placental reference preprocessing

This directory documents the upstream construction of the placental reference used for RCTD/annotation transfer.

The original server script is archived as:

```text
provenance/archived_scripts/upstream/meta_placenta_ref_v1.R
```

The public-facing scaffold is:

```text
scripts/upstream_preprocessing/00_build_placental_reference_from_h5ad.R
```

The h5ad input manifest template is:

```text
config/placental_reference_h5ad_manifest.csv
```

The archived workflow downloaded public h5ad files, converted each to a Seurat object with `zellkonverter::readH5AD()`, retained public metadata, merged selected placental reference datasets, harmonized gene names and cell-type labels, and saved a large placental reference Seurat object for downstream annotation.

The full reference object is too large for GitHub. Reproduction requires downloading the public h5ad files listed in the manifest, then running:

```bash
Rscript scripts/upstream_preprocessing/00_build_placental_reference_from_h5ad.R
```

Downstream RCTD and chunk/batch processing scripts will be added in the next upstream-preprocessing layer.
