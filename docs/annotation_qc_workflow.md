# Annotation QC and origin-support workflow

This workflow documents the annotation-support analyses used for Supplementary Figures 1-3.

## What can be regenerated from the public object

The public object retains integrated metadata columns used for reviewer-facing annotation QC, including cell-type annotations, platform/sample identifiers, forced-annotation support columns when present, and rule-based maternal/fetal origin assignments. The public workflow therefore regenerates:

- Supplementary Fig. 1b-style UMAP split by platform.
- Supplementary Fig. 1c-d-style platform composition and forced-label support summaries.
- Supplementary Fig. 2a-b-style maternal/fetal origin UMAP and origin-by-annotation summaries.
- Supplementary Fig. 3a-b-style sample/platform annotation support summaries.

Run:

```r
source("scripts/02_load_public_object.R")
source("scripts/10_annotation_qc_support.R")
```

Outputs are written to:

```text
outputs/tables/annotation_qc/
outputs/figures/annotation_qc/
```

## What is provenance-only

Supplementary Fig. 1a required the original image-bearing Seurat objects, because images and the RNA assay were dropped during batch/chunk processing. The repository therefore archives the original image-remap script as provenance rather than claiming that SF1a can be regenerated from the public object alone.

The archived scripts also document how the final public metadata was produced:

- `make_obj_with_RNA_assay_fullcells.R`: reconstructs an RNA assay from SCT data for downstream compatibility. These RNA counts are pseudo-counts derived from SCT log-normalized data, not original raw RNA counts.
- `force_annot.R`: performs KNN-based rescue of low-confidence RCTD labels using UMAP neighborhoods and stores support values.
- `maternal_fetal_origin_v2.R`: assigns conservative rule-based maternal/fetal origin labels and adds sex-linked XIST/Y support where available.
- `obj_image_remap_v2_for_annotation_remap.R`: remaps annotations onto original image-bearing sample objects and exports annotation CSVs for visualization in platform-native viewers.
