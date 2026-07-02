# Spatial MNP placenta analysis code

This repository contains code to reproduce the downstream analyses, figures, and tables for the manuscript:

**Placental Micro- and Nanoplastic Burden is Associated with Spatial Immune–Vascular Endotypes in Term and Preterm Birth**

## Repository scope

The primary reproducibility workflow starts from the public processed Seurat object:

```text
data/external/spatial_mnp_publication_processed_seurat_v1.rds
```

This object is generated from the internal `obj_forced` object using:

```text
scripts/01_build_public_object.R
```

Upstream preprocessing, RCTD annotation, SCT conversion, and cross-platform integration required high-memory batch execution. Those steps are documented for provenance in `provenance/`, but are not required to regenerate manuscript-level outputs from the public processed object.

## Data availability

Raw sequencing data, platform-native outputs, images, metadata needed for manual alignment, source tables, and the public processed Seurat object will be deposited in GEO under accession `GSEXXXXXX`.

## Quick start

1. Clone this repository.
2. Place the public processed object at:

```text
data/external/spatial_mnp_publication_processed_seurat_v1.rds
```

3. Restore the R environment:

```r
install.packages("renv")
renv::restore()
```

4. Run the base workflow:

```r
source("scripts/00_run_all.R")
```

Generated outputs are written to:

```text
outputs/figures/
outputs/tables/
outputs/logs/
```

## Current workflow status

The starter workflow currently includes:

- public object generation from `obj_forced`;
- public object metadata validation;
- model-ready sample/run index generation;
- draft subject-level cell-type composition summaries;
- draft subject-level platform-normalized module score summaries;
- a figure/table manifest for tracking remaining figure-specific code cleanup.

Figure-specific scripts should be added incrementally as final working analysis code is cleaned and validated.

## Recommended development pattern

For each figure or table:

1. Identify the final working code block from historical scripts.
2. Move reusable logic to `R/` helper files.
3. Create or update one figure/table script in `scripts/`.
4. Export both the final figure and its source table.
5. Update `config/figure_manifest.csv`.
6. Confirm the figure/table can be regenerated from `spatial_mnp_publication_processed_seurat_v1.rds`.

## Citation

Please cite the associated manuscript and archived software release DOI.

## Migration from working analysis scripts

The original notebook-like working script `meta_update_2026-06-12.R` is archived under `provenance/archived_scripts/`. The public workflow extracts the downstream, reviewer-facing analyses into modular scripts. The first migrated analyses are cell-type composition and platform-normalized module-score analyses.

Run all currently migrated downstream analyses with:

```r
source("scripts/00_run_all.R")
```

For faster testing of the permutation section in `scripts/05_module_scores.R`, set the number of permutations before running:

```r
Sys.setenv(N_PERM = 100)
source("scripts/05_module_scores.R")
```

## Communication analyses

NicheNet, CellChat, and CellPhoneDB exploratory/rebuild scripts are archived in `provenance/archived_scripts/`. The public workflow uses curated exported source tables in `data/derived/communication/` to regenerate communication figures without requiring reviewers to re-run high-memory ligand-receptor jobs. See `docs/communication_workflow.md`.


## Fig. 6 / Xenium Explorer note

Figure 6 was generated in Xenium Explorer 4 from native 10x Genomics Xenium output directories. The repository provides `config/figure6_xenium_explorer_manifest.csv`, `docs/figure6_xenium_explorer_workflow.md`, and `scripts/08_xenium_explorer_exports.R` to reproduce the sample mapping and annotation CSVs needed for Explorer import. The exact Explorer screenshots require the native Xenium output bundles listed in the manifest.

## Optional pseudobulk DE and pathway enrichment

The pseudobulk DE/pathway-enrichment workflow is available as an opt-in step because it requires Bioconductor annotation/enrichment packages and can be slower than the figure/table scripts that only summarize existing metadata.

```r
Sys.setenv(RUN_PSEUDOBULK = "true")
source("scripts/00_run_all.R")
```

Or run directly after loading the public object:

```r
source("scripts/02_load_public_object.R")
source("scripts/09_pseudobulk_de_enrichment.R")
```

Outputs are written to `outputs/tables/pseudobulk_de/` and `outputs/figures/pseudobulk_de/`.


## Exploratory polymer-specific module associations

The repository includes `scripts/07_polymer_associations.R`, which regenerates exploratory polymer-module slope tables, leave-one-subject-out summaries, and Supplementary Fig. 4c-d style heatmap/scatter outputs from the public Seurat object and subject-level platform-normalized module scores. These analyses are hypothesis-generating and should be interpreted as prioritization rather than definitive polymer-specific mechanism testing.

## Annotation QC and origin support

The repository includes a public-object workflow for Supplementary Fig. 1b-d, Supplementary Fig. 2a-b, and Supplementary Fig. 3-style annotation support summaries:

```r
source("scripts/02_load_public_object.R")
source("scripts/10_annotation_qc_support.R")
```

Supplementary Fig. 1a is documented as provenance-only because it requires the original image-bearing Seurat sample objects that were not retained in the batch/chunk-processed public object.


## Supplementary Table 1 workbook

The repository includes a public-facing Supplementary Table 1 builder:

```r
source("scripts/11_build_supplementary_table1.R")
```

This script reads migrated output source tables from `outputs/tables/` and writes:

```text
outputs/tables/supplementary_table1/Supplementary_Table_1_MNP_spatial_analysis_outputs.xlsx
outputs/tables/supplementary_table1/csv_source_tables/
```

The original working script is archived as `provenance/archived_scripts/supp_table1.R`. The public builder intentionally avoids relying on objects left in the R session after `meta_update.R` and uses a conservative metadata allowlist so direct identifiers, exact dates, free-text clinical fields, and PHI-like columns are not included in public metadata exports.

## Upstream preprocessing provenance

This repository now includes an upstream preprocessing layer for the raw-to-public-object workflow. Server-specific historical scripts are archived under `provenance/archived_scripts/upstream/`, while public-facing scaffolds and manifests are provided under `scripts/upstream_preprocessing/`, `workflows/00_reference_preprocessing/`, and `workflows/01_spaceranger/`.

These upstream scripts document placental reference construction from public h5ad files and Space Ranger count preprocessing for Visium v1 and Visium HD samples. They require local edits to raw-data paths before execution.


## Upstream RCTD and chunk-integration provenance

The repository includes an upstream provenance layer for RCTD annotation transfer, sample/barcode harmonization, SCT chunking, mapped-chunk integration, and final UMAP/object construction. These scripts are archived under:

```text
provenance/archived_scripts/upstream/rctd_chunk_integration/
```

The full upstream workflow requires large GEO-linked files and server-scale memory. The default public analysis starts from the processed public Seurat object; the upstream scripts document how that object was generated and how RCTD/chunk-processing decisions were made.

See:

```text
docs/rctd_chunk_integration_workflow.md
workflows/02_rctd_annotation/README.md
workflows/03_sct_chunking_integration/README.md
config/upstream_rctd_chunk_integration_script_inventory.csv
```

## Path sanitization

Private absolute server paths in archived provenance scripts have been replaced with generic placeholders before public release. See `docs/path_sanitization_note.md` for details.

