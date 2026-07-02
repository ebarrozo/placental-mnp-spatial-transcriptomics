# Preprocessing provenance overview

This directory documents the upstream preprocessing and integration workflow used to create the public analysis object.

## Scope

These scripts and notes are provided for transparency and provenance. They are not required to regenerate the manuscript figures and tables from the public processed object.

## Major steps

1. Platform-specific object creation from Space Ranger and Xenium outputs.
2. Quality control and metadata harmonization.
3. Reference-based cell-type prediction / RCTD annotation.
4. Batch execution because full integration exceeded server memory.
5. Conversion to SCT assays to stabilize integration.
6. Cross-platform object integration.
7. Final forced annotation object, `obj_forced`.
8. Public object cleanup using `scripts/01_build_public_object.R`.

## Archived scripts

Place sanitized historical scripts in `provenance/archived_scripts/` only when they clarify the provenance of the final object. Do not put abandoned exploratory code in the primary workflow.

## v10 upstream preprocessing additions

The upstream provenance now includes:

- archived Visium v1 Space Ranger scripts (`run_count*.sh`, `visiumv1_spaceranger_v1.sh`);
- archived Visium HD Space Ranger script (`visiumhd_spaceranger_v1.sh`);
- archived placental reference construction script (`meta_placenta_ref_v1.R`);
- a parsed Space Ranger manifest (`config/spaceranger_count_manifest_from_archived_scripts.csv`);
- a reusable dry-run wrapper for Space Ranger commands;
- a public-facing placental-reference build scaffold from h5ad inputs.

These files are intended to document the raw-to-processed path. Exact reproduction requires raw FASTQs, histology/CytAssist images, manual Loupe alignment JSONs, and public h5ad reference files.

## RCTD, chunk processing, harmonization, and final integration provenance

The repository now includes archived upstream scripts documenting the RCTD annotation-transfer and high-memory SCT chunking/integration workflow. These scripts are retained under:

```text
provenance/archived_scripts/upstream/rctd_chunk_integration/
```

They are not part of the default public analysis run because they require large GEO-linked raw/intermediate files and server-specific paths. See:

```text
docs/rctd_chunk_integration_workflow.md
workflows/02_rctd_annotation/README.md
workflows/03_sct_chunking_integration/README.md
config/upstream_rctd_chunk_integration_script_inventory.csv
```
