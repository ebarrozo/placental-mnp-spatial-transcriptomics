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
