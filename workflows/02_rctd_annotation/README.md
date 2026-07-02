# Workflow 02: RCTD annotation transfer

This folder documents the RCTD annotation-transfer stage used upstream of the public analysis object. The full runs depend on large Space Ranger/Xenium/Visium HD outputs, the placental reference object, and intermediate Seurat objects that are not stored in GitHub. The historical scripts are preserved under:

```text
provenance/archived_scripts/upstream/rctd_chunk_integration/
```

## Purpose

1. Load spatial Seurat objects from Visium v1, Xenium, and Visium HD.
2. Extract count matrices and spatial coordinates in a platform-aware way.
3. Build `spacexr::SpatialRNA` query objects.
4. Build an RCTD reference from the placental reference object using `cell_type_expanded` labels.
5. Run RCTD by platform/sample group.
6. Harmonize barcodes and write predicted annotation columns back to Seurat objects.

## Key archived scripts

| Script | Role |
|---|---|
| `RCTD_spatialmnp.R` | Broad RCTD run for Visium v1/v2 and initial post-processing. |
| `RCTD_spatialmnp_visiumv1.R` | Revised Visium v1 RCTD run with dynamic UMI filtering and gene-overlap checks. |
| `RCTD_spatialmnp_xen_hd.R` | Xenium and Visium HD RCTD run with dynamic UMI filtering. |
| `RCTD_spatialmnp_xen_hd_HP13b_v1.R` | Focused HP13b / Visium v1 rescue run. |
| `HP2_hd.R` | Special handling for duplicate `HP2` identifiers across Xenium and Visium HD. |
| `weee_harmony.R` | Annotation metadata remapping back into source Seurat objects. |

## Public reproducibility note

The GitHub repository records the workflow and scripts, but exact re-execution requires the GEO-linked raw/intermediate data, including Space Ranger outputs, Xenium outputs, image-bearing Seurat objects, RCTD result RDS files, and the placental reference object. The downstream public object and source tables are the reproducible analysis entry point for manuscript figures.
