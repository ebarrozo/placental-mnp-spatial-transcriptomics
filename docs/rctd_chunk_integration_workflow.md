# RCTD, chunk processing, harmonization, and final integration workflow

This document summarizes the upstream scripts used to move from raw/platform-specific spatial objects to the integrated Seurat object used by the public analysis scripts.

## Workflow map

1. **Space Ranger / Xenium output generation**  
   Generated outside R using the archived Space Ranger shell scripts and native Xenium outputs.

2. **Placental reference construction**  
   `meta_placenta_ref_v1.R` builds and harmonizes a placental reference object from public h5ad files and existing placental/decidual references.

3. **RCTD annotation transfer**  
   The `RCTD_spatialmnp*.R` scripts load spatial Seurat objects, extract counts and coordinates, build `SpatialRNA` queries, build a `spacexr::Reference`, and run RCTD by platform/sample group.

4. **Annotation harmonization**  
   `weee_harmony.R` and related scripts resolve sample names, duplicate sample IDs, barcode formats, and predicted annotation columns.

5. **SCT chunking and mapping**  
   `sctransform_spatial_harmonized.R` chunks large objects, applies SCTransform, drops heavy image/reduction slots, and maps chunks into a shared low-dimensional representation.

6. **Visualization object rebuild**  
   `harmonized_vis_v1.R`, `viz_meta_remap.R`, and `rebuild_viz_embed_only.R` reconstruct lightweight metadata/embedding objects for UMAP visualization and QC.

7. **Final integration**  
   `obj_integrated_final.R` performs final Seurat v5-safe SCT merges by platform and writes the integrated object used for downstream annotation, composition, module-score, macrophage-state, communication, and supplementary-table workflows.

## Why this is provenance-only

The upstream path uses large raw/intermediate data: Space Ranger outputs, Xenium outputs, image-bearing Seurat objects, RCTD result objects, per-chunk SCT RDS files, and mapped chunk files. These are expected to be linked through GEO rather than stored in GitHub. The GitHub repo preserves the exact historical scripts and provides a public processed object workflow for routine reproducibility.

## Script inventory

See:

```text
config/upstream_rctd_chunk_integration_script_inventory.csv
```
