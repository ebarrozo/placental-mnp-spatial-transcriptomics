# Workflow 03: SCT chunking, harmonization, integration, and final UMAP

This folder documents the memory-safe SCT chunking and final integration stage used to generate the cross-platform public analysis object.

## Purpose

1. Load RCTD-annotated source Seurat objects.
2. Select the appropriate source assay for each platform.
3. Split large objects into manageable chunks.
4. Run `SCTransform` per chunk while retaining all genes needed for downstream mapping.
5. Map chunks into a shared reference/integration space.
6. Rebuild lightweight visualization objects from mapped chunk embeddings and metadata.
7. Harmonize sample, chunk, study, platform, and annotation metadata.
8. Merge platform-specific SCT objects into the final integrated object and UMAP.

## Key archived scripts

| Script | Role |
|---|---|
| `sctransform_spatial_harmonized.R` | Main chunking/SCT/mapping workflow; drops images and heavy slots for memory control. |
| `harmonized_vis_v1.R` | Builds a lightweight visualization skeleton from mapped chunk embeddings. |
| `viz_meta_remap.R` | Remaps chunk/sample/platform metadata onto the visualization object. |
| `rebuild_viz_embed_only.R` | Reconstructs an embedding-only visualization object from mapped chunks. |
| `obj_integrated_final.R` | Final Seurat v5-safe platform merge and dimensional-reduction workflow. |

## Important constraints

This workflow was developed for high-memory servers and large RDS files. The public repository therefore keeps these scripts as provenance rather than making them part of the default `00_run_all.R` path. Reviewers should start from the public processed object unless they specifically need to audit the full raw-to-processed path.
