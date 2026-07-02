# Raw-to-processed roadmap

This document summarizes how raw/platform-native data were transformed into the public processed object.

1. Raw/platform-native data deposited in GEO.
2. Platform-specific preprocessing:
   - Visium v1: Space Ranger outputs.
   - Visium HD: Space Ranger / Visium HD outputs.
   - Xenium 500-gene panel: Xenium output bundle.
   - Xenium 5,000-gene panel: Xenium output bundle with segmentation.
3. Reference-based annotation using placenta sc/snRNA-seq reference and RCTD/spacexr.
4. Batch execution to avoid memory crashes.
5. SCT assay conversion for integration stability.
6. Harmonized metadata and final public object generation using `scripts/01_build_public_object.R`.

The default manuscript reproducibility workflow does not require rerunning these high-memory preprocessing steps.

## Upstream preprocessing layer added in v10

The upstream preprocessing layer now includes archived provenance and public-facing scaffolds for:

1. **Placental reference construction** from public h5ad files used for RCTD/annotation transfer.
2. **Space Ranger count preprocessing** for Visium v1 and Visium HD runs.

Key files:

```text
workflows/00_reference_preprocessing/README.md
workflows/01_spaceranger/README.md
scripts/upstream_preprocessing/00_build_placental_reference_from_h5ad.R
scripts/upstream_preprocessing/run_spaceranger_count_from_manifest.sh
config/placental_reference_h5ad_manifest.csv
config/spaceranger_count_manifest_from_archived_scripts.csv
provenance/archived_scripts/upstream/
```

The archived scripts preserve server-specific paths and exact historical commands. The public-facing scripts are intended as reproducible scaffolds that require local path edits and access to raw/public h5ad, FASTQ, image, and Loupe alignment files.

Future upstream additions should include:

- custom reference build commands;
- RCTD training and annotation transfer;
- Xenium / Visium object import;
- chunked SCT/batch integration;
- final public-object construction handoff.

## Added upstream provenance: RCTD and chunk integration

The raw-to-processed roadmap now includes the following additional upstream stages:

1. RCTD annotation transfer from a placental reference using `spacexr`.
2. Barcode/sample harmonization across Visium v1, Xenium, and Visium HD.
3. Memory-safe SCT chunking and mapping.
4. Lightweight UMAP/metadata reconstruction from mapped chunks.
5. Final Seurat v5-safe platform merge and integrated object generation.

Historical scripts for these stages are archived as provenance. Exact re-execution requires GEO-linked raw/intermediate files.
