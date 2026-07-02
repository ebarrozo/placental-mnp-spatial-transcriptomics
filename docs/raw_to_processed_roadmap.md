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
