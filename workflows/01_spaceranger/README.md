# Space Ranger preprocessing

This directory documents the Visium v1 and Visium HD `spaceranger count` preprocessing used upstream of the public Seurat object.

The archived server scripts are stored in:

```text
provenance/archived_scripts/upstream/
```

The parsed command manifest is:

```text
config/spaceranger_count_manifest_from_archived_scripts.csv
```

Run commands from the manifest with:

```bash
DRY_RUN=true bash scripts/upstream_preprocessing/run_spaceranger_count_from_manifest.sh
```

To run selected samples after editing local paths:

```bash
DRY_RUN=false RUN_IDS=S11231,S13627 WORKDIR=/path/to/results/spaceranger \
  SPACERANGER_BIN=/path/to/spaceranger \
  bash scripts/upstream_preprocessing/run_spaceranger_count_from_manifest.sh
```

Notes:

- Visium v1 runs used Space Ranger 3.1.3 and a custom reference that included GRCh38 plus pathogen genomes.
- Visium HD runs used Space Ranger 4.0.1, GRCh38-2024-A, and the Visium Human Transcriptome Probe Set v2.0.
- Manual Loupe alignment JSON files are required to reproduce the exact spatial registration.
- BAM creation differed by run; this is preserved in the manifest.
- Server-specific paths in archived scripts are provenance and should be edited before reuse.
