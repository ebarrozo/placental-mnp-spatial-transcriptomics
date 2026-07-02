# Xenium Explorer annotation exports

This directory is populated by `scripts/08_xenium_explorer_exports.R`.

The manuscript Fig. 6 panels were generated manually in Xenium Explorer 4 using the 10x Genomics Xenium output directories listed in `config/figure6_xenium_explorer_manifest.csv`. Because the complete Xenium output bundles are too large for GitHub, this repository provides the sample/panel mapping and reproducible CSV exports of cell-level annotations that can be imported into Xenium Explorer.

Expected generated files include:

- `Fig6_XeniumExplorer_panel_manifest_resolved.csv`
- `<sample_id_std>_predicted_celltype_forXeniumExplorer.csv`
- `<sample_id_std>_nn_predicted_celltype_forXeniumExplorer.csv`
- `<sample_id_std>_celltype_display_forXeniumExplorer.csv`, when available
- `<sample_id_std>_full_public_metadata_forXeniumExplorer.csv`

Each annotation-import CSV contains two columns: `cell_id` and `group`.
