# Figure 6 Xenium Explorer workflow

Fig. 6 was generated in Xenium Explorer 4 from native 10x Genomics Xenium output directories, not directly from the public Seurat object. The public Seurat object is sufficient to reproduce the annotation-export CSVs and sample mapping, but it is not sufficient to regenerate the exact Explorer screenshots because Xenium Explorer requires the corresponding platform-native output bundle.

## Panel mapping

| Panel | Delivery class | MNP burden | Marker | Public sample ID | Xenium output directory | Region label |
|---|---|---|---|---|---|---|
| A | Preterm | Low | HOXB7 | LP1_xen | `output-XETG00320__0023238__Region_1__20250108__200831` | Belfort1-region1 |
| B | Preterm | High | HOXB7 | HP1_xen | `output-XETG00320__0023238__Region_2__20250108__200832` | Belfort1-region2 |
| C | Term | Low | PECAM1 | LT2_xen | `output-XETG00320__0058012__Region_2__20250116__182319` | Belfort2-region2 |
| D | Term | High | PECAM1 | HT2b_xen | `output-XETG00320__0058014__Region_2__20250116__182319` | Belfort2-region2 |

## Reproducible repository component

Run:

```r
source("scripts/02_load_public_object.R")
source("scripts/08_xenium_explorer_exports.R")
```

This writes Xenium Explorer-compatible annotation CSVs to:

```text
data/derived/xenium_explorer_annotations/
```

Each import file has the two-column format expected for custom cell-group import:

```text
cell_id,group
```

The script removes trailing Seurat suffixes such as `_1` from cell names to better match Xenium Explorer cell IDs, following the original provenance script.

## Manual Explorer reconstruction

For each panel:

1. Open the listed Xenium output directory in Xenium Explorer 4.
2. Load the corresponding annotation CSV if cell-type or niche overlays are needed.
3. Visualize the marker listed for the panel: HOXB7 for preterm panels and PECAM1 for term panels.
4. Use the region/sample mapping above to identify the correct panel.
5. Export the final view from Xenium Explorer.

## Provenance

The original notebook-style script that built Seurat niche assays and exported Xenium Explorer annotations is archived at:

```text
provenance/archived_scripts/obj_image_remap_v2.R
```

The original script used `BuildNicheAssay()`, renamed Seurat cell barcodes by removing trailing suffixes, and exported `predicted.celltype` and `niches` CSVs for Xenium Explorer import.


Note: The manuscript Fig. 6 marker panels do not use niche overlays. Niche export code is retained only in archived provenance scripts from earlier exploratory spatial visualizations.
