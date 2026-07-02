# Communication-analysis workflow

The manuscript communication analyses used NicheNet and CellChat to prioritize ligand-receptor and pathway-level communication programs associated with placental MNP burden. The original exploratory scripts are preserved under `provenance/archived_scripts/`.

For the public reproducibility workflow, communication figure generation is source-table-first. This avoids requiring reviewers to re-run high-memory CellChat/NicheNet jobs, download external priors, or match package-version-specific APIs. Place the curated exported CSVs in:

```text
data/derived/communication/
```

Expected files:

```text
nichenet_ligand_network_edges.csv
nichenet_top_ligands.csv
cellchat_core_pathway_sender_contributions.csv
```

Then run:

```r
source("scripts/07_communication_figures.R")
```

The script writes a manifest of expected source files and generates available communication tables/figures under:

```text
outputs/tables/communication/
outputs/figures/communication/
```

Full re-execution from the public Seurat object is intentionally not part of `00_run_all.R` because the original communication analysis required very large spatial objects, SCT data, NicheNet priors, CellChat databases, and CellPhoneDB/Python-side dependencies. Those scripts are retained as provenance rather than as the default reviewer-facing workflow.
