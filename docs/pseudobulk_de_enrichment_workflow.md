# Pseudobulk DE and pathway enrichment workflow

This workflow reproduces the manuscript pseudobulk differential-expression and pathway-enrichment source tables from the public processed Seurat object.

The public-facing script is:

```r
scripts/09_pseudobulk_de_enrichment.R
```

It performs, separately within term and preterm delivery classes:

1. subsetting to matched high/low MNP samples;
2. subject-aware pseudobulk aggregation by predicted cell type, MNP burden, and `study_id`;
3. limma high-versus-low MNP differential expression with duplicate-correlation support when repeated study-level blocks are present;
4. gene filtering to remove lncRNA/ENSG-like entries, Y-chromosome genes from enrichment, and pathogen-panel features;
5. GO Biological Process and MSigDB Hallmark over-representation analysis using `clusterProfiler`;
6. family-wise BH correction within delivery class, pathway family, cell type, and direction;
7. export of source tables for Fig. 3a-c and Supplementary Table S1e.

Original exploratory scripts are archived in `provenance/archived_scripts/`.

Important outputs:

```text
outputs/tables/pseudobulk_de/Fig3_pseudobulk_DE_ALL_delivery.tsv
outputs/tables/pseudobulk_de/Fig3a_pseudobulk_logFC_heatmap_source.tsv
outputs/tables/pseudobulk_de/Fig3_pathway_enrichment_ALL_delivery.tsv
outputs/tables/pseudobulk_de/Fig3bc_pathway_enrichment_plot_source.tsv
outputs/figures/pseudobulk_de/Fig3a_pseudobulk_logFC_heatmap_draft.png
outputs/figures/pseudobulk_de/Fig3bc_pathway_enrichment_draft.png
```

The generated PNGs are reviewer-facing draft reproductions. Final figure assembly may still use the manuscript-specific layout.
