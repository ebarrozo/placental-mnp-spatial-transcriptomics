## save_audra_macrophage_obj.R

rmd_obj <- readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/audra/mnp_obj_w_nn_predicted_celltype_fullcells_with_RNAassay.rds")

macrophages <- subset(
  rmd_obj,
  subset = (nn.predicted.celltype == "Macrophage" | nn.predicted.celltype == "Hofbauer Cell") &
    (Cumulative_MNPs_ug_per_g > 0)
)

DefaultAssay(macrophages) <- "RNA"
macrophages <- NormalizeData(macrophages, verbose = FALSE)
macrophages <- FindVariableFeatures(macrophages, verbose = FALSE)
macrophages_final <- ScaleData(macrophages, verbose = FALSE, vars.to.regress = c("platform"))
macrophages_final <- RunPCA(macrophages_final, verbose = FALSE)
macrophages_final <- FindNeighbors(macrophages_final, dims = 1:6)
macrophages_final <- FindClusters(macrophages_final, resolution = 0.5)
macrophages_final <- RunUMAP(macrophages_final, dims = 1:6)

macrophages_final$polarization <- dplyr::case_when(
  macrophages_final$seurat_clusters %in% c("0", "1", "2") ~ "M0",
  macrophages_final$seurat_clusters %in% c("4", "5", "14") ~ "M2",
  macrophages_final$seurat_clusters %in% c("6", "7", "8", "9", "10", "11", "12", "13") ~ "M1",
  macrophages_final$seurat_clusters == "3" ~ "Intermediate",
  TRUE ~ "Unassigned"
)

saveRDS(
  macrophages_final,
  "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/audra/macrophages_final_rebuilt_from_Rmd.rds"
)









library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(tibble)
library(pheatmap)
library(stringr)
library(forcats)

DefaultAssay(macrophages_final) <- "RNA"

# If not already normalized/scaled enough for plotting/scoring
macrophages_final <- NormalizeData(macrophages_final, verbose = FALSE)

mac_programs <- list(
  Inflammatory = c("IL1B","TNF","CXCL8","CXCL2","CCL2","CCL3","CCL4","NFKBIA","JUN","FOS"),
  Phagocytic_FcR = c("FCGR1A","FCGR2A","FCGR3A","MSR1","MARCO","LST1","TYROBP","CTSS","CTSB","LYZ"),
  Complement_Hofbauer = c("C1QA","C1QB","C1QC","APOE","TREM2","VSIG4","CD163","MRC1","FOLR2"),
  Vascular_ECM = c("SPP1","VEGFA","MMP9","MMP14","TIMP1","TIMP3","LGALS3","FN1","THBS1","SERPINE1"),
  Immunoregulatory = c("IL10","TGFB1","MERTK","MAF","SOCS3","CD274","HLA-DRA","HLA-DRB1","VSIR","LILRB1")
)

# score programs
old_meta <- colnames(macrophages_final@meta.data)

for (nm in names(mac_programs)) {
  genes_use <- intersect(mac_programs[[nm]], rownames(macrophages_final))
  if (length(genes_use) >= 2) {
    macrophages_final <- AddModuleScore(
      macrophages_final,
      features = list(genes_use),
      assay = "RNA",
      name = paste0(nm, "_score"),
      search = FALSE
    )
  } else {
    message("Too few genes present for: ", nm)
  }
}

# standardize AddModuleScore column names
new_cols <- setdiff(colnames(macrophages_final@meta.data), old_meta)

for (nm in names(mac_programs)) {
  cand <- grep(paste0("^", nm, "_score"), colnames(macrophages_final@meta.data), value = TRUE)
  if (length(cand) >= 1) {
    macrophages_final[[paste0(nm, "_score")]] <- macrophages_final@meta.data[[cand[1]]]
  }
}

score_cols <- paste0(names(mac_programs), "_score")

# optional z-score within macrophage object
for (sc in score_cols) {
  macrophages_final[[paste0(sc, "_z")]] <- as.numeric(scale(macrophages_final@meta.data[[sc]]))
}

score_cols_z <- paste0(score_cols, "_z")








mac_meta <- macrophages_final@meta.data %>%
  rownames_to_column("cell") %>%
  as_tibble()

# assign state by highest z-scored program
state_mat <- mac_meta %>%
  select(all_of(score_cols_z)) %>%
  as.matrix()

state_names <- gsub("_score_z$", "", colnames(state_mat))

max_state_idx <- apply(state_mat, 1, which.max)

mac_meta$mac_state <- state_names[max_state_idx]

# nicer labels
mac_meta$mac_state <- recode(
  mac_meta$mac_state,
  Inflammatory = "Inflammatory",
  Phagocytic_FcR = "Phagocytic / FcR",
  Complement_Hofbauer = "Complement / Hofbauer-like",
  Vascular_ECM = "Vascular / ECM-remodeling",
  Immunoregulatory = "Immunoregulatory"
)

macrophages_final$mac_state <- mac_meta$mac_state




mac_meta %>%
  dplyr::count(mac_state, sort = TRUE)

mac_meta %>%
  dplyr::count(seurat_clusters, mac_state) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::arrange(seurat_clusters, dplyr::desc(frac))
############################### inspect    
write.csv(mac_meta, "mac_state_counts.csv", row.names = FALSE)

p5a <- DimPlot(
  macrophages_final,
  reduction = "umap",
  group.by = "mac_state",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "Arial"),
    legend.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p5a


if (!"MNPBurden_Delivery_class" %in% colnames(macrophages_final@meta.data)) {
  macrophages_final$MNPBurden_Delivery_class <- paste(
    macrophages_final$MNPBurden,
    macrophages_final$Delivery_class,
    sep = "_"
  )
}

p5a_split <- DimPlot(
  macrophages_final,
  reduction = "umap",
  group.by = "mac_state",
  split.by = "MNPBurden_Delivery_class",
  label = FALSE,
  raster = FALSE,
  ncol = 2
) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "Arial"),
    legend.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p5a_split





state_frac_subject <- macrophages_final@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  tibble::as_tibble() %>%
  dplyr::distinct(cell, study_id, sample_id_std, Delivery_class, MNPBurden, mac_state) %>%
  dplyr::count(study_id, Delivery_class, MNPBurden, mac_state, name = "n_cells") %>%
  dplyr::group_by(study_id, Delivery_class, MNPBurden) %>%
  dplyr::mutate(frac = n_cells / sum(n_cells)) %>%
  dplyr::ungroup()
state_frac_subject %>%
  dplyr::arrange(Delivery_class, MNPBurden, study_id, mac_state) %>%
  print(n = Inf)

state_frac_stats <- state_frac_subject %>%
  dplyr::group_by(mac_state, Delivery_class) %>%
  dplyr::summarise(
    mean_low = mean(frac[MNPBurden == "Low"], na.rm = TRUE),
    mean_high = mean(frac[MNPBurden == "High"], na.rm = TRUE),
    estimate_high_minus_low = mean_high - mean_low,
    p_value = tryCatch(t.test(frac ~ MNPBurden)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  dplyr::mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(p_value)

state_frac_stats

state_frac_subject %>%
  dplyr::count(Delivery_class, MNPBurden, study_id, name = "n_states_present") %>%
  dplyr::arrange(Delivery_class, MNPBurden, study_id)
library(dplyr)
find("count")

readr::write_csv(
  state_frac_stats,
  "Fig5_macrophage_state_fraction_stats.csv"
)

p5b <- ggplot(
  state_frac_subject,
  aes(x = MNPBurden, y = frac)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
  geom_point(
    aes(group = study_id),
    position = position_jitter(width = 0.08, height = 0),
    size = 2.0,
    alpha = 0.9
  ) +
  facet_grid(mac_state ~ Delivery_class, scales = "free_y") +
  labs(
    x = "MNP burden",
    y = "Macrophage state fraction"
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 7.5, color = "black"),
    strip.text.x = element_text(face = "bold", color = "black"),
    legend.position = "none"
  )

p5b




genes_heat <- unique(c(
  "C1QA","C1QB","C1QC","APOE","TREM2","VSIG4","CD163","MRC1","FOLR2",
  "FCGR1A","FCGR2A","FCGR3A","CTSB","CTSS","MSR1","TYROBP","LST1",
  "SPP1","VEGFA","FN1","THBS1","SERPINE1","MMP14","TIMP3","LGALS3",
  "IL1B","TNF","CCL2","CXCL8","NFKBIA","JUN","FOS",
  "IL10","TGFB1","MERTK","SOCS3","CD274","HLA-DRA","VSIR"
))

genes_heat <- intersect(genes_heat, rownames(macrophages_final))

avg_expr_state <- AverageExpression(
  macrophages_final,
  assays = "RNA",
  features = genes_heat,
  group.by = "mac_state",
  slot = "data"
)$RNA

# z-score by gene
avg_expr_state_z <- t(scale(t(avg_expr_state)))

pheatmap(
  avg_expr_state_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize_row = 8,
  fontsize_col = 9,
  main = "Macrophage state-associated genes"
)


png("Fig5c_macrophage_state_gene_heatmap.png", width = 2200, height = 1800, res = 300)
pheatmap(
  avg_expr_state_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize_row = 8,
  fontsize_col = 9,
  main = "Macrophage state-associated genes"
)
dev.off()

genes_rep <- intersect(
  c("C1QA","FCGR3A","SPP1","FN1","IL1B","IL10","VEGFA","VSIG4"),
  rownames(macrophages_final)
)

p5d_genes <- DotPlot(
  macrophages_final,
  features = genes_rep,
  group.by = "mac_state",
  assay = "RNA"
) +
  RotatedAxis() +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

p5d_genes




subject_program_scores <- macrophages_final@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  tibble::as_tibble() %>%
  dplyr::select(
    cell,
    study_id,
    Delivery_class,
    MNPBurden,
    dplyr::all_of(score_cols_z)
  ) %>%
  dplyr::distinct() %>%
  dplyr::group_by(study_id, Delivery_class, MNPBurden) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(score_cols_z), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(score_cols_z),
    names_to = "program",
    values_to = "score"
  ) %>%
  dplyr::mutate(
    program = stringr::str_remove(program, "_score_z$"),
    program = dplyr::recode(
      program,
      Inflammatory = "Inflammatory",
      Phagocytic_FcR = "Phagocytic / FcR",
      Complement_Hofbauer = "Complement / Hofbauer-like",
      Vascular_ECM = "Vascular / ECM-remodeling",
      Immunoregulatory = "Immunoregulatory"
    )
  )

p5d_programs <- ggplot(
  subject_program_scores,
  aes(x = MNPBurden, y = score)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
  geom_point(
    aes(group = study_id),
    position = position_jitter(width = 0.08, height = 0),
    size = 2.0,
    alpha = 0.9
  ) +
  facet_grid(program ~ Delivery_class, scales = "free_y") +
  labs(
    x = "MNP burden",
    y = "Mean macrophage program score"
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 7.5, color = "black"),
    strip.text.x = element_text(face = "bold", color = "black"),
    legend.position = "none"
  )

p5d_programs



ggsave("Fig5a_macrophage_state_umap.png", p5a, width = 6.2, height = 4.8, dpi = 600)
ggsave("Fig5b_macrophage_state_fractions.png", p5b, width = 7.2, height = 6.5, dpi = 600)
ggsave("Fig5d_macrophage_program_scores.png", p5d_programs, width = 7.2, height = 6.5, dpi = 600)
ggsave("Fig5d_macrophage_representative_genes_dotplot.png", p5d_genes, width = 6.4, height = 4.8, dpi = 600)



polarization_qc <- macrophages_final@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  tibble::as_tibble() %>%
  dplyr::select(cell, study_id, Delivery_class, MNPBurden, seurat_clusters, polarization, mac_state) %>%
  dplyr::filter(!is.na(polarization), !is.na(mac_state))

polarization_qc %>%
  dplyr::count(polarization, sort = TRUE)

polarization_qc %>%
  dplyr::count(seurat_clusters, polarization, mac_state, name = "n") %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::arrange(seurat_clusters, dplyr::desc(frac)) %>%
  print(n = Inf)

readr::write_csv(
  polarization_qc,
  "Fig5_polarization_qc.csv"
)


p5a_polarization <- DimPlot(
  macrophages_final,
  reduction = "umap",
  group.by = "polarization",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) + 
  theme_classic(base_size = 11) +
  labs(title = NULL) +
  theme(
    text = element_text(family = "Arial"),
    legend.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p5a_polarization

p5b_compact <- ggplot(
  state_frac_subject,
  aes(x = MNPBurden, y = frac, group = study_id)
) +
  geom_line(alpha = 0.35, linewidth = 0.35) +
  geom_point(size = 2.1, alpha = 0.9) +
  stat_summary(
    aes(group = MNPBurden),
    fun = mean,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.35,
    color = "black"
  ) +
  facet_grid(mac_state ~ Delivery_class, scales = "free_y") +
  labs(
    x = "MNP burden",
    y = "Macrophage state fraction"
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 7.5, color = "black"),
    strip.text.x = element_text(face = "bold", color = "black"),
    legend.position = "none"
  )

p5b_compact


genes_rep <- intersect(
  c(
    "C1QA","C1QB","C1QC","APOE","VSIG4","CD163",
    "FCGR1A","FCGR3A","TYROBP","LST1","CTSS",
    "SPP1","FN1","LGALS3","MMP14","TIMP3",
    "IL1B","CCL2","CXCL8","NFKBIA",
    "MERTK","SOCS3","HLA-DRA","CD274"
  ),
  rownames(macrophages_final)
)

p5d_genes <- DotPlot(
  macrophages_final,
  features = genes_rep,
  group.by = "mac_state",
  assay = "RNA"
) +
  RotatedAxis() +
  theme_classic(base_size = 9) +
  labs(x = NULL, y = NULL) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

p5d_genes
p5c_genes_vertical <- DotPlot(
  macrophages_final,
  features = genes_rep,
  group.by = "mac_state",
  assay = "RNA"
) +
  coord_flip() +
  theme_classic(base_size = 9) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right"
  )

p5c_genes_vertical
ggsave(
  "Fig5c_macrophage_representative_genes_dotplot_vertical.png",
  p5c_genes_vertical,
  width = 7.2,
  height = 4.2,
  dpi = 600
)
ggsave("Fig5a_macrophage_state_umap_polarization.png", p5a_polarization, width = 6.2, height = 4.8, dpi = 600)
ggsave("Fig5b_macrophage_state_fractions_compact.png", p5b_compact, width = 7.2, height = 6.5, dpi = 600)
#ggsave("Fig5d_macrophage_program_scores.png", p5d_programs, width = 7.2, height = 6.5, dpi = 600)
ggsave("Fig5d_macrophage_representative_genes_dotplot.png", p5d_genes, width = 6.4, height = 4.8, dpi = 600)







genes_rep_ordered <- intersect(
  rev(c(
    # Complement / Hofbauer-like
    "C1QA","C1QB","C1QC","APOE","VSIG4","CD163",

    # Phagocytic / Fc receptor
    "FCGR1A","FCGR3A","TYROBP","LST1","CTSS",

    # Vascular / ECM-remodeling
    "SPP1","FN1","LGALS3","MMP14","TIMP3",

    # Inflammatory
    "IL1B","CCL2","CXCL8","NFKBIA",

    # Immunoregulatory
    "MERTK","SOCS3","HLA-DRA","CD274"
  )),
  rownames(macrophages_final)
)

macrophages_final$mac_state_short <- dplyr::recode(
  macrophages_final$mac_state,
  "Complement / Hofbauer-like" = "Complement",
  "Phagocytic / FcR" = "Phagocytic/FcR",
  "Vascular / ECM-remodeling" = "Vascular/ECM",
  "Inflammatory" = "Inflammatory",
  "Immunoregulatory" = "Immunoregulatory"
)

macrophages_final$mac_state_short <- factor(
  macrophages_final$mac_state_short,
  levels = c(
    "Complement",
    "Phagocytic/FcR",
    "Vascular/ECM",
    "Inflammatory",
    "Immunoregulatory"
  )
)

genes_rep_ordered_clean <- intersect(
  rev(c(
    "C1QA","APOE","VSIG4","CD163",
    "FCGR1A","FCGR3A","TYROBP","CTSS",
    "SPP1","FN1","MMP14","TIMP3",
    "IL1B","CCL2","NFKBIA",
    "MERTK","SOCS3","HLA-DRA"
  )),
  rownames(macrophages_final)
)


p5c_genes_vertical <- DotPlot(
  macrophages_final,
  features = genes_rep_ordered_clean,
  group.by = "mac_state_short",
  assay = "RNA"
) +
  coord_flip() +
  theme_classic(base_size = 9) +
  labs(x = NULL, y = NULL) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right"
  )

p5c_genes_vertical

p5c_genes_vertical
ggsave(
  "Fig5c_macrophage_representative_genes_dotplot_vertical.png",
  p5c_genes_vertical,
  width = 7.2,
  height = 4.2,
  dpi = 600
)


library(Seurat)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)

DefaultAssay(macrophages_final) <- "RNA"

# Make sure mac_state exists and is clean
macrophages_final$mac_state <- as.character(macrophages_final$mac_state)

Idents(macrophages_final) <- "mac_state"

mac_state_markers <- FindAllMarkers(
  macrophages_final,
  assay = "RNA",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25
) %>%
  as_tibble() %>%
  dplyr::arrange(cluster, p_val_adj, dplyr::desc(avg_log2FC))

readr::write_csv(
  mac_state_markers,
  "Fig5_macrophage_state_FindAllMarkers.csv"
)

mac_state_markers %>%
  dplyr::count(cluster)



top_mac_state_markers <- mac_state_markers %>%
  dplyr::filter(
    p_val_adj < 0.05,
    pct.1 >= 0.15
  ) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 3, with_ties = FALSE) %>%
  dplyr::ungroup()

top_mac_state_markers %>%
  dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  print(n = Inf)

readr::write_csv(
  top_mac_state_markers,
  "Fig5_macrophage_state_top_DEGs_for_dotplot.csv"
)




state_order <- c(
  "Complement / Hofbauer-like",
  "Phagocytic / FcR",
  "Vascular / ECM-remodeling",
  "Inflammatory",
  "Immunoregulatory"
)

top_mac_state_markers <- top_mac_state_markers %>%
  dplyr::mutate(cluster = factor(cluster, levels = state_order)) %>%
  dplyr::arrange(cluster, dplyr::desc(avg_log2FC))

genes_deg_ordered <- top_mac_state_markers %>%
  dplyr::pull(gene) %>%
  unique()

# coord_flip reverses the y-axis order, so reverse gene order for nicer display
genes_deg_ordered_flip <- rev(genes_deg_ordered)

macrophages_final$mac_state <- factor(
  macrophages_final$mac_state,
  levels = state_order
)



genes_deg_curated_from_markers <- c(
  # Complement / Hofbauer-like
  "VSIG4", "LYVE1", "MPEG1",

  # Phagocytic / FcR
  "FCGR3A", "FCGR2A", "MSR1",

  # Vascular / ECM-remodeling
  "FN1", "SPP1", "SERPINE1", "COL3A1",

  # Inflammatory or inflammatory-like
  "SPARC", "COL3A1", "STAB1",

  # Immunoregulatory
  "MAF", "LILRB1", "STAB1"
)



p5c_deg_dotplot <- DotPlot(
  macrophages_final,
  features = genes_deg_ordered_flip,
  group.by = "mac_state",
  assay = "RNA"
) +
  coord_flip() +
  theme_classic(base_size = 9) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right"
  )

p5c_deg_dotplot

ggsave(
  "Fig5c_macrophage_state_DEG_dotplot_vertical.png",
  p5c_deg_dotplot,
  width = 7.2,
  height = 4.8,
  dpi = 600
)












exclude_pattern <- "^(MT-|RPL|RPS|MALAT1|XIST|Y_RNA|HBB|HBA|HBM|HBD|HBE|HBQ)"

top_mac_state_markers_clean <- mac_state_markers %>%
  dplyr::filter(
    p_val_adj < 0.05,
    pct.1 >= 0.15,
    !grepl(exclude_pattern, gene, ignore.case = FALSE)
  ) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 3, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cluster = factor(cluster, levels = state_order)) %>%
  dplyr::arrange(cluster, dplyr::desc(avg_log2FC))

genes_deg_ordered_clean <- top_mac_state_markers_clean %>%
  dplyr::pull(gene) %>%
  unique()

p5c_deg_dotplot_clean <- DotPlot(
  macrophages_final,
  features = rev(genes_deg_ordered_clean),
  group.by = "mac_state",
  assay = "RNA"
) +
  coord_flip() +
  theme_classic(base_size = 9) +
  labs(x = NULL, y = NULL) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "right"
  )

p5c_deg_dotplot_clean







genes_deg_curated_from_markers <- c(
  # Complement / Hofbauer-like
  "VSIG4", "LYVE1", "MPEG1",

  # Phagocytic / FcR
  "FCGR3A", "FCGR2A", "MSR1",

  # Vascular / ECM-remodeling
  "FN1", "SPP1", "SERPINE1", "COL3A1",

  # Inflammatory or inflammatory-like
  "SPARC", "COL3A1", "STAB1",

  # Immunoregulatory
  "MAF", "LILRB1", "STAB1"
)

genes_deg_curated_from_markers <- unique(genes_deg_curated_from_markers)
genes_deg_curated_from_markers <- intersect(genes_deg_curated_from_markers, rownames(macrophages_final))
genes_deg_curated_from_markers <- unique(genes_deg_curated_from_markers)
genes_deg_curated_from_markers <- intersect(genes_deg_curated_from_markers, rownames(macrophages_final))
genes_deg_final <- intersect(
  rev(c(
    "VSIG4", "LYVE1", "MPEG1",
    "FCGR3A", "FCGR2A", "MSR1",
    "FN1", "SPP1", "SERPINE1", "SPARC", "COL3A1",
    "MAF", "LILRB1", "RYR1"
  )),
  rownames(macrophages_final)
)

p5c_deg_dotplot_final <- DotPlot(
  macrophages_final,
  features = genes_deg_final,
  group.by = "mac_state",
  assay = "RNA"
) +
  coord_flip() +
  theme_classic(base_size = 9) +
  labs(x = NULL, y = NULL) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "right"
  )

p5c_deg_dotplot_final

ggsave(
  "Fig5c_macrophage_state_DEG_dotplot_final.png",
  p5c_deg_dotplot_final,
  width = 7.0,
  height = 4.6,
  dpi = 600
)


macrophages_final$mac_state_plot <- dplyr::recode(
  macrophages_final$mac_state,
  "Inflammatory" = "Activated remodeling",
  "Complement / Hofbauer-like" = "Complement / Hofbauer-like",
  "Phagocytic / FcR" = "Phagocytic / FcR",
  "Vascular / ECM-remodeling" = "Vascular / ECM-remodeling",
  "Immunoregulatory" = "Immunoregulatory"
)

macrophages_final$mac_state_plot <- factor(
  macrophages_final$mac_state_plot,
  levels = c(
    "Complement / Hofbauer-like",
    "Phagocytic / FcR",
    "Vascular / ECM-remodeling",
    "Activated remodeling",
    "Immunoregulatory"
  )
)


genes_deg_final <- intersect(
  rev(c(
    "VSIG4", "LYVE1", "MPEG1",
    "FCGR3A", "FCGR2A", "MSR1",
    "FN1", "SPP1", "SERPINE1", "SPARC", "COL3A1",
    "MAF", "LILRB1"
  )),
  rownames(macrophages_final)
)


p5c_deg_dotplot_final <- DotPlot(
  macrophages_final,
  features = genes_deg_final,
  group.by = "mac_state_plot",
  assay = "RNA"
) +
  coord_flip() +
  theme_classic(base_size = 9) +
  labs(x = NULL, y = NULL) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

p5c_deg_dotplot_final


p5c_deg_dotplot_final

ggsave(
  "Fig5c_macrophage_state_DEG_dotplot_final.png",
  p5c_deg_dotplot_final,
  width = 7.0,
  height = 4.6,
  dpi = 600
)



macrophages_final$mac_state_plot <- dplyr::recode(
  as.character(macrophages_final$mac_state),
  "Inflammatory" = "Activated remodeling",
  "Complement / Hofbauer-like" = "Complement / Hofbauer-like",
  "Phagocytic / FcR" = "Phagocytic / FcR",
  "Vascular / ECM-remodeling" = "Vascular / ECM-remodeling",
  "Immunoregulatory" = "Immunoregulatory"
)

state_order <- c(
  "Complement / Hofbauer-like",
  "Phagocytic / FcR",
  "Vascular / ECM-remodeling",
  "Activated remodeling",
  "Immunoregulatory"
)

macrophages_final$mac_state_plot <- factor(
  macrophages_final$mac_state_plot,
  levels = state_order
)


state_frac_subject <- macrophages_final@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  tibble::as_tibble() %>%
  dplyr::select(cell, study_id, sample_id_std, Delivery_class, MNPBurden, mac_state_plot) %>%
  dplyr::distinct() %>%
  dplyr::count(study_id, Delivery_class, MNPBurden, mac_state_plot, name = "n_cells") %>%
  dplyr::group_by(study_id, Delivery_class, MNPBurden) %>%
  dplyr::mutate(frac = n_cells / sum(n_cells)) %>%
  dplyr::ungroup()

state_frac_stats <- state_frac_subject %>%
  dplyr::group_by(mac_state_plot, Delivery_class) %>%
  dplyr::summarise(
    n_subjects = dplyr::n_distinct(study_id),
    n_low = dplyr::n_distinct(study_id[MNPBurden == "Low"]),
    n_high = dplyr::n_distinct(study_id[MNPBurden == "High"]),
    mean_low = mean(frac[MNPBurden == "Low"], na.rm = TRUE),
    mean_high = mean(frac[MNPBurden == "High"], na.rm = TRUE),
    estimate_high_minus_low = mean_high - mean_low,
    p_value = ifelse(
      n_low >= 2 & n_high >= 2,
      tryCatch(t.test(frac ~ MNPBurden)$p.value, error = function(e) NA_real_),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(Delivery_class, p_value)

state_frac_stats

readr::write_csv(
  state_frac_stats,
  "Fig5b_macrophage_state_fraction_stats_activated_remodeling.csv"
)

p5b_state_fractions <- ggplot(
  state_frac_subject,
  aes(x = MNPBurden, y = frac)
) +
  geom_point(
    aes(group = study_id),
    position = position_jitter(width = 0.08, height = 0),
    size = 2.1,
    alpha = 0.9
  ) +
  stat_summary(
    aes(group = MNPBurden),
    fun = mean,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.35,
    color = "black"
  ) +
  facet_grid(mac_state_plot ~ Delivery_class, scales = "free_y") +
  labs(
    x = "MNP burden",
    y = "Macrophage state fraction"
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 7.5, color = "black"),
    strip.text.x = element_text(face = "bold", color = "black"),
    legend.position = "none"
  )

p5b_state_fractions

ggsave(
  "Fig5b_macrophage_state_fractions_activated_remodeling.png",
  p5b_state_fractions,
  width = 7.2,
  height = 6.2,
  dpi = 600
)





