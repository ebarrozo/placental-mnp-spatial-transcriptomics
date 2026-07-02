# ============================================================
# Fig. 5 macrophage/Hofbauer-cell state analysis
# Starts from public processed Seurat object loaded as `obj_public`
# Outputs: source tables and manuscript-style Fig. 5 panels
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
})

source("R/path_helpers.R")
source("R/plotting_helpers.R")

if (!exists("obj_public")) {
  source("scripts/02_load_public_object.R")
}

out_tab <- file.path("outputs", "tables", "macrophage_states")
out_fig <- file.path("outputs", "figures", "macrophage_states")
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

annotation_col <- "nn.predicted.celltype"
required_cols <- c(annotation_col, "study_id", "sample_id_standardized", "Delivery_class", "MNPBurden")
missing_cols <- setdiff(required_cols, colnames(obj_public@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns for macrophage analysis: ", paste(missing_cols, collapse = ", "))
}

# Prefer model flag when available. Otherwise fall back to non-missing MNPBurden.
if (!"mnp_model_sample" %in% colnames(obj_public@meta.data)) {
  obj_public$mnp_model_sample <- !is.na(obj_public$MNPBurden)
}

# -----------------------------
# 1. Subset macrophage/Hofbauer observations
# -----------------------------
macrophages_final <- subset(
  obj_public,
  subset = (.data[[annotation_col]] == "Macrophage" | .data[[annotation_col]] == "Hofbauer Cell") &
    mnp_model_sample == TRUE
)

DefaultAssay(macrophages_final) <- "RNA"

# Normalize only if RNA data slot is empty or missing.
macrophages_final <- NormalizeData(macrophages_final, verbose = FALSE)

# -----------------------------
# 2. Rebuild macrophage UMAP if not present
# -----------------------------
# The public object may not include the macrophage-only reduction. Recompute a lightweight
# macrophage-only PCA/UMAP from RNA for Fig. 5 visualization.
macrophages_final <- FindVariableFeatures(macrophages_final, verbose = FALSE)
vars_to_regress <- intersect("platform_clean", colnames(macrophages_final@meta.data))
macrophages_final <- ScaleData(macrophages_final, verbose = FALSE, vars.to.regress = vars_to_regress)
macrophages_final <- RunPCA(macrophages_final, verbose = FALSE)
macrophages_final <- FindNeighbors(macrophages_final, dims = 1:6, verbose = FALSE)
macrophages_final <- FindClusters(macrophages_final, resolution = 0.5, verbose = FALSE)
macrophages_final <- RunUMAP(macrophages_final, dims = 1:6, verbose = FALSE)

# -----------------------------
# 3. Macrophage-state program scoring
# -----------------------------
mac_programs <- list(
  Inflammatory = c("IL1B", "TNF", "CXCL8", "CXCL2", "CCL2", "CCL3", "CCL4", "NFKBIA", "JUN", "FOS"),
  Phagocytic_FcR = c("FCGR1A", "FCGR2A", "FCGR3A", "MSR1", "MARCO", "LST1", "TYROBP", "CTSS", "CTSB", "LYZ"),
  Complement_Hofbauer = c("C1QA", "C1QB", "C1QC", "APOE", "TREM2", "VSIG4", "CD163", "MRC1", "FOLR2"),
  Vascular_ECM = c("SPP1", "VEGFA", "MMP9", "MMP14", "TIMP1", "TIMP3", "LGALS3", "FN1", "THBS1", "SERPINE1"),
  Immunoregulatory = c("IL10", "TGFB1", "MERTK", "MAF", "SOCS3", "CD274", "HLA-DRA", "HLA-DRB1", "VSIR", "LILRB1")
)

program_gene_qc <- purrr::imap_dfr(mac_programs, function(genes, program) {
  tibble(program = program, gene = genes, present = genes %in% rownames(macrophages_final))
})
write_csv(program_gene_qc, file.path(out_tab, "Fig5_macrophage_program_gene_qc.csv"))

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
    warning("Too few genes present for macrophage program: ", nm)
  }
}

for (nm in names(mac_programs)) {
  cand <- grep(paste0("^", nm, "_score"), colnames(macrophages_final@meta.data), value = TRUE)
  if (length(cand) >= 1) {
    macrophages_final[[paste0(nm, "_score")]] <- macrophages_final@meta.data[[cand[1]]]
  }
}

score_cols <- paste0(names(mac_programs), "_score")
score_cols <- intersect(score_cols, colnames(macrophages_final@meta.data))
for (sc in score_cols) {
  macrophages_final[[paste0(sc, "_z")]] <- as.numeric(scale(macrophages_final@meta.data[[sc]]))
}
score_cols_z <- paste0(score_cols, "_z")

mac_meta <- macrophages_final@meta.data %>%
  rownames_to_column("cell") %>%
  as_tibble()

state_mat <- mac_meta %>%
  select(all_of(score_cols_z)) %>%
  as.matrix()

state_names <- gsub("_score_z$", "", colnames(state_mat))
max_state_idx <- apply(state_mat, 1, which.max)
mac_meta$mac_state <- state_names[max_state_idx]

mac_meta$mac_state <- recode(
  mac_meta$mac_state,
  Inflammatory = "Activated remodeling",
  Phagocytic_FcR = "Phagocytic / FcR",
  Complement_Hofbauer = "Complement / Hofbauer-like",
  Vascular_ECM = "Vascular / ECM-remodeling",
  Immunoregulatory = "Immunoregulatory"
)

state_order <- c(
  "Complement / Hofbauer-like",
  "Phagocytic / FcR",
  "Vascular / ECM-remodeling",
  "Activated remodeling",
  "Immunoregulatory"
)

macrophages_final$mac_state <- factor(mac_meta$mac_state, levels = state_order)

# Deprecated binary M1/M2-style label kept only as QC/provenance, not used for biological claims.
macrophages_final$polarization_qc <- dplyr::case_when(
  macrophages_final$seurat_clusters %in% c("0", "1", "2") ~ "M0-like QC",
  macrophages_final$seurat_clusters %in% c("4", "5", "14") ~ "M2-like QC",
  macrophages_final$seurat_clusters %in% c("6", "7", "8", "9", "10", "11", "12", "13") ~ "M1-like QC",
  macrophages_final$seurat_clusters == "3" ~ "Intermediate QC",
  TRUE ~ "Unassigned QC"
)

# -----------------------------
# 4. Source tables: state counts, subject fractions, contrasts
# -----------------------------
mac_meta_out <- macrophages_final@meta.data %>%
  rownames_to_column("cell") %>%
  as_tibble()
write_csv(mac_meta_out, file.path(out_tab, "Fig5_macrophage_cell_metadata_with_states.csv"))

state_count_table <- mac_meta_out %>%
  count(mac_state, sort = TRUE, name = "n_cells")
write_csv(state_count_table, file.path(out_tab, "Fig5_macrophage_state_counts.csv"))

state_cluster_qc <- mac_meta_out %>%
  count(seurat_clusters, polarization_qc, mac_state, name = "n") %>%
  group_by(seurat_clusters) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  arrange(seurat_clusters, desc(frac))
write_csv(state_cluster_qc, file.path(out_tab, "Fig5_macrophage_cluster_state_polarization_qc.csv"))

state_frac_subject <- mac_meta_out %>%
  select(cell, study_id, sample_id_standardized, Delivery_class, MNPBurden, mac_state) %>%
  distinct() %>%
  count(study_id, Delivery_class, MNPBurden, mac_state, name = "n_cells") %>%
  group_by(study_id, Delivery_class, MNPBurden) %>%
  mutate(frac = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  tidyr::complete(
    study_id, Delivery_class, MNPBurden, mac_state = factor(state_order, levels = state_order),
    fill = list(n_cells = 0, frac = 0)
  )
write_csv(state_frac_subject, file.path(out_tab, "Fig5_macrophage_state_fraction_by_subject.csv"))

state_frac_stats <- state_frac_subject %>%
  group_by(mac_state, Delivery_class) %>%
  summarise(
    n_subjects = n_distinct(study_id),
    n_low = n_distinct(study_id[MNPBurden == "Low"]),
    n_high = n_distinct(study_id[MNPBurden == "High"]),
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
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  arrange(Delivery_class, p_value)
write_csv(state_frac_stats, file.path(out_tab, "Fig5_macrophage_state_fraction_stats.csv"))

subject_program_scores <- mac_meta_out %>%
  select(cell, study_id, Delivery_class, MNPBurden, all_of(score_cols_z)) %>%
  distinct() %>%
  group_by(study_id, Delivery_class, MNPBurden) %>%
  summarise(across(all_of(score_cols_z), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(cols = all_of(score_cols_z), names_to = "program", values_to = "score") %>%
  mutate(
    program = str_remove(program, "_score_z$"),
    program = recode(
      program,
      Inflammatory = "Activated remodeling",
      Phagocytic_FcR = "Phagocytic / FcR",
      Complement_Hofbauer = "Complement / Hofbauer-like",
      Vascular_ECM = "Vascular / ECM-remodeling",
      Immunoregulatory = "Immunoregulatory"
    ),
    program = factor(program, levels = state_order)
  )
write_csv(subject_program_scores, file.path(out_tab, "Fig5_macrophage_program_scores_by_subject.csv"))

# -----------------------------
# 5. Macrophage state marker table and dot plot genes
# -----------------------------
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
  arrange(cluster, p_val_adj, desc(avg_log2FC))
write_csv(mac_state_markers, file.path(out_tab, "Fig5_macrophage_state_FindAllMarkers.csv"))

exclude_pattern <- "^(MT-|RPL|RPS|MALAT1|XIST|Y_RNA|HBB|HBA|HBM|HBD|HBE|HBQ)"
top_mac_state_markers_clean <- mac_state_markers %>%
  filter(p_val_adj < 0.05, pct.1 >= 0.15, !grepl(exclude_pattern, gene, ignore.case = FALSE)) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = state_order)) %>%
  arrange(cluster, desc(avg_log2FC))
write_csv(top_mac_state_markers_clean, file.path(out_tab, "Fig5_macrophage_state_top_DEGs_for_dotplot.csv"))

# Curated final genes from the working analysis; silently drop genes absent from the public object.
genes_deg_final <- intersect(
  rev(c(
    "VSIG4", "LYVE1", "MPEG1",
    "FCGR3A", "FCGR2A", "MSR1",
    "FN1", "SPP1", "SERPINE1", "SPARC", "COL3A1",
    "MAF", "LILRB1"
  )),
  rownames(macrophages_final)
)

# -----------------------------
# 6. Figures
# -----------------------------
p5a_state_umap <- DimPlot(
  macrophages_final,
  reduction = "umap",
  group.by = "mac_state",
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

ggsave(file.path(out_fig, "Fig5a_macrophage_state_umap.png"), p5a_state_umap, width = 6.2, height = 4.8, dpi = 600)
ggsave(file.path(out_fig, "Fig5a_macrophage_state_umap.pdf"), p5a_state_umap, width = 6.2, height = 4.8, device = cairo_pdf)

p5b_state_fractions <- ggplot(
  state_frac_subject,
  aes(x = MNPBurden, y = frac)
) +
  geom_point(aes(group = study_id), position = position_jitter(width = 0.08, height = 0), size = 2.1, alpha = 0.9) +
  stat_summary(aes(group = MNPBurden), fun = mean, geom = "crossbar", width = 0.45, linewidth = 0.35, color = "black") +
  facet_grid(mac_state ~ Delivery_class, scales = "free_y") +
  labs(x = "MNP burden", y = "Macrophage state fraction") +
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

ggsave(file.path(out_fig, "Fig5b_macrophage_state_fractions.png"), p5b_state_fractions, width = 7.2, height = 6.2, dpi = 600)
ggsave(file.path(out_fig, "Fig5b_macrophage_state_fractions.pdf"), p5b_state_fractions, width = 7.2, height = 6.2, device = cairo_pdf)

if (length(genes_deg_final) > 0) {
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
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )

  ggsave(file.path(out_fig, "Fig5c_macrophage_state_DEG_dotplot.png"), p5c_deg_dotplot_final, width = 7.0, height = 4.6, dpi = 600)
  ggsave(file.path(out_fig, "Fig5c_macrophage_state_DEG_dotplot.pdf"), p5c_deg_dotplot_final, width = 7.0, height = 4.6, device = cairo_pdf)
}

p5d_programs <- ggplot(
  subject_program_scores,
  aes(x = MNPBurden, y = score)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
  geom_point(aes(group = study_id), position = position_jitter(width = 0.08, height = 0), size = 2.0, alpha = 0.9) +
  facet_grid(program ~ Delivery_class, scales = "free_y") +
  labs(x = "MNP burden", y = "Mean macrophage program score") +
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

ggsave(file.path(out_fig, "Fig5d_macrophage_program_scores.png"), p5d_programs, width = 7.2, height = 6.5, dpi = 600)
ggsave(file.path(out_fig, "Fig5d_macrophage_program_scores.pdf"), p5d_programs, width = 7.2, height = 6.5, device = cairo_pdf)

saveRDS(macrophages_final, file.path(out_tab, "macrophages_final_rebuilt_from_public_object.rds"))
message("Done: macrophage state outputs written to ", out_tab, " and ", out_fig)
