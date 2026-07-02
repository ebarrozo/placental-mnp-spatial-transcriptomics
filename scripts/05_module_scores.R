# ============================================================
# Module score analysis from public Seurat object
# Reproduces subject-level platform-normalized module-score tables,
# within-delivery high-vs-low contrasts, LOO robustness, exact
# delivery-preserving permutation tests, and a Fig. 3-style module plot.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(emmeans)
  library(broom)
  library(ggbeeswarm)
})

if (!exists("obj")) {
  source("scripts/02_load_public_object.R")
}

outdir <- file.path("outputs", "tables", "module_scores")
figdir <- file.path("outputs", "figures", "module_scores")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

required_cols <- c(
  "sample_id_standardized", "study_id", "platform_clean",
  "Delivery_class", "MNPBurden", "MNPBurden_Delivery_class",
  "cumulative_mnp_ug_per_g", "cumulative_mnp_log10ugpg",
  "gestational_age_weeks_continuous", "mnp_model_sample"
)
missing_required <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_required) > 0) {
  stop("Missing required metadata columns: ", paste(missing_required, collapse = ", "))
}

# Module definitions used if scores are absent from the public object.
module_list <- list(
  ECM_Remodeling = c(
    "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2",
    "FN1", "SPARC", "TIMP3", "MMP14", "SERPINE1"
  ),
  Vascular_Interface = c(
    "PECAM1", "VWF", "KDR", "FLT1", "ENG", "ESAM", "NRP1", "PGF"
  ),
  Macrophage_Remodeling = c(
    "CD68", "CD163", "VSIG4", "FCGR1A", "FCGR3A",
    "C1QA", "C1QB", "C1QC", "MRC1", "MSR1"
  ),
  NK_Cytotoxic_Regulatory = c(
    "NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "KLRB1", "XCL1", "XCL2"
  ),
  Trophoblast_Remodeling = c(
    "KRT7", "GATA3", "HLA-G", "FLT1", "ENG", "TFPI", "LGALS13", "ERVFRD-1"
  )
)

module_score_cols <- paste0(names(module_list), "_score")

module_gene_check <- imap_dfr(module_list, function(genes, module) {
  tibble(module = module, gene = genes, present = genes %in% rownames(obj))
})
write_csv(module_gene_check, file.path(outdir, "module_gene_presence_qc.csv"))

# Prefer scores already stored in the public object. Recompute only if missing.
if (!all(module_score_cols %in% colnames(obj@meta.data))) {
  module_list_present <- map(module_list, ~ intersect(.x, rownames(obj)))
  module_list_present <- module_list_present[lengths(module_list_present) >= 2]
  if (length(module_list_present) == 0) stop("No module gene sets have >=2 genes present.")

  obj <- AddModuleScore(
    object = obj,
    features = module_list_present,
    name = names(module_list_present),
    assay = DefaultAssay(obj)
  )
  new_score_cols <- tail(colnames(obj@meta.data), length(module_list_present))
  clean_score_cols <- paste0(names(module_list_present), "_score")
  colnames(obj@meta.data)[match(new_score_cols, colnames(obj@meta.data))] <- clean_score_cols
  module_score_cols <- clean_score_cols
}

obj_mnp <- subset(obj, subset = mnp_model_sample == TRUE)

sample_scores <- obj_mnp@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  as_tibble() %>%
  group_by(
    sample_id_standardized,
    study_id,
    platform_clean,
    Delivery_class,
    MNPBurden,
    MNPBurden_Delivery_class,
    cumulative_mnp_ug_per_g,
    cumulative_mnp_log10ugpg,
    gestational_age_weeks_continuous
  ) %>%
  summarise(
    n_cells_or_spots = n(),
    across(all_of(module_score_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

sample_scores_z <- sample_scores %>%
  group_by(platform_clean) %>%
  mutate(across(all_of(module_score_cols), ~ as.numeric(scale(.x)), .names = "{.col}_platform_z")) %>%
  ungroup()

module_z_cols <- paste0(module_score_cols, "_platform_z")

subject_scores_z <- sample_scores_z %>%
  group_by(
    study_id,
    Delivery_class,
    MNPBurden,
    MNPBurden_Delivery_class,
    cumulative_mnp_ug_per_g,
    cumulative_mnp_log10ugpg,
    gestational_age_weeks_continuous
  ) %>%
  summarise(
    n_samples = n(),
    platforms = paste(sort(unique(as.character(platform_clean))), collapse = ";"),
    across(all_of(module_z_cols), ~ weighted.mean(.x, w = n_cells_or_spots, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

run_subject_model_z <- function(score_col, df = subject_scores_z) {
  fit <- lm(as.formula(paste0(score_col, " ~ MNPBurden * Delivery_class")), data = df)
  broom::tidy(fit, conf.int = TRUE) %>%
    mutate(score = score_col, model = "subject_level_platform_z_MNPBurden_x_Delivery")
}

subject_model_results_z <- map_dfr(module_z_cols, run_subject_model_z)

run_within_delivery_contrasts <- function(score_col, df = subject_scores_z) {
  fit <- lm(as.formula(paste0(score_col, " ~ MNPBurden * Delivery_class")), data = df)
  emm <- emmeans(fit, ~ MNPBurden | Delivery_class)
  contrast(emm, method = "revpairwise") %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(score = score_col, model = "within_delivery_high_vs_low_subject_platform_z")
}

within_delivery_contrasts <- map_dfr(module_z_cols, run_within_delivery_contrasts) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  arrange(score, Delivery_class)

run_loo <- function(left_out_id) {
  df_loo <- subject_scores_z %>% filter(study_id != left_out_id)
  map_dfr(module_z_cols, function(score_col) {
    run_within_delivery_contrasts(score_col, df_loo) %>%
      mutate(left_out_study_id = left_out_id)
  })
}

loo_contrasts <- map_dfr(unique(subject_scores_z$study_id), run_loo) %>%
  mutate(q.value = p.adjust(p.value, method = "BH"))

set.seed(1463)
permute_once <- function(i) {
  df_perm <- subject_scores_z %>%
    group_by(Delivery_class) %>%
    mutate(MNPBurden_perm = sample(MNPBurden)) %>%
    ungroup()
  map_dfr(module_z_cols, function(score_col) {
    fit <- lm(as.formula(paste0(score_col, " ~ MNPBurden_perm * Delivery_class")), data = df_perm)
    emm <- emmeans(fit, ~ MNPBurden_perm | Delivery_class)
    contrast(emm, method = "revpairwise") %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(score = score_col, permutation = i)
  })
}

n_perm <- as.integer(Sys.getenv("N_PERM", "10000"))
permutation_results <- map_dfr(seq_len(n_perm), permute_once)

observed_for_perm <- within_delivery_contrasts %>%
  select(score, Delivery_class, observed_estimate = estimate)

permutation_summary <- permutation_results %>%
  left_join(observed_for_perm, by = c("score", "Delivery_class")) %>%
  group_by(score, Delivery_class) %>%
  summarise(
    observed_estimate = first(observed_estimate),
    n_perm = n(),
    two_sided_p = mean(abs(estimate) >= abs(observed_estimate), na.rm = TRUE),
    one_sided_positive_p = mean(estimate >= observed_estimate, na.rm = TRUE),
    one_sided_negative_p = mean(estimate <= observed_estimate, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(sample_scores, file.path(outdir, "sample_level_module_scores.csv"))
write_csv(sample_scores_z, file.path(outdir, "sample_level_module_scores_platform_z.csv"))
write_csv(subject_scores_z, file.path(outdir, "subject_level_module_scores_platform_z.csv"))
write_csv(subject_model_results_z, file.path(outdir, "subject_level_platform_z_mnpburden_interaction_models.csv"))
write_csv(within_delivery_contrasts, file.path(outdir, "within_delivery_high_vs_low_platform_z_contrasts.csv"))
write_csv(loo_contrasts, file.path(outdir, "loo_primary_mnpburden_within_delivery_contrasts.csv"))
write_csv(permutation_summary, file.path(outdir, "delivery_preserving_permutation_summary.csv"))

module_labels <- c(
  ECM_Remodeling_score_platform_z = "ECM\nremodeling",
  Vascular_Interface_score_platform_z = "Vascular\ninterface",
  Macrophage_Remodeling_score_platform_z = "Macrophage\nremodeling",
  NK_Cytotoxic_Regulatory_score_platform_z = "NK cytotoxic /\nregulatory",
  Trophoblast_Remodeling_score_platform_z = "Trophoblast\nremodeling"
)
module_order <- names(module_labels)

plot_df <- subject_scores_z %>%
  select(study_id, Delivery_class, MNPBurden, n_samples, platforms, all_of(module_z_cols)) %>%
  pivot_longer(all_of(module_z_cols), names_to = "module", values_to = "module_score_platform_z") %>%
  mutate(
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High")),
    module = factor(module, levels = module_order, labels = module_labels[module_order])
  )

p_annot <- within_delivery_contrasts %>%
  filter(score %in% module_z_cols) %>%
  mutate(
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    module = factor(score, levels = module_order, labels = module_labels[module_order]),
    p_label = paste0("p = ", signif(p.value, 2))
  ) %>%
  left_join(
    plot_df %>% group_by(Delivery_class, module) %>% summarise(y.position = max(module_score_platform_z, na.rm = TRUE) + 0.35, .groups = "drop"),
    by = c("Delivery_class", "module")
  ) %>%
  filter(p.value < 0.10)

p_module <- ggplot(plot_df, aes(x = module, y = module_score_platform_z, color = MNPBurden, fill = MNPBurden)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey75") +
  geom_boxplot(width = 0.52, alpha = 0.18, outlier.shape = NA, linewidth = 0.45, position = position_dodge(width = 0.68)) +
  ggbeeswarm::geom_quasirandom(aes(shape = MNPBurden), dodge.width = 0.68, width = 0.08, size = 2.8, alpha = 0.92, stroke = 0.45) +
  geom_text(data = p_annot, aes(x = module, y = y.position, label = p_label), inherit.aes = FALSE, size = 3.2, color = "black") +
  facet_wrap(~ Delivery_class, nrow = 1) +
  scale_color_manual(values = c(Low = "#64748B", High = "#D97706")) +
  scale_fill_manual(values = c(Low = "#64748B", High = "#D97706")) +
  scale_shape_manual(values = c(Low = 21, High = 24)) +
  labs(x = NULL, y = "Platform-normalized module score", color = "MNP burden", fill = "MNP burden", shape = "MNP burden") +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(color = "black", size = 9.5, lineheight = 0.9),
    axis.text.y = element_text(color = "black", size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    legend.position = "top",
    panel.spacing.x = unit(1.0, "lines")
  )

ggsave(file.path(figdir, "Fig3_module_scores_platform_z_subject_level.png"), p_module, width = 9, height = 3.8, units = "in", dpi = 600)
ggsave(file.path(figdir, "Fig3_module_scores_platform_z_subject_level.pdf"), p_module, width = 9, height = 3.8, units = "in")

message("Module-score outputs written to: ", outdir)
