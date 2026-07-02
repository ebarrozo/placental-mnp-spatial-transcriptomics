# ============================================================
# Exploratory polymer-specific module associations
# Reproduces Supplementary Fig. 4c-d style source tables/plots
# from the public processed Seurat object and module-score outputs.
#
# Notes:
# - These analyses are hypothesis-generating because the nested spatial
#   cohort is small, platforms are heterogeneous, and polymer burdens are
#   sparse/collinear.
# - Primary module-score robustness remains in scripts/05_module_scores.R.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(broom)
})

if (!exists("obj")) {
  source("scripts/02_load_public_object.R")
}

module_outdir <- file.path("outputs", "tables", "module_scores")
module_file <- file.path(module_outdir, "subject_level_module_scores_platform_z.csv")
if (!file.exists(module_file)) {
  message("Module-score subject table not found; running scripts/05_module_scores.R first.")
  source("scripts/05_module_scores.R")
}

outdir <- file.path("outputs", "tables", "polymer_associations")
figdir <- file.path("outputs", "figures", "polymer_associations")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

subject_scores_z <- read_csv(module_file, show_col_types = FALSE) %>%
  mutate(
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

polymer_ug_cols <- c(
  PE = "pe_ug_per_g",
  SBR = "sbr_ug_per_g",
  PP = "pp_ug_per_g",
  PVC = "pvc_ug_per_g",
  N66 = "n66_ug_per_g",
  PET = "pet_ug_per_g",
  N6 = "n6_ug_per_g",
  PMMA = "pmma_ug_per_g",
  PC = "pc_ug_per_g",
  PU = "pu_ug_per_g",
  ABS = "abs_ug_per_g",
  PS = "ps_ug_per_g"
)

polymer_ug_cols <- polymer_ug_cols[polymer_ug_cols %in% colnames(obj@meta.data)]
if (length(polymer_ug_cols) == 0) {
  stop("No polymer ug/g columns found in obj@meta.data. Expected columns like pet_ug_per_g or pu_ug_per_g.")
}

module_z_cols <- grep("_score_platform_z$", colnames(subject_scores_z), value = TRUE)
if (length(module_z_cols) == 0) stop("No platform-z module score columns found in subject_scores_z.")

module_labels <- c(
  ECM_Remodeling_score_platform_z = "ECM remodeling",
  Vascular_Interface_score_platform_z = "Vascular interface",
  Macrophage_Remodeling_score_platform_z = "Macrophage remodeling",
  NK_Cytotoxic_Regulatory_score_platform_z = "NK cytotoxic/regulatory",
  Trophoblast_Remodeling_score_platform_z = "Trophoblast remodeling"
)
module_labels <- module_labels[names(module_labels) %in% module_z_cols]
module_order <- unname(module_labels)

polymer_subject <- obj@meta.data %>%
  as_tibble(rownames = "cell_barcode") %>%
  filter(.data$mnp_model_sample %in% TRUE) %>%
  group_by(study_id, Delivery_class, MNPBurden) %>%
  summarise(
    across(all_of(unname(polymer_ug_cols)), ~ {
      vals <- unique(.x[!is.na(.x)])
      if (length(vals) == 0) NA_real_ else vals[1]
    }),
    .groups = "drop"
  ) %>%
  mutate(
    across(all_of(unname(polymer_ug_cols)), ~ suppressWarnings(as.numeric(.x)))
  )

# Add log2(ug/g + 1) columns for modeling. This avoids depending on private
# intermediate metadata columns such as log2_mnp_pet_ug_g_plus1.
for (polymer in names(polymer_ug_cols)) {
  src <- unname(polymer_ug_cols[[polymer]])
  dst <- paste0("log2_mnp_", tolower(polymer), "_ug_g_plus1")
  polymer_subject[[dst]] <- log2(polymer_subject[[src]] + 1)
}

polymer_log2_cols <- grep("^log2_mnp_.*_ug_g_plus1$", colnames(polymer_subject), value = TRUE)

model_df <- subject_scores_z %>%
  left_join(
    polymer_subject %>% select(study_id, all_of(polymer_log2_cols), all_of(unname(polymer_ug_cols))),
    by = "study_id"
  )

polymer_missingness <- model_df %>%
  summarise(across(all_of(polymer_log2_cols), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "polymer_col", values_to = "n_missing") %>%
  mutate(
    n_subjects = nrow(model_df),
    n_nonmissing = n_subjects - n_missing,
    polymer = toupper(gsub("^log2_mnp_|_ug_g_plus1$", "", polymer_col))
  ) %>%
  arrange(polymer)

write_csv(polymer_missingness, file.path(outdir, "polymer_subject_level_missingness.csv"))
write_csv(model_df, file.path(outdir, "polymer_module_model_input_subject_level.csv"))

fit_polymer_module <- function(polymer_col, module_col, delivery, df = model_df) {
  dat <- df %>%
    filter(Delivery_class == delivery) %>%
    select(study_id, Delivery_class, MNPBurden, all_of(polymer_col), all_of(module_col)) %>%
    rename(polymer_value = all_of(polymer_col), module_score = all_of(module_col)) %>%
    drop_na(polymer_value, module_score)

  polymer <- toupper(gsub("^log2_mnp_|_ug_g_plus1$", "", polymer_col))
  module <- unname(module_labels[[module_col]] %||% module_col)

  if (nrow(dat) < 4 || length(unique(dat$polymer_value)) < 2) {
    return(tibble(
      polymer_col = polymer_col, polymer = polymer,
      module_col = module_col, module = module,
      Delivery_class = delivery,
      n = nrow(dat), slope = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_,
      model = "module_score_platform_z__polymer_log2ugpg_plus1"
    ))
  }

  fit <- lm(module_score ~ polymer_value, data = dat)
  tidy(fit) %>%
    filter(term == "polymer_value") %>%
    transmute(
      polymer_col = polymer_col,
      polymer = polymer,
      module_col = module_col,
      module = module,
      Delivery_class = delivery,
      n = nrow(dat),
      slope = estimate,
      std.error = std.error,
      statistic = statistic,
      p.value = p.value,
      model = "module_score_platform_z__polymer_log2ugpg_plus1"
    )
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a)) a else b

polymer_module_slopes <- tidyr::crossing(
  polymer_col = polymer_log2_cols,
  module_col = module_z_cols,
  Delivery_class = c("Term", "Preterm")
) %>%
  mutate(result = pmap(list(polymer_col, module_col, Delivery_class), fit_polymer_module)) %>%
  select(result) %>%
  unnest(result) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  arrange(Delivery_class, p.value)

write_csv(polymer_module_slopes, file.path(outdir, "polymer_module_within_delivery_slopes.csv"))

fit_polymer_loo <- function(polymer_col, module_col, delivery, left_out_id, df = model_df) {
  dat <- df %>% filter(study_id != left_out_id)
  fit_polymer_module(polymer_col, module_col, delivery, dat) %>%
    mutate(left_out_study_id = left_out_id)
}

loo_polymer_module <- tidyr::crossing(
  polymer_col = polymer_log2_cols,
  module_col = module_z_cols,
  Delivery_class = c("Term", "Preterm"),
  left_out_study_id = unique(model_df$study_id)
) %>%
  mutate(result = pmap(list(polymer_col, module_col, Delivery_class, left_out_study_id), fit_polymer_loo)) %>%
  select(result) %>%
  unnest(result)

loo_summary <- loo_polymer_module %>%
  group_by(polymer_col, polymer, module_col, module, Delivery_class) %>%
  summarise(
    n_iterations = n(),
    n_valid = sum(!is.na(slope)),
    median_slope = median(slope, na.rm = TRUE),
    min_slope = suppressWarnings(min(slope, na.rm = TRUE)),
    max_slope = suppressWarnings(max(slope, na.rm = TRUE)),
    prop_positive = mean(slope > 0, na.rm = TRUE),
    prop_negative = mean(slope < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    min_slope = ifelse(is.infinite(min_slope), NA_real_, min_slope),
    max_slope = ifelse(is.infinite(max_slope), NA_real_, max_slope)
  )

write_csv(loo_polymer_module, file.path(outdir, "polymer_module_leave_one_subject_out_iterations.csv"))
write_csv(loo_summary, file.path(outdir, "polymer_module_leave_one_subject_out_summary.csv"))

# -----------------------------
# Heatmap: within-delivery polymer-module slopes
# -----------------------------
polymer_order <- c("PE", "SBR", "PP", "PVC", "N66", "PET", "N6", "PMMA", "PC", "PU", "ABS", "PS")
polymer_order <- intersect(polymer_order, unique(polymer_module_slopes$polymer))
module_order <- unname(module_labels[module_z_cols])

heatmap_df <- polymer_module_slopes %>%
  mutate(
    polymer = factor(polymer, levels = polymer_order),
    module = factor(module, levels = module_order),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    sig_label = case_when(!is.na(p.value) & p.value < 0.05 ~ "●", TRUE ~ "")
  )

max_abs <- max(abs(heatmap_df$slope), na.rm = TRUE)
if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1

p_poly_heatmap <- ggplot(heatmap_df, aes(x = polymer, y = module, fill = slope)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3.2, color = "black") +
  facet_wrap(~ Delivery_class, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-max_abs, max_abs), name = "Slope"
  ) +
  labs(x = "Polymer burden", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 9),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11, color = "black"),
    legend.position = "right"
  )

ggsave(file.path(figdir, "SuppFig4c_polymer_module_within_delivery_slopes_heatmap.png"), p_poly_heatmap, width = 9.5, height = 3.8, units = "in", dpi = 600)
ggsave(file.path(figdir, "SuppFig4c_polymer_module_within_delivery_slopes_heatmap.pdf"), p_poly_heatmap, width = 9.5, height = 3.8, units = "in")

# -----------------------------
# Focused scatter/QC plot
# -----------------------------
focus_polymers <- intersect(
  c("log2_mnp_pet_ug_g_plus1", "log2_mnp_pu_ug_g_plus1", "log2_mnp_pvc_ug_g_plus1", "log2_mnp_sbr_ug_g_plus1", "log2_mnp_abs_ug_g_plus1"),
  polymer_log2_cols
)
focus_modules <- intersect(
  c("Macrophage_Remodeling_score_platform_z", "NK_Cytotoxic_Regulatory_score_platform_z", "Vascular_Interface_score_platform_z"),
  module_z_cols
)

if (length(focus_polymers) > 0 && length(focus_modules) > 0) {
  focus_long <- model_df %>%
    select(study_id, Delivery_class, MNPBurden, all_of(focus_polymers), all_of(focus_modules)) %>%
    pivot_longer(cols = all_of(focus_polymers), names_to = "polymer_col", values_to = "polymer_value") %>%
    pivot_longer(cols = all_of(focus_modules), names_to = "module_col", values_to = "module_score") %>%
    mutate(
      polymer = toupper(gsub("^log2_mnp_|_ug_g_plus1$", "", polymer_col)),
      polymer = factor(polymer, levels = c("PET", "PU", "PVC", "SBR", "ABS")),
      module = unname(module_labels[module_col]),
      module = factor(module, levels = c("Macrophage remodeling", "NK cytotoxic/regulatory", "Vascular interface")),
      Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm"))
    )

  p_poly_scatter <- ggplot(focus_long, aes(x = polymer_value, y = module_score, color = Delivery_class)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
    facet_grid(module ~ polymer, scales = "free_x") +
    scale_color_manual(values = c(Term = "#0F766E", Preterm = "#7C3AED")) +
    labs(x = "Polymer burden, log2(µg/g + 1)", y = "Platform-normalized module score", color = "Delivery class") +
    theme_classic(base_size = 10) +
    theme(
      text = element_text(family = "Arial"),
      axis.text.x = element_text(color = "black", size = 8),
      axis.text.y = element_text(color = "black", size = 8),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 8.5, color = "black"),
      legend.position = "top",
      panel.spacing.x = unit(0.65, "lines"),
      panel.spacing.y = unit(0.65, "lines")
    )

  ggsave(file.path(figdir, "SuppFig4d_polymer_module_focused_scatter_QC.png"), p_poly_scatter, width = 10.5, height = 5.2, units = "in", dpi = 600)
  ggsave(file.path(figdir, "SuppFig4d_polymer_module_focused_scatter_QC.pdf"), p_poly_scatter, width = 10.5, height = 5.2, units = "in")
}

message("Polymer-specific module outputs written to: ", outdir)
