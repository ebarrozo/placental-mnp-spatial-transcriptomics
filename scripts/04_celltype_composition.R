# ============================================================
# Fig. 2 composition summaries from public Seurat object
# Summarizes predicted cell-type composition across platforms/subjects,
# estimates high-vs-low MNP contrasts within delivery class, and exports
# source tables plus draft figure panels.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggbeeswarm)
  library(emmeans)
  library(purrr)
})

if (!exists("obj")) {
  source("scripts/02_load_public_object.R")
}

outdir <- file.path("outputs", "tables", "composition")
figdir <- file.path("outputs", "figures", "composition")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

annotation_col <- if ("nn.predicted.celltype" %in% colnames(obj@meta.data)) "nn.predicted.celltype" else "celltype_display"
required_cols <- c("sample_id_standardized", "study_id", "platform_clean", "Delivery_class", "MNPBurden", "mnp_model_sample", annotation_col)
missing_required <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_required) > 0) stop("Missing required columns: ", paste(missing_required, collapse = ", "))

meta_all <- obj@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  as_tibble() %>%
  mutate(
    celltype_main = as.character(.data[[annotation_col]]),
    celltype_main = if_else(is.na(celltype_main) | celltype_main == "", "Unassigned", celltype_main),
    platform_clean = factor(platform_clean, levels = c("VisiumV1", "Xenium", "VisiumHD")),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

annotation_counts_all <- meta_all %>%
  count(platform_clean, celltype_main, name = "n_cells_or_spots") %>%
  group_by(platform_clean) %>%
  mutate(prop = n_cells_or_spots / sum(n_cells_or_spots)) %>%
  ungroup() %>%
  arrange(platform_clean, desc(prop))
write_csv(annotation_counts_all, file.path(outdir, "annotation_counts_all_by_platform.csv"))

annotation_counts_mnp <- meta_all %>%
  filter(mnp_model_sample == TRUE) %>%
  count(platform_clean, Delivery_class, MNPBurden, celltype_main, name = "n_cells_or_spots") %>%
  group_by(platform_clean, Delivery_class, MNPBurden) %>%
  mutate(prop = n_cells_or_spots / sum(n_cells_or_spots)) %>%
  ungroup()
write_csv(annotation_counts_mnp, file.path(outdir, "annotation_counts_mnp_by_platform_delivery_burden.csv"))

sample_comp_long <- meta_all %>%
  group_by(sample_id_standardized, study_id, platform_clean, Delivery_class, MNPBurden, mnp_model_sample, celltype_main) %>%
  summarise(n_celltype = n(), .groups = "drop") %>%
  group_by(sample_id_standardized) %>%
  mutate(n_total = sum(n_celltype), prop = n_celltype / n_total, pct = prop * 100) %>%
  ungroup()

all_sample_keys <- meta_all %>% distinct(sample_id_standardized, study_id, platform_clean, Delivery_class, MNPBurden, mnp_model_sample)
all_celltypes <- sample_comp_long %>% distinct(celltype_main)

sample_comp_long <- all_sample_keys %>%
  tidyr::crossing(all_celltypes) %>%
  left_join(sample_comp_long, by = c("sample_id_standardized", "study_id", "platform_clean", "Delivery_class", "MNPBurden", "mnp_model_sample", "celltype_main")) %>%
  group_by(sample_id_standardized) %>%
  mutate(
    n_celltype = replace_na(n_celltype, 0L),
    n_total = if_else(is.na(n_total), sum(n_celltype), n_total),
    prop = replace_na(prop, 0),
    pct = replace_na(pct, 0)
  ) %>%
  ungroup()

subject_comp_long_weighted <- sample_comp_long %>%
  filter(mnp_model_sample == TRUE) %>%
  group_by(study_id, Delivery_class, MNPBurden, celltype_main) %>%
  summarise(
    n_celltype = sum(n_celltype),
    n_total = sum(n_total),
    prop = n_celltype / n_total,
    pct = prop * 100,
    n_samples = n_distinct(sample_id_standardized),
    platforms = paste(sort(unique(as.character(platform_clean))), collapse = ";"),
    .groups = "drop"
  )

subject_comp_long_unweighted <- sample_comp_long %>%
  filter(mnp_model_sample == TRUE) %>%
  group_by(study_id, Delivery_class, MNPBurden, celltype_main) %>%
  summarise(
    prop = mean(prop, na.rm = TRUE),
    pct = prop * 100,
    n_samples = n_distinct(sample_id_standardized),
    platforms = paste(sort(unique(as.character(platform_clean))), collapse = ";"),
    .groups = "drop"
  )

subject_comp_wide <- subject_comp_long_weighted %>%
  select(study_id, Delivery_class, MNPBurden, celltype_main, prop) %>%
  mutate(celltype_main = make.names(celltype_main)) %>%
  pivot_wider(names_from = celltype_main, values_from = prop, names_prefix = "prop_", values_fill = 0)

write_csv(sample_comp_long, file.path(outdir, "sample_level_celltype_composition_long.csv"))
write_csv(subject_comp_long_weighted, file.path(outdir, "subject_level_celltype_composition_weighted_long.csv"))
write_csv(subject_comp_long_unweighted, file.path(outdir, "subject_level_celltype_composition_unweighted_long.csv"))
write_csv(subject_comp_wide, file.path(outdir, "subject_level_celltype_composition.csv"))

top_celltypes <- subject_comp_long_weighted %>%
  group_by(celltype_main) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_prop)) %>%
  slice_head(n = 12) %>%
  pull(celltype_main)

subject_comp_model <- subject_comp_long_weighted %>%
  mutate(
    prop_offset = (n_celltype + 0.5) / (n_total + 1),
    logit_prop = qlogis(prop_offset),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

run_celltype_contrast <- function(ct, df = subject_comp_model) {
  dat <- df %>% filter(celltype_main == ct)
  if (n_distinct(dat$MNPBurden) < 2 || n_distinct(dat$Delivery_class) < 2) return(tibble())
  fit <- lm(logit_prop ~ MNPBurden * Delivery_class, data = dat)
  emm <- emmeans(fit, ~ MNPBurden | Delivery_class)
  contrast(emm, method = list("High - Low" = c(-1, 1))) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(celltype_main = ct, model = "subject_level_logit_prop_high_vs_low_within_delivery")
}

celltype_composition_contrasts <- map_dfr(sort(unique(subject_comp_model$celltype_main)), run_celltype_contrast) %>%
  group_by(Delivery_class) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(Delivery_class, p.value)
write_csv(celltype_composition_contrasts, file.path(outdir, "within_delivery_celltype_composition_contrasts.csv"))

# Draft composition by sample/platform panel.
plot_sample_comp <- sample_comp_long %>%
  filter(celltype_main %in% top_celltypes) %>%
  mutate(sample_label = paste0(sample_id_standardized, "\n", platform_clean), sample_label = fct_reorder(sample_label, as.numeric(platform_clean)))

p_platform <- ggplot(plot_sample_comp, aes(x = sample_label, y = pct, fill = celltype_main)) +
  geom_col(width = 0.85, color = NA) +
  facet_grid(mnp_model_sample ~ Delivery_class, scales = "free_x", space = "free_x", labeller = labeller(mnp_model_sample = c(`TRUE` = "Matched MNP subset", `FALSE` = "Spatial annotation / validation"))) +
  labs(x = NULL, y = "Predicted cell-type composition (%)", fill = "Predicted cell type") +
  theme_classic(base_size = 10) +
  theme(text = element_text(family = "Arial"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10), legend.position = "right")

ggsave(file.path(figdir, "Fig2_celltype_composition_by_platform.png"), p_platform, width = 11, height = 5.5, units = "in", dpi = 600)

composition_celltypes_focus <- c("Placental Endothelial", "Macrophage", "NK Cell", "ILC3", "Plasma Cell", "Villous Cytotrophoblast", "Extravillous Trophoblast")
plot_focus_comp <- subject_comp_long_weighted %>%
  filter(celltype_main %in% composition_celltypes_focus) %>%
  mutate(celltype_main = factor(celltype_main, levels = composition_celltypes_focus), Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")), MNPBurden = factor(MNPBurden, levels = c("Low", "High")))

p_focus <- ggplot(plot_focus_comp, aes(x = MNPBurden, y = pct, color = MNPBurden, fill = MNPBurden)) +
  geom_boxplot(width = 0.48, alpha = 0.18, outlier.shape = NA, linewidth = 0.4) +
  ggbeeswarm::geom_quasirandom(width = 0.08, size = 2.8, alpha = 0.95, stroke = 0.4) +
  facet_grid(Delivery_class ~ celltype_main, scales = "free_y") +
  scale_color_manual(values = c(Low = "#64748B", High = "#D97706")) +
  scale_fill_manual(values = c(Low = "#64748B", High = "#D97706")) +
  labs(x = NULL, y = "Predicted composition (%)", color = "MNP burden", fill = "MNP burden") +
  theme_classic(base_size = 10) +
  theme(text = element_text(family = "Arial"), strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9), legend.position = "top")

ggsave(file.path(figdir, "Fig2_targeted_celltype_composition.png"), p_focus, width = 10.5, height = 4.8, units = "in", dpi = 600)

message("Composition outputs written to: ", outdir)
