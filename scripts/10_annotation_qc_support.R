# 10_annotation_qc_support.R
# Purpose: Recreate reviewer-facing annotation QC/support outputs from the public Seurat object.
# Covers: Supplementary Fig. 1b-d, Supplementary Fig. 2a-b, and SF3-style annotation/support summaries
# where the required metadata/reductions are present in obj_public.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

if (!exists("obj_public")) {
  if (file.exists("scripts/02_load_public_object.R")) {
    source("scripts/02_load_public_object.R")
  } else {
    stop("obj_public not found and scripts/02_load_public_object.R is unavailable.")
  }
}

obj <- obj_public

out_tables <- file.path("outputs", "tables", "annotation_qc")
out_figs <- file.path("outputs", "figures", "annotation_qc")
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figs, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helpers
# -----------------------------
has_reduction <- function(object, reduction = "umap") {
  reduction %in% names(object@reductions)
}

save_plot <- function(p, filename, width = 8, height = 5, dpi = 600) {
  ggsave(filename = filename, plot = p, width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)
}

make_dimplot <- function(object, group.by, split.by = NULL, filename_stub, title = NULL,
                         width = 10, height = 5, raster = TRUE) {
  if (!has_reduction(object, "umap")) {
    warning("UMAP reduction is not present; skipping ", filename_stub)
    return(invisible(NULL))
  }
  if (!group.by %in% colnames(object@meta.data)) {
    warning("Metadata column not present: ", group.by, "; skipping ", filename_stub)
    return(invisible(NULL))
  }
  if (!is.null(split.by) && !split.by %in% colnames(object@meta.data)) {
    warning("Split metadata column not present: ", split.by, "; skipping ", filename_stub)
    return(invisible(NULL))
  }

  p <- DimPlot(
    object,
    reduction = "umap",
    group.by = group.by,
    split.by = split.by,
    label = FALSE,
    raster = raster
  ) +
    NoAxes() +
    theme_classic(base_size = 10) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank()
    )

  if (!is.null(title)) p <- p + ggtitle(title)

  save_plot(p, file.path(out_figs, paste0(filename_stub, ".png")), width = width, height = height)
  save_plot(p, file.path(out_figs, paste0(filename_stub, ".pdf")), width = width, height = height, dpi = 300)
  invisible(p)
}

# -----------------------------
# Core metadata table
# -----------------------------
meta <- obj@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  as_tibble()

write_csv(
  tibble(
    metadata_column = colnames(meta),
    n_missing = vapply(meta, function(x) sum(is.na(x)), integer(1)),
    n_unique = vapply(meta, function(x) length(unique(x)), integer(1))
  ),
  file.path(out_tables, "annotation_qc_metadata_column_summary.csv")
)

# -----------------------------
# Supplementary Fig. 1b: UMAP split by platform
# -----------------------------
platform_col <- dplyr::case_when(
  "platform_clean" %in% colnames(meta) ~ "platform_clean",
  "platform" %in% colnames(meta) ~ "platform",
  TRUE ~ NA_character_
)
annotation_col <- dplyr::case_when(
  "celltype_display" %in% colnames(meta) ~ "celltype_display",
  "nn.predicted.celltype" %in% colnames(meta) ~ "nn.predicted.celltype",
  "predicted.celltype" %in% colnames(meta) ~ "predicted.celltype",
  TRUE ~ NA_character_
)

if (!is.na(annotation_col) && !is.na(platform_col)) {
  make_dimplot(
    obj,
    group.by = annotation_col,
    split.by = platform_col,
    filename_stub = "SuppFig1b_umap_annotation_split_by_platform",
    title = "Integrated atlas annotation by platform",
    width = 14,
    height = 5
  )
}

# Platform/cell-type support tables used for SF1c/SF1d-style summaries.
if (!is.na(annotation_col) && !is.na(platform_col)) {
  annotation_by_platform <- meta %>%
    mutate(
      annotation = as.character(.data[[annotation_col]]),
      platform = as.character(.data[[platform_col]])
    ) %>%
    filter(!is.na(annotation), !is.na(platform)) %>%
    count(platform, annotation, name = "n_cells_or_spots") %>%
    group_by(platform) %>%
    mutate(fraction = n_cells_or_spots / sum(n_cells_or_spots)) %>%
    ungroup() %>%
    arrange(platform, desc(fraction))

  write_csv(annotation_by_platform, file.path(out_tables, "SuppFig1c_annotation_counts_by_platform.csv"))

  p_platform_comp <- annotation_by_platform %>%
    group_by(annotation) %>%
    mutate(mean_fraction = mean(fraction, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(mean_fraction > 0.005) %>%
    ggplot(aes(x = platform, y = fraction, fill = annotation)) +
    geom_col(position = "fill", width = 0.85) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = NULL, y = "Predicted composition", fill = "Annotation") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

  save_plot(p_platform_comp, file.path(out_figs, "SuppFig1c_annotation_composition_by_platform.png"), width = 9, height = 5)
  save_plot(p_platform_comp, file.path(out_figs, "SuppFig1c_annotation_composition_by_platform.pdf"), width = 9, height = 5, dpi = 300)
}

# Annotation support / forced-label QC, when columns are available.
support_cols <- intersect(c("predicted.celltype_support", "nn.predicted.celltype_support"), colnames(meta))
if (!is.na(annotation_col) && length(support_cols) > 0) {
  support_col <- support_cols[1]

  support_summary <- meta %>%
    mutate(
      annotation = as.character(.data[[annotation_col]]),
      support = suppressWarnings(as.numeric(.data[[support_col]]))
    ) %>%
    group_by(annotation) %>%
    summarise(
      n_cells_or_spots = n(),
      median_support = median(support, na.rm = TRUE),
      mean_support = mean(support, na.rm = TRUE),
      pct_support_ge_0_40 = mean(support >= 0.40, na.rm = TRUE),
      pct_support_ge_0_60 = mean(support >= 0.60, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_cells_or_spots))

  write_csv(support_summary, file.path(out_tables, "SuppFig1d_forced_annotation_support_by_celltype.csv"))

  plot_support <- meta %>%
    mutate(
      annotation = as.character(.data[[annotation_col]]),
      support = suppressWarnings(as.numeric(.data[[support_col]]))
    ) %>%
    filter(!is.na(annotation), !is.na(support)) %>%
    semi_join(support_summary %>% slice_max(n_cells_or_spots, n = 16), by = "annotation") %>%
    ggplot(aes(x = reorder(annotation, support, FUN = median), y = support)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
    coord_flip() +
    labs(x = NULL, y = "KNN annotation support") +
    theme_classic(base_size = 10)

  save_plot(plot_support, file.path(out_figs, "SuppFig1d_forced_annotation_support.png"), width = 7, height = 5)
  save_plot(plot_support, file.path(out_figs, "SuppFig1d_forced_annotation_support.pdf"), width = 7, height = 5, dpi = 300)
}

# -----------------------------
# Supplementary Fig. 2a-b: maternal/fetal origin support
# -----------------------------
origin_col <- if ("origin_rule_based" %in% colnames(meta)) "origin_rule_based" else NA_character_
if (!is.na(origin_col)) {
  make_dimplot(
    obj,
    group.by = origin_col,
    split.by = origin_col,
    filename_stub = "SuppFig2a_umap_split_by_origin_rule_based",
    title = "Rule-based maternal/fetal origin assignment",
    width = 12,
    height = 4.5
  )

  origin_by_celltype <- meta %>%
    mutate(
      origin = as.character(.data[[origin_col]]),
      annotation = if (!is.na(annotation_col)) as.character(.data[[annotation_col]]) else NA_character_
    ) %>%
    filter(!is.na(origin), !is.na(annotation)) %>%
    count(annotation, origin, name = "n_cells_or_spots") %>%
    group_by(annotation) %>%
    mutate(fraction = n_cells_or_spots / sum(n_cells_or_spots)) %>%
    ungroup() %>%
    arrange(annotation, desc(fraction))

  write_csv(origin_by_celltype, file.path(out_tables, "SuppFig2b_origin_by_annotation.csv"))

  p_origin <- origin_by_celltype %>%
    ggplot(aes(x = annotation, y = fraction, fill = origin)) +
    geom_col(width = 0.85) +
    scale_y_continuous(labels = scales::percent_format()) +
    coord_flip() +
    labs(x = NULL, y = "Origin assignment fraction", fill = "Origin") +
    theme_classic(base_size = 10)

  save_plot(p_origin, file.path(out_figs, "SuppFig2b_origin_by_annotation.png"), width = 8, height = 5.5)
  save_plot(p_origin, file.path(out_figs, "SuppFig2b_origin_by_annotation.pdf"), width = 8, height = 5.5, dpi = 300)
}

# Sex-linked support tables, if RNA-assay-derived scores were retained.
sex_support_cols <- intersect(c("XIST_score", "Y_score", "origin_sex_supported", "origin_confidence", "origin_reason"), colnames(meta))
if (length(sex_support_cols) > 0) {
  sex_support_summary <- meta %>%
    select(any_of(c("cell_barcode", annotation_col, origin_col, sex_support_cols)))
  write_csv(sex_support_summary, file.path(out_tables, "maternal_fetal_origin_sex_linked_support_cell_level.csv"))
}

# -----------------------------
# SF3-style annotation/sample/platform support summaries
# -----------------------------
if (!is.na(annotation_col)) {
  sample_col <- if ("sample_id_standardized" %in% colnames(meta)) "sample_id_standardized" else if ("sample_id_std" %in% colnames(meta)) "sample_id_std" else NA_character_
  if (!is.na(sample_col)) {
    annot_by_sample <- meta %>%
      mutate(
        sample_id = as.character(.data[[sample_col]]),
        annotation = as.character(.data[[annotation_col]]),
        platform = if (!is.na(platform_col)) as.character(.data[[platform_col]]) else NA_character_
      ) %>%
      filter(!is.na(sample_id), !is.na(annotation)) %>%
      count(sample_id, platform, annotation, name = "n_cells_or_spots") %>%
      group_by(sample_id) %>%
      mutate(fraction = n_cells_or_spots / sum(n_cells_or_spots)) %>%
      ungroup()

    write_csv(annot_by_sample, file.path(out_tables, "SuppFig3a_annotation_counts_by_sample.csv"))

    p_sample_comp <- annot_by_sample %>%
      group_by(annotation) %>%
      mutate(mean_fraction = mean(fraction, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(mean_fraction > 0.005) %>%
      ggplot(aes(x = sample_id, y = fraction, fill = annotation)) +
      geom_col(position = "fill", width = 0.85) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(x = NULL, y = "Predicted composition", fill = "Annotation") +
      theme_classic(base_size = 9) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "right")

    save_plot(p_sample_comp, file.path(out_figs, "SuppFig3a_annotation_composition_by_sample.png"), width = 12, height = 5.5)
    save_plot(p_sample_comp, file.path(out_figs, "SuppFig3a_annotation_composition_by_sample.pdf"), width = 12, height = 5.5, dpi = 300)
  }
}

message("Annotation QC outputs written to: ", out_tables, " and ", out_figs)
