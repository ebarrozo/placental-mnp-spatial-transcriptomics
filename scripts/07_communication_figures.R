# ============================================================
# Communication analyses: NicheNet / CellChat source-table workflow
# Purpose:
#   Reproduce manuscript communication figure/source outputs from exported
#   NicheNet and CellChat source tables when available.
#
# Scope:
#   The original high-memory CellChat/NicheNet exploratory scripts are archived
#   in provenance/archived_scripts/. This public-facing script is intentionally
#   source-table-first because full re-execution requires large objects,
#   external interaction priors, package-version-specific APIs, and substantial
#   memory/runtime.
#
# Expected primary inputs, if available:
#   data/derived/communication/nichenet_ligand_network_edges.csv
#   data/derived/communication/nichenet_top_ligands.csv
#   data/derived/communication/cellchat_core_pathway_sender_contributions.csv
#
# Outputs:
#   outputs/tables/communication/*
#   outputs/figures/communication/*
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(forcats)
})

source("R/path_helpers.R")

communication_input_dir <- file.path("data", "derived", "communication")
communication_table_dir <- safe_dir(file.path("outputs", "tables", "communication"))
communication_figure_dir <- safe_dir(file.path("outputs", "figures", "communication"))

expected_sources <- tibble::tribble(
  ~file_name, ~purpose, ~required_for,
  "nichenet_ligand_network_edges.csv", "Long-form NicheNet sender-ligand-receptor-target or sender-ligand-target edges", "Fig4a/Fig4b network summaries",
  "nichenet_top_ligands.csv", "Top NicheNet ligand activity table by receiver, delivery class, and MNP contrast", "Fig4 ligand prioritization/source table",
  "cellchat_core_pathway_sender_contributions.csv", "CellChat pathway-level sender contributions for selected ECM/adhesion/chemokine/cytokine pathways", "Supplementary Fig5 alluvial/source table"
) %>%
  mutate(
    path = file.path(communication_input_dir, file_name),
    present = file.exists(path)
  )

readr::write_csv(
  expected_sources,
  file.path(communication_table_dir, "communication_expected_source_tables.csv")
)

message("Communication source-table availability:")
print(expected_sources %>% select(file_name, present, required_for), n = Inf)

# ------------------------------------------------------------
# Helper plotting functions
# ------------------------------------------------------------

standardize_name <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

plot_nichenet_top_ligands <- function(top_ligands, out_prefix = "Fig4_nichenet_top_ligands") {
  needed <- c("receiver", "ligand")
  if (!all(needed %in% colnames(top_ligands))) {
    warning("nichenet_top_ligands.csv must contain at least receiver and ligand columns.")
    return(invisible(NULL))
  }

  score_col <- dplyr::case_when(
    "aupr_corrected" %in% colnames(top_ligands) ~ "aupr_corrected",
    "activity_score" %in% colnames(top_ligands) ~ "activity_score",
    "score" %in% colnames(top_ligands) ~ "score",
    TRUE ~ NA_character_
  )

  if (is.na(score_col)) {
    top_ligands <- top_ligands %>%
      group_by(receiver) %>%
      mutate(score_for_plot = rev(seq_len(n()))) %>%
      ungroup()
    score_col <- "score_for_plot"
  }

  plot_df <- top_ligands %>%
    mutate(
      receiver = standardize_name(receiver),
      ligand = as.character(ligand)
    ) %>%
    group_by(receiver) %>%
    slice_max(order_by = .data[[score_col]], n = 15, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(ligand = fct_reorder(ligand, .data[[score_col]]))

  p <- ggplot(plot_df, aes(x = .data[[score_col]], y = ligand)) +
    geom_col(width = 0.75) +
    facet_wrap(~ receiver, scales = "free_y") +
    labs(x = score_col, y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      text = element_text(family = "Arial"),
      axis.text = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black")
    )

  ggsave(file.path(communication_figure_dir, paste0(out_prefix, ".png")), p, width = 8, height = 5, dpi = 600)
  ggsave(file.path(communication_figure_dir, paste0(out_prefix, ".pdf")), p, width = 8, height = 5)
  invisible(p)
}

plot_nichenet_edge_summary <- function(edges, out_prefix = "Fig4_nichenet_edge_summary") {
  possible_sender <- intersect(c("sender", "sender_celltype", "celltype", "source"), colnames(edges))[1]
  possible_ligand <- intersect(c("ligand", "from", "test_ligand"), colnames(edges))[1]
  possible_receiver <- intersect(c("receiver", "receiver_celltype", "target_celltype"), colnames(edges))[1]
  possible_weight <- intersect(c("weight", "regulatory_potential", "score", "aupr_corrected", "delta"), colnames(edges))[1]

  if (any(is.na(c(possible_sender, possible_ligand)))) {
    warning("nichenet_ligand_network_edges.csv needs sender/source and ligand columns for summary plotting.")
    return(invisible(NULL))
  }

  plot_df <- edges %>%
    mutate(
      sender_for_plot = standardize_name(.data[[possible_sender]]),
      ligand_for_plot = as.character(.data[[possible_ligand]]),
      receiver_for_plot = if (!is.na(possible_receiver)) standardize_name(.data[[possible_receiver]]) else "All receivers",
      weight_for_plot = if (!is.na(possible_weight)) abs(as.numeric(.data[[possible_weight]])) else 1
    ) %>%
    filter(!is.na(sender_for_plot), !is.na(ligand_for_plot)) %>%
    group_by(receiver_for_plot, sender_for_plot, ligand_for_plot) %>%
    summarise(edge_weight = sum(weight_for_plot, na.rm = TRUE), n_edges = n(), .groups = "drop") %>%
    group_by(receiver_for_plot) %>%
    slice_max(order_by = edge_weight, n = 30, with_ties = FALSE) %>%
    ungroup()

  p <- ggplot(plot_df, aes(x = sender_for_plot, y = ligand_for_plot, size = edge_weight)) +
    geom_point(alpha = 0.8) +
    facet_wrap(~ receiver_for_plot, scales = "free") +
    labs(x = "Sender cell type", y = "Prioritized ligand", size = "Edge weight") +
    theme_classic(base_size = 10) +
    theme(
      text = element_text(family = "Arial"),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black")
    )

  ggsave(file.path(communication_figure_dir, paste0(out_prefix, ".png")), p, width = 9, height = 5.5, dpi = 600)
  ggsave(file.path(communication_figure_dir, paste0(out_prefix, ".pdf")), p, width = 9, height = 5.5)
  invisible(p)
}

plot_cellchat_sender_contributions <- function(cellchat_tbl, out_prefix = "SuppFig5_cellchat_sender_contributions") {
  needed <- c("pathway", "sender")
  if (!all(needed %in% colnames(cellchat_tbl))) {
    warning("cellchat_core_pathway_sender_contributions.csv must contain at least pathway and sender columns.")
    return(invisible(NULL))
  }

  value_col <- dplyr::case_when(
    "interaction_probability" %in% colnames(cellchat_tbl) ~ "interaction_probability",
    "prob" %in% colnames(cellchat_tbl) ~ "prob",
    "weight" %in% colnames(cellchat_tbl) ~ "weight",
    "contribution" %in% colnames(cellchat_tbl) ~ "contribution",
    TRUE ~ NA_character_
  )

  if (is.na(value_col)) {
    cellchat_tbl <- cellchat_tbl %>% mutate(contribution_for_plot = 1)
    value_col <- "contribution_for_plot"
  }

  grouping_cols <- intersect(c("Delivery_class", "MNPBurden", "condition", "comparison"), colnames(cellchat_tbl))

  plot_df <- cellchat_tbl %>%
    mutate(
      pathway = standardize_name(pathway),
      sender = standardize_name(sender),
      value_for_plot = as.numeric(.data[[value_col]])
    ) %>%
    group_by(across(all_of(c(grouping_cols, "pathway", "sender")))) %>%
    summarise(value_for_plot = sum(value_for_plot, na.rm = TRUE), .groups = "drop")

  if (length(grouping_cols) == 0) {
    p <- ggplot(plot_df, aes(x = pathway, y = value_for_plot, fill = sender)) +
      geom_col(width = 0.8) +
      labs(x = NULL, y = "Summed communication contribution", fill = "Sender")
  } else {
    plot_df <- plot_df %>% unite("facet_group", all_of(grouping_cols), sep = " / ", remove = FALSE)
    p <- ggplot(plot_df, aes(x = pathway, y = value_for_plot, fill = sender)) +
      geom_col(width = 0.8) +
      facet_wrap(~ facet_group, scales = "free_y") +
      labs(x = NULL, y = "Summed communication contribution", fill = "Sender")
  }

  p <- p +
    theme_classic(base_size = 10) +
    theme(
      text = element_text(family = "Arial"),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      legend.position = "right"
    )

  ggsave(file.path(communication_figure_dir, paste0(out_prefix, ".png")), p, width = 10, height = 5.5, dpi = 600)
  ggsave(file.path(communication_figure_dir, paste0(out_prefix, ".pdf")), p, width = 10, height = 5.5)
  invisible(p)
}

# ------------------------------------------------------------
# Source-table-driven execution
# ------------------------------------------------------------

if (file.exists(file.path(communication_input_dir, "nichenet_top_ligands.csv"))) {
  nichenet_top_ligands <- readr::read_csv(file.path(communication_input_dir, "nichenet_top_ligands.csv"), show_col_types = FALSE)
  readr::write_csv(nichenet_top_ligands, file.path(communication_table_dir, "nichenet_top_ligands.csv"))
  plot_nichenet_top_ligands(nichenet_top_ligands)
}

if (file.exists(file.path(communication_input_dir, "nichenet_ligand_network_edges.csv"))) {
  nichenet_edges <- readr::read_csv(file.path(communication_input_dir, "nichenet_ligand_network_edges.csv"), show_col_types = FALSE)
  readr::write_csv(nichenet_edges, file.path(communication_table_dir, "nichenet_ligand_network_edges.csv"))
  plot_nichenet_edge_summary(nichenet_edges)
}

if (file.exists(file.path(communication_input_dir, "cellchat_core_pathway_sender_contributions.csv"))) {
  cellchat_contrib <- readr::read_csv(file.path(communication_input_dir, "cellchat_core_pathway_sender_contributions.csv"), show_col_types = FALSE)
  readr::write_csv(cellchat_contrib, file.path(communication_table_dir, "cellchat_core_pathway_sender_contributions.csv"))
  plot_cellchat_sender_contributions(cellchat_contrib)
}

# ------------------------------------------------------------
# Optional heavy re-run notice
# ------------------------------------------------------------
# Full NicheNet/CellChat re-execution is not called by default. The original
# scripts relied on high-memory objects, SCT matrices, package-version-specific
# CellChat/NicheNet behavior, and external priors/databases. They are archived at:
#   provenance/archived_scripts/MNPBurden_ligand_analysis_v1.R
#   provenance/archived_scripts/MNPBurden_ligand_analysis_v2.R
#   provenance/archived_scripts/nichenet_v1.R
#   provenance/archived_scripts/cellphonedb_v1.R

if (!any(expected_sources$present)) {
  message("No communication source tables found. Wrote expected-source manifest and skipped plotting.")
  message("Place exported communication source CSVs in: ", communication_input_dir)
}
