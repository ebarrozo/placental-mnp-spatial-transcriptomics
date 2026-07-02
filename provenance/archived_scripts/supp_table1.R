## supp_table1.R
## run after meta_update.R (might be missing some helpers from MNPBurden_ligand_analysis_v2.R)

# ============================================================
# Supplementary Table 1 inventory
# Run this first and paste the console output
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
})

supp_dir <- "/home/ebarrozo/visiumhd_mnp/results/audra/supplementary_tables_final"
csv_dir  <- file.path(supp_dir, "csv_source_tables")

dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n===== CURRENT WORKING DIRECTORY =====\n")
print(getwd())

cat("\n===== OBJECTS IN CURRENT R SESSION =====\n")
objs <- ls(envir = .GlobalEnv)
print(objs)

object_inventory <- tibble(
  object_name = objs,
  class = vapply(objs, function(x) paste(class(get(x, envir = .GlobalEnv)), collapse = "; "), character(1)),
  nrow = vapply(objs, function(x) {
    obj <- get(x, envir = .GlobalEnv)
    if (is.data.frame(obj) || tibble::is_tibble(obj)) nrow(obj) else NA_integer_
  }, integer(1)),
  ncol = vapply(objs, function(x) {
    obj <- get(x, envir = .GlobalEnv)
    if (is.data.frame(obj) || tibble::is_tibble(obj)) ncol(obj) else NA_integer_
  }, integer(1))
)

print(object_inventory, n = Inf)

readr::write_csv(
  object_inventory,
  file.path(csv_dir, "R_session_object_inventory.csv")
)

cat("\n===== CANDIDATE CSV/TSV FILES UNDER AUDRA RESULTS =====\n")

candidate_roots <- c(
  "/home/ebarrozo/visiumhd_mnp/results/audra",
  "/home/ebarrozo/visiumhd_mnp/results/seurat_all_v1/seurat_all_ref_merged"
)

candidate_files <- purrr::map_dfr(candidate_roots, function(root) {
  if (!dir.exists(root)) return(tibble())
  tibble(
    root = root,
    file = list.files(
      root,
      pattern = "\\.(csv|tsv|txt)$",
      recursive = TRUE,
      full.names = TRUE
    )
  )
}) %>%
  mutate(
    basename = basename(file),
    likely_category = case_when(
      str_detect(basename, regex("metadata|model_ready|meta", ignore_case = TRUE)) ~ "metadata",
      str_detect(basename, regex("platform|celltype|composition|origin|maternal|NK", ignore_case = TRUE)) ~ "annotation_origin_NK",
      str_detect(basename, regex("module|contrast|loo|permutation", ignore_case = TRUE)) ~ "module_models",
      str_detect(basename, regex("DE|limma|pathway|enrich|GO|hallmark", ignore_case = TRUE)) ~ "limma_pathway",
      str_detect(basename, regex("macrophage|mac_state|marker", ignore_case = TRUE)) ~ "macrophage",
      str_detect(basename, regex("nichenet|ligand|target|CellChat|alluvial|communication", ignore_case = TRUE)) ~ "communication",
      str_detect(basename, regex("polymer|neighborhood|coenrichment|VisiumHD", ignore_case = TRUE)) ~ "exploratory_sensitivity",
      TRUE ~ "other"
    )
  )

candidate_files_short <- candidate_files %>%
  filter(likely_category != "other") %>%
  arrange(likely_category, basename)

print(candidate_files_short %>% select(likely_category, basename, file), n = 200)

readr::write_csv(
  candidate_files,
  file.path(csv_dir, "candidate_source_files_inventory.csv")
)

cat("\nInventory written to:\n")
cat(file.path(csv_dir, "R_session_object_inventory.csv"), "\n")
cat(file.path(csv_dir, "candidate_source_files_inventory.csv"), "\n")



# ============================================================
# Build Supplementary Table 1 workbook from current R session
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
  library(openxlsx)
})

# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------

supp_dir <- "/home/ebarrozo/visiumhd_mnp/results/audra/supplementary_tables_final"
csv_dir  <- file.path(supp_dir, "csv_source_tables")

dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)

xlsx_out <- file.path(
  supp_dir,
  "Supplementary_Table_1_MNP_spatial_analysis_outputs.xlsx"
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

get_obj <- function(object_name) {
  if (exists(object_name, envir = .GlobalEnv)) {
    get(object_name, envir = .GlobalEnv)
  } else {
    NULL
  }
}

as_table <- function(x, table_id, description, source_figure = NA_character_, source = NA_character_) {
  if (is.null(x)) {
    return(tibble(
      table_id = table_id,
      description = description,
      source_figure = source_figure,
      source = source,
      note = "OBJECT OR FILE NOT FOUND"
    ))
  }

  if (!is.data.frame(x) && !tibble::is_tibble(x)) {
    return(tibble(
      table_id = table_id,
      description = description,
      source_figure = source_figure,
      source = source,
      object_class = paste(class(x), collapse = "; "),
      note = "Object was not a data frame or tibble"
    ))
  }

  x <- x %>%
    as_tibble() %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::everything(),
        .fns = ~ {
          if (inherits(.x, "Date")) {
            as.character(.x)
          } else if (inherits(.x, "POSIXct") || inherits(.x, "POSIXt")) {
            as.character(.x)
          } else if (is.factor(.x)) {
            as.character(.x)
          } else if (is.list(.x)) {
            vapply(.x, function(z) paste(z, collapse = "; "), character(1))
          } else {
            as.character(.x)
          }
        }
      )
    )

  x %>%
    dplyr::mutate(
      table_id = table_id,
      description = description,
      source_figure = source_figure,
      source = source,
      .before = 1
    )
}

safe_read_table <- function(path, table_id, description, source_figure = NA_character_) {
  if (!file.exists(path)) {
    return(as_table(NULL, table_id, description, source_figure, path))
  }

  ext <- tolower(tools::file_ext(path))

  x <- if (ext == "tsv" || ext == "txt") {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    readr::read_csv(path, show_col_types = FALSE)
  }

  as_table(
    x,
    table_id = table_id,
    description = description,
    source_figure = source_figure,
    source = path
  )
}

safe_write_csv <- function(x, filename) {
  readr::write_csv(x, file.path(csv_dir, filename))
  invisible(file.path(csv_dir, filename))
}

write_sheet <- function(wb, sheet_name, df) {
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, df)
  openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)
  openxlsx::addFilter(wb, sheet_name, row = 1, cols = seq_len(ncol(df)))

  # Avoid insane widths on huge text columns
  widths <- rep(14, ncol(df))
  names(widths) <- colnames(df)
  widths[colnames(df) %in% c("description", "source", "source_file", "pathway", "term", "gene", "ligand", "target")] <- 28
  widths <- pmin(widths, 40)

  openxlsx::setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = widths)
}

# Avoid accidentally exporting 3.8M-cell metadata tables
compact_metadata <- function(x) {
  x <- as_tibble(x)

  keep_cols <- intersect(
    c(
      "sample_id_std", "study_id", "platform", "platform_clean", "platform_detail",
      "Delivery_class", "MNPBurden", "mnp_model_sample",
      "cumulative_abundance_ug_g", "cumulative_MNP_ug_g", "MNP_total_ug_g",
      "MNP_total", "MNPBurden_threshold", "MassivePVF", "Microcalc",
      "SDI_score", "SVI_score", "BMI", "maternal_age", "Prior_PTB",
      "HTN", "Hypertensive_Disorder", "Gestational_Age", "GA", "PTB",
      "Ethnicity_Race", "hospital"
    ),
    colnames(x)
  )

  if (length(keep_cols) == 0) return(x)

  x %>%
    select(all_of(keep_cols), everything())
}

# ------------------------------------------------------------
# 1a. Model-ready metadata
# ------------------------------------------------------------

metadata_obj <- dplyr::case_when(
  exists("meta_new", envir = .GlobalEnv) ~ "meta_new",
  exists("meta_new_subset", envir = .GlobalEnv) ~ "meta_new_subset",
  exists("meta_new_raw", envir = .GlobalEnv) ~ "meta_new_raw",
  TRUE ~ NA_character_
)

st1a <- if (!is.na(metadata_obj)) {
  get_obj(metadata_obj) %>%
    compact_metadata() %>%
    as_table(
      table_id = "ST1a",
      description = "Model-ready sample-level metadata used for spatial MNP analyses",
      source_figure = "Fig. 1 and all downstream analyses",
      source = paste0("R object: ", metadata_obj)
    )
} else {
  as_table(NULL, "ST1a", "Model-ready metadata", "Fig. 1", "R object not found")
}

safe_write_csv(st1a, "ST1a_model_ready_metadata.csv")

# ------------------------------------------------------------
# 1b. Cohort / platform / annotation support
# ------------------------------------------------------------

st1b <- bind_rows(
  as_table(
    get_obj("all_sample_keys"),
    "ST1b.1",
    "Sample-level key linking sample_id_std, study_id, platform, and related identifiers",
    "Fig. 1",
    "R object: all_sample_keys"
  ),
  as_table(
    get_obj("platform_celltype_counts_detail"),
    "ST1b.2",
    "Predicted cell-type counts and fractions by platform and assay version",
    "Supplementary Fig. 1",
    "R object: platform_celltype_counts_detail"
  ),
  as_table(
    get_obj("platform_cor_detail_df"),
    "ST1b.3",
    "Spearman correlations of predicted cell-type fractions across platform versions",
    "Supplementary Fig. 1",
    "R object: platform_cor_detail_df"
  ),
  as_table(
    get_obj("qa_platform_mismatch"),
    "ST1b.4",
    "Platform metadata mismatch QA report",
    "Supplementary Fig. 1",
    "R object: qa_platform_mismatch"
  )
)

safe_write_csv(st1b, "ST1b_cohort_platform_annotation.csv")

# ------------------------------------------------------------
# 1c. Cell-type composition, origin, and NK support
# ------------------------------------------------------------

st1c <- bind_rows(
  as_table(
    get_obj("celltype_composition_contrasts"),
    "ST1c.1",
    "Within-delivery high-versus-low MNP predicted cell-type composition contrasts",
    "Fig. 2",
    "R object: celltype_composition_contrasts"
  ),
  as_table(
    get_obj("module_celltype_correlations"),
    "ST1c.2",
    "Correlations between module scores and predicted cell-type composition",
    "Fig. 2 and Supplementary Fig. 1",
    "R object: module_celltype_correlations"
  ),
  as_table(
    get_obj("origin_by_celltype"),
    "ST1c.3",
    "Predicted maternal/fetal-placental origin counts and fractions by predicted cell type",
    "Supplementary Fig. 2",
    "R object: origin_by_celltype"
  ),
  as_table(
    get_obj("immune_origin_focus"),
    "ST1c.4",
    "Focused maternal-origin summary for NK and T-cell annotations",
    "Supplementary Fig. 2",
    "R object: immune_origin_focus"
  ),
  as_table(
    get_obj("nk_subject_origin"),
    "ST1c.5",
    "Subject-level maternal-origin fractions among NK-cell observations",
    "Supplementary Fig. 2",
    "R object: nk_subject_origin"
  ),
  as_table(
    get_obj("nk_scores_subject_mnp"),
    "ST1c.6",
    "Subject-level NK marker-program scores among MNP-model samples",
    "Supplementary Fig. 3",
    "R object: nk_scores_subject_mnp"
  ),
  as_table(
    get_obj("nk_program_stats"),
    "ST1c.7",
    "Within-delivery high-versus-low MNP contrasts for NK marker-program scores",
    "Supplementary Fig. 3",
    "R object: nk_program_stats"
  )
)

safe_write_csv(st1c, "ST1c_celltype_origin_NK_support.csv")

# ------------------------------------------------------------
# 1d. Module score models and robustness
# ------------------------------------------------------------

st1d <- bind_rows(
  as_table(
    get_obj("sample_scores_z"),
    "ST1d.1",
    "Sample-level platform-normalized module scores",
    "Fig. 3",
    "R object: sample_scores_z"
  ),
  as_table(
    get_obj("subject_program_scores"),
    "ST1d.2",
    "Subject-level macrophage program scores or spatial module-program scores, depending on analysis context",
    "Fig. 3 and Fig. 5",
    "R object: subject_program_scores"
  ),
  as_table(
    get_obj("within_delivery_contrasts"),
    "ST1d.3",
    "Within-delivery high-versus-low MNP module-score contrasts",
    "Fig. 3",
    "R object: within_delivery_contrasts"
  ),
  as_table(
    get_obj("subject_model_results_z"),
    "ST1d.4",
    "Subject-level platform-normalized MNP burden interaction model results",
    "Fig. 3",
    "R object: subject_model_results_z"
  ),
  as_table(
    get_obj("loo_mnp_summary"),
    "ST1d.5",
    "Leave-one-subject-out summary for primary MNP-burden module contrasts",
    "Supplementary Fig. 4",
    "R object: loo_mnp_summary"
  ),
  as_table(
    get_obj("loo_mnp_contrasts"),
    "ST1d.6",
    "Leave-one-subject-out high-versus-low MNP module contrasts",
    "Supplementary Fig. 4",
    "R object: loo_mnp_contrasts"
  ),
  as_table(
    get_obj("exact_perm_summary"),
    "ST1d.7",
    "Exact delivery-preserving permutation summary for primary module contrasts",
    "Supplementary Fig. 4",
    "R object: exact_perm_summary"
  ),
  as_table(
    get_obj("exact_perm_results"),
    "ST1d.8",
    "Exact permutation results for primary module contrasts",
    "Supplementary Fig. 4",
    "R object: exact_perm_results"
  ),
  as_table(
    get_obj("module_gene_check"),
    "ST1d.9",
    "Module gene presence and scoring QA",
    "Fig. 3",
    "R object: module_gene_check"
  )
)

safe_write_csv(st1d, "ST1d_module_score_models.csv")

# ------------------------------------------------------------
# 1e. Limma pseudobulk and pathway enrichment
# ------------------------------------------------------------

# Keep full de_all_clean out of Excel unless it is reasonably sized.
# Export full clean DE as CSV only, then include a filtered/top version in Excel.

if (exists("de_all_clean", envir = .GlobalEnv)) {
  de_full <- get_obj("de_all_clean")
  readr::write_csv(
    de_full,
    file.path(csv_dir, "ST1e_full_de_all_clean_subject_aware_limma.csv")
  )

  de_export <- de_full %>%
    as_tibble() %>%
    mutate(
      sig_q05 = if ("adj.P.Val" %in% colnames(.)) adj.P.Val < 0.05 else NA,
      sig_q10 = if ("adj.P.Val" %in% colnames(.)) adj.P.Val < 0.10 else NA
    ) %>%
    filter(
      (sig_q05 %in% TRUE) |
        (sig_q10 %in% TRUE & abs(logFC) >= 0.25)
    )

  if (nrow(de_export) > 50000) {
    de_export <- de_export %>%
      group_by(Delivery_class, predicted.celltype) %>%
      arrange(adj.P.Val, desc(abs(logFC)), .by_group = TRUE) %>%
      slice_head(n = 1000) %>%
      ungroup()
  }
} else {
  de_export <- NULL
}

# Try to read key pathway enrichment files directly from known folder
base_de <- "/home/ebarrozo/visiumhd_mnp/results/seurat_all_v1/seurat_all_ref_merged/de_MNPBurden_Delivery_class_bulk_sensitive_go"

pathway_files <- c(
  file.path(base_de, "pb_Term_MNPBurden_limma_SCT/cp_results/GO_BP_UP_NK Cell.csv"),
  file.path(base_de, "pb_Term_MNPBurden_limma_SCT/cp_results/MSIG_HALLMARK_UP_NK Cell.csv"),
  file.path(base_de, "pb_Preterm_MNPBurden_limma_SCT/cp_results/GO_BP_UP_Macrophage.csv"),
  file.path(base_de, "pb_Preterm_MNPBurden_limma_SCT/cp_results/MSIG_HALLMARK_UP_Macrophage.csv")
)

pathway_tables <- purrr::imap_dfr(
  pathway_files,
  ~ safe_read_table(
    .x,
    table_id = paste0("ST1e.pathway.", .y),
    description = "Pathway enrichment table used for Fig. 3b/c when available",
    source_figure = "Fig. 3b/c"
  )
)

st1e <- bind_rows(
  as_table(
    get_obj("fig3a_long"),
    "ST1e.1",
    "Long-format Fig. 3a heatmap input with limma log2 fold-change values",
    "Fig. 3a",
    "R object: fig3a_long"
  ),
  as_table(
    de_export,
    "ST1e.2",
    "Filtered subject-aware limma pseudobulk DE results; full table exported as CSV source table when available",
    "Fig. 3a-c",
    "R object: de_all_clean; filtered for Excel"
  ),
  as_table(
    get_obj("term_de"),
    "ST1e.3",
    "Term subject-aware limma pseudobulk DE results",
    "Fig. 3",
    "R object: term_de"
  ),
  as_table(
    get_obj("preterm_de") %>% { if (!is.null(.)) head(., 50000) else NULL },
    "ST1e.4",
    "Preterm subject-aware limma pseudobulk DE results; truncated to first 50,000 rows for Excel if full table is large",
    "Fig. 3",
    "R object: preterm_de"
  ),
  pathway_tables
)

safe_write_csv(st1e, "ST1e_limma_pathway_outputs.csv")

# ------------------------------------------------------------
# 1f. Macrophage state analysis
# ------------------------------------------------------------

st1f <- bind_rows(
  as_table(
    get_obj("state_frac_subject"),
    "ST1f.1",
    "Subject-level macrophage state fractions by delivery class and MNP burden",
    "Fig. 5",
    "R object: state_frac_subject"
  ),
  as_table(
    get_obj("state_frac_stats"),
    "ST1f.2",
    "Within-delivery high-versus-low MNP contrasts for macrophage state fractions",
    "Fig. 5",
    "R object: state_frac_stats"
  ),
  as_table(
    get_obj("top_mac_state_markers"),
    "ST1f.3",
    "Selected positive marker genes for macrophage states used in dot plot",
    "Fig. 5",
    "R object: top_mac_state_markers"
  ),
  as_table(
    get_obj("mac_state_markers"),
    "ST1f.4",
    "FindAllMarkers results for macrophage state annotations",
    "Fig. 5",
    "R object: mac_state_markers"
  ),
  as_table(
    get_obj("polarization_qc"),
    "ST1f.5",
    "Macrophage polarization/state annotation QC table",
    "Fig. 5",
    "R object: polarization_qc"
  )
)

safe_write_csv(st1f, "ST1f_macrophage_state_analysis.csv")

# ------------------------------------------------------------
# 1g. Communication analysis
# ------------------------------------------------------------

# Use final CellChat alluvial CSVs if present.
cellchat_final_dir <- "/home/ebarrozo/visiumhd_mnp/results/audra/MNP_interactomes_CellChat_v12_full_min500_png/Figure_addons/alluvial_AB_core12_ordered_publication_final"

st1g <- bind_rows(
  safe_read_table(
    file.path(cellchat_final_dir, "final_alluvial_A_sender_by_pathway_Preterm_High.csv"),
    "ST1g.1",
    "CellChat sender-to-pathway strengths for preterm high-MNP placentae",
    "Supplementary Fig. 5a"
  ),
  safe_read_table(
    file.path(cellchat_final_dir, "final_alluvial_A_sender_by_pathway_Term_High.csv"),
    "ST1g.2",
    "CellChat sender-to-pathway strengths for term high-MNP placentae",
    "Supplementary Fig. 5b"
  ),
  safe_read_table(
    file.path(cellchat_final_dir, "final_alluvial_B_pathway_by_receiver_Preterm_High.csv"),
    "ST1g.3",
    "CellChat pathway-to-receiver strengths for preterm high-MNP placentae; exported for provenance",
    "Supplementary Fig. 5 source"
  ),
  safe_read_table(
    file.path(cellchat_final_dir, "final_alluvial_B_pathway_by_receiver_Term_High.csv"),
    "ST1g.4",
    "CellChat pathway-to-receiver strengths for term high-MNP placentae; exported for provenance",
    "Supplementary Fig. 5 source"
  )
)

safe_write_csv(st1g, "ST1g_communication_analysis.csv")

# ------------------------------------------------------------
# 1h. Exploratory and sensitivity analyses
# ------------------------------------------------------------

st1h <- bind_rows(
  as_table(
    get_obj("polymer_missingness"),
    "ST1h.1",
    "Polymer missingness and availability summary",
    "Supplementary Fig. 4",
    "R object: polymer_missingness"
  ),
  as_table(
    get_obj("polymer_module_slopes"),
    "ST1h.2",
    "Exploratory polymer-specific module-score slopes by delivery class",
    "Supplementary Fig. 4",
    "R object: polymer_module_slopes"
  ),
  as_table(
    get_obj("polymer_module_model_terms"),
    "ST1h.3",
    "Exploratory polymer-specific module model terms",
    "Supplementary Fig. 4",
    "R object: polymer_module_model_terms"
  ),
  as_table(
    get_obj("loo_polymer_summary"),
    "ST1h.4",
    "Leave-one-out summary for polymer-specific module associations",
    "Supplementary Fig. 4",
    "R object: loo_polymer_summary"
  ),
  as_table(
    get_obj("loo_polymer_slopes"),
    "ST1h.5",
    "Leave-one-out polymer-specific slopes",
    "Supplementary Fig. 4",
    "R object: loo_polymer_slopes"
  ),
  as_table(
    get_obj("within_delivery_neighborhood_contrasts"),
    "ST1h.6",
    "Focused Xenium neighborhood high-versus-low MNP contrasts",
    "Supplementary sensitivity analysis",
    "R object: within_delivery_neighborhood_contrasts"
  ),
  as_table(
    get_obj("within_delivery_coenrichment_contrasts"),
    "ST1h.7",
    "Visium HD local module co-enrichment high-versus-low MNP contrasts",
    "Supplementary sensitivity analysis",
    "R object: within_delivery_coenrichment_contrasts"
  ),
  as_table(
    get_obj("subject_level_coenrichment"),
    "ST1h.8",
    "Subject-level Visium HD local module co-enrichment summaries",
    "Supplementary sensitivity analysis",
    "R object: subject_level_coenrichment"
  )
)

safe_write_csv(st1h, "ST1h_exploratory_sensitivity.csv")

# ------------------------------------------------------------
# 1i. Data dictionary
# ------------------------------------------------------------

st1i <- tibble::tribble(
  ~sheet, ~table_id, ~description,
  "1a_Model_ready_metadata", "ST1a", "Model-ready metadata used for sample/subject inclusion and downstream MNP burden analyses.",
  "1b_Cohort_platform_annotation", "ST1b.1-ST1b.4", "Cohort, platform, and cross-platform annotation support tables.",
  "1c_Celltype_origin_NK", "ST1c.1-ST1c.7", "Predicted cell-type composition, predicted maternal-origin assignment, and NK marker-program support.",
  "1d_Module_score_models", "ST1d.1-ST1d.9", "Spatial module scores, high-versus-low MNP contrasts, and robustness analyses.",
  "1e_Limma_pathway_outputs", "ST1e.1-ST1e.pathway", "Subject-aware limma pseudobulk DE and pathway enrichment outputs supporting Fig. 3.",
  "1f_Macrophage_state_analysis", "ST1f.1-ST1f.5", "Macrophage state fractions, high-versus-low contrasts, and state marker genes.",
  "1g_Communication_analysis", "ST1g.1-ST1g.4", "CellChat alluvial source tables and related communication provenance.",
  "1h_Exploratory_sensitivity", "ST1h.1-ST1h.8", "Exploratory polymer-specific, neighborhood, and Visium HD co-enrichment sensitivity analyses.",
  "1i_Data_dictionary", "ST1i", "This data dictionary."
)

safe_write_csv(st1i, "ST1i_data_dictionary.csv")

# ------------------------------------------------------------
# Build Excel workbook
# ------------------------------------------------------------

wb <- openxlsx::createWorkbook()

write_sheet(wb, "1a_Model_ready_metadata", st1a)
write_sheet(wb, "1b_Cohort_platform", st1b)
write_sheet(wb, "1c_Celltype_origin_NK", st1c)
write_sheet(wb, "1d_Module_score_models", st1d)
write_sheet(wb, "1e_Limma_pathway", st1e)
write_sheet(wb, "1f_Macrophage_states", st1f)
write_sheet(wb, "1g_Communication", st1g)
write_sheet(wb, "1h_Exploratory_sens", st1h)
write_sheet(wb, "1i_Data_dictionary", st1i)

openxlsx::saveWorkbook(wb, xlsx_out, overwrite = TRUE)

# ------------------------------------------------------------
# Missing object/file report
# ------------------------------------------------------------

missing_report <- bind_rows(st1a, st1b, st1c, st1d, st1e, st1f, st1g, st1h) %>%
  filter(str_detect(coalesce(note, ""), "NOT FOUND|NULL|not a data"))

readr::write_csv(
  missing_report,
  file.path(csv_dir, "missing_or_placeholder_tables_report.csv")
)

cat("\nSaved workbook:\n", xlsx_out, "\n")
cat("\nSaved source CSVs in:\n", csv_dir, "\n")
cat("\nMissing/placeholder report:\n", file.path(csv_dir, "missing_or_placeholder_tables_report.csv"), "\n")






# ============================================================
# QC Supplementary Table 1 source CSVs
# ============================================================

library(readr)
library(dplyr)
library(purrr)
library(stringr)

csv_dir <- "/home/ebarrozo/visiumhd_mnp/results/audra/supplementary_tables_final/csv_source_tables"

csv_files <- list.files(csv_dir, pattern = "^ST1.*\\.csv$", full.names = TRUE)

qc <- purrr::map_dfr(csv_files, function(f) {
  x <- readr::read_csv(f, show_col_types = FALSE)
  tibble(
    file = basename(f),
    n_rows = nrow(x),
    n_cols = ncol(x),
    table_ids = paste(unique(x$table_id), collapse = "; ")
  )
})

print(qc, n = Inf)

readr::write_csv(
  qc,
  file.path(csv_dir, "ST1_source_table_QC_summary.csv")
)


# ============================================================
# Patch ST1e: make Excel sheet lean, keep full DE as CSV only
# ============================================================

library(dplyr)
library(readr)
library(stringr)
library(tibble)

csv_dir <- "/home/ebarrozo/visiumhd_mnp/results/audra/supplementary_tables_final/csv_source_tables"

# Full DE remains here:
# ST1e_full_de_all_clean_subject_aware_limma.csv

de_full <- readr::read_csv(
  file.path(csv_dir, "ST1e_full_de_all_clean_subject_aware_limma.csv"),
  show_col_types = FALSE
)

# Make a compact Excel-facing DE table
de_excel <- de_full %>%
  as_tibble() %>%
  mutate(
    abs_logFC = abs(logFC),
    sig_rank = if_else(!is.na(adj.P.Val), adj.P.Val, P.Value)
  ) %>%
  group_by(Delivery_class, predicted.celltype) %>%
  arrange(sig_rank, desc(abs_logFC), .by_group = TRUE) %>%
  slice_head(n = 250) %>%
  ungroup() %>%
  mutate(
    table_id = "ST1e.2",
    description = "Top subject-aware limma pseudobulk DE results by delivery class and predicted cell type; full DE table is provided as a separate CSV source file",
    source_figure = "Fig. 3a-c",
    source = "ST1e_full_de_all_clean_subject_aware_limma.csv",
    .before = 1
  )

# Pull Fig. 3a and pathway rows from existing ST1e table, if already made
st1e_old <- readr::read_csv(
  file.path(csv_dir, "ST1e_limma_pathway_outputs.csv"),
  show_col_types = FALSE
)

st1e_keep <- st1e_old %>%
  filter(
    table_id == "ST1e.1" |
      str_detect(table_id, "^ST1e\\.pathway")
  )

st1e_new <- bind_rows(
  st1e_keep,
  de_excel %>% mutate(across(everything(), as.character))
)

readr::write_csv(
  st1e_new,
  file.path(csv_dir, "ST1e_limma_pathway_outputs_COMPACT_FOR_EXCEL.csv")
)

cat("Compact ST1e rows:", nrow(st1e_new), "\n")