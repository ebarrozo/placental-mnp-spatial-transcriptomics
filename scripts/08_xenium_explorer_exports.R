# ============================================================
# Fig. 6 Xenium Explorer annotation exports
# Purpose:
#   Export public-object annotations for the four Xenium samples used
#   in Fig. 6. The figure itself was generated manually in Xenium
#   Explorer 4 from native 10x Xenium output directories.
#
# Inputs:
#   - obj_public loaded by scripts/02_load_public_object.R
#   - config/figure6_xenium_explorer_manifest.csv
#
# Outputs:
#   - data/derived/xenium_explorer_annotations/*_forXeniumExplorer.csv
#   - resolved Fig. 6 panel manifest
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

if (!exists("obj_public")) {
  stop("obj_public not found. Run scripts/02_load_public_object.R first.")
}

manifest_path <- file.path("config", "figure6_xenium_explorer_manifest.csv")
outdir <- file.path("data", "derived", "xenium_explorer_annotations")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(manifest_path)) {
  stop("Missing Fig. 6 manifest: ", manifest_path)
}

fig6_manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)

required_manifest_cols <- c(
  "panel", "delivery_class", "mnp_burden", "interpretation_marker",
  "sample_id_std", "xenium_output_directory", "xenium_region_label"
)
missing_manifest_cols <- setdiff(required_manifest_cols, colnames(fig6_manifest))
if (length(missing_manifest_cols) > 0) {
  stop("Fig. 6 manifest missing columns: ", paste(missing_manifest_cols, collapse = ", "))
}

md <- obj_public@meta.data %>%
  tibble::rownames_to_column("seurat_cell_id") %>%
  tibble::as_tibble()

if (!"sample_id_std" %in% colnames(md) && "sample_id_standardized" %in% colnames(md)) {
  md <- md %>% dplyr::mutate(sample_id_std = sample_id_standardized)
}

if (!"sample_id_std" %in% colnames(md)) {
  stop("Public object lacks sample_id_std or sample_id_standardized.")
}

# Xenium Explorer uses native cell IDs. In the original workflow, Seurat
# suffixes were removed with gsub("_[0-9]+$", "", Cells(object)).
normalize_xenium_cell_id <- function(x) {
  gsub("_[0-9]+$", "", x)
}

# Prefer cleaner labels when present, but export multiple compatible groupings.
annotation_cols <- c(
  "predicted.celltype",
  "nn.predicted.celltype",
  "celltype_display",
  "niches",
  "MNPBurden",
  "Delivery_class"
)
annotation_cols <- intersect(annotation_cols, colnames(md))

export_one_annotation <- function(df, sample_id, group_col) {
  out <- df %>%
    transmute(
      cell_id = normalize_xenium_cell_id(seurat_cell_id),
      group = as.character(.data[[group_col]])
    ) %>%
    filter(!is.na(cell_id), !is.na(group), group != "") %>%
    distinct(cell_id, .keep_all = TRUE)

  safe_group <- gsub("[^A-Za-z0-9]+", "_", group_col)
  outfile <- file.path(outdir, paste0(sample_id, "_", safe_group, "_forXeniumExplorer.csv"))
  readr::write_csv(out, outfile)

  tibble(
    sample_id_std = sample_id,
    annotation_column = group_col,
    exported_file = outfile,
    n_cells_exported = nrow(out),
    n_groups = dplyr::n_distinct(out$group)
  )
}

export_records <- list()
full_meta_records <- list()

for (sample_id in fig6_manifest$sample_id_std) {
  sample_df <- md %>% filter(sample_id_std == sample_id)

  if (nrow(sample_df) == 0) {
    warning("No cells found for Fig. 6 sample: ", sample_id)
    next
  }

  # Full public metadata table for audit/reuse.
  full_meta <- sample_df %>%
    mutate(cell_id = normalize_xenium_cell_id(seurat_cell_id)) %>%
    relocate(cell_id, seurat_cell_id)

  full_meta_file <- file.path(outdir, paste0(sample_id, "_full_public_metadata_forXeniumExplorer.csv"))
  readr::write_csv(full_meta, full_meta_file)

  full_meta_records[[sample_id]] <- tibble(
    sample_id_std = sample_id,
    full_metadata_file = full_meta_file,
    n_cells = nrow(full_meta)
  )

  for (group_col in annotation_cols) {
    export_records[[paste(sample_id, group_col, sep = "__")]] <-
      export_one_annotation(sample_df, sample_id, group_col)
  }
}

export_summary <- bind_rows(export_records)
full_meta_summary <- bind_rows(full_meta_records)

resolved_manifest <- fig6_manifest %>%
  left_join(
    md %>%
      count(sample_id_std, name = "n_public_object_cells"),
    by = "sample_id_std"
  ) %>%
  left_join(full_meta_summary, by = "sample_id_std")

readr::write_csv(
  resolved_manifest,
  file.path(outdir, "Fig6_XeniumExplorer_panel_manifest_resolved.csv")
)

readr::write_csv(
  export_summary,
  file.path(outdir, "Fig6_XeniumExplorer_annotation_export_summary.csv")
)

message("Fig. 6 Xenium Explorer annotation exports written to: ", outdir)
message("Reminder: exact Fig. 6 screenshots require opening the listed native Xenium output directories in Xenium Explorer 4.")
