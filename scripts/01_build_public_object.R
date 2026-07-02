# ============================================================
# Build GEO/publication-facing Seurat object from obj_forced
# SAFE: does not modify obj_forced in place
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

# Usage options:
# 1) Interactive: create obj_forced in your R session, then source this script.
# 2) File-based: set OBJ_FORCED_RDS=/path/to/obj_forced.rds before running.
# Optional: set OUTDIR=/path/to/output_dir. Default: data/external.

obj_forced_rds <- Sys.getenv("OBJ_FORCED_RDS", unset = NA_character_)
outdir <- Sys.getenv("OUTDIR", unset = file.path("data", "external"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!exists("obj_forced")) {
  if (is.na(obj_forced_rds) || !nzchar(obj_forced_rds) || !file.exists(obj_forced_rds)) {
    stop("Provide obj_forced in the R session or set OBJ_FORCED_RDS to an existing .rds file.")
  }
  obj_forced <- readRDS(obj_forced_rds)
}

stopifnot(is(obj_forced, "Seurat"))

# ----------------------------
# 1. Make explicit backup before any derived object is created
# ----------------------------
backup_path <- file.path(outdir, "obj_forced_BACKUP_before_geo_cleanup.rds")
if (!file.exists(backup_path)) {
  saveRDS(obj_forced, backup_path, compress = "xz")
}

# ----------------------------
# 2. Define curated GEO/public metadata columns
# ----------------------------
keep_cols <- c(
  # Object / sample identity
  "orig.ident", "sample", "sample_id_std", "sample_id", "study_id", "seurat_id", "seq_id",
  "Subject.ID", "Pregnancy.ID",

  # Count QC
  "nCount_RNA", "nFeature_RNA",

  # Platform / sequencing metadata
  "platform", "platform_clean", "platform_inferred",
  "xen_id_seurat", "vis_v1_seurat", "vis_hd_seurat",
  "xenium_id_seq", "visium_v1_seq", "visium_hd_seq",
  "xenium_platform", "visium_platform", "novel_spatial", "dataset", "multiple_spatial_ro_is",
  "xenium5k_run_id", "xenium5k_run_date", "xenium5k_slide_id", "xenium5k_block",
  "xenium5k_slide", "xenium5k_run",

  # Cell annotation
  "predicted.celltype", "predicted.celltype.score", "predicted.celltype_orig", "nn.predicted.celltype",
  "predicted.celltype_support", "celltype_display", "cell_type", "celltype_clean", "label_source",
  "origin_rule_based", "origin_confidence", "origin_reason", "origin_sex_supported",

  # Clinical / phenotype metadata
  "GA", "GA_continuous", "Delivery_class", "GA_group", "birth_weight", "Birth.weight",
  "infant_gender", "Infant.gender", "Ethnicity_Race", "Preeclampsia_bin", "Infections_any",
  "Prior_PTB", "Prior_PTB_TRUE", "Maternal_age_years", "BMI", "SDI_score_from_zip", "SVI_score",
  "Py_GC_MS", "mnp_model_sample", "mnp_model_subject", "histology_model_sample",

  # MNP burden, harmonized lowercase/public columns
  "MNPBurden", "MNPBurden_final", "MNPBurden_Delivery_class", "high_mnp_binary", "mnp_threshold_ug_g",
  "cumulative_mnp_ug_per_g", "cumulative_mnp_log10ugpg",
  "pe_ug_per_g", "sbr_ug_per_g", "pp_ug_per_g", "pvc_ug_per_g", "n66_ug_per_g", "pet_ug_per_g",
  "n6_ug_per_g", "pmma_ug_per_g", "pc_ug_per_g", "pu_ug_per_g", "abs_ug_per_g", "ps_ug_per_g",
  "pe_log10ugpg", "sbr_log10ugpg", "pp_log10ugpg", "pvc_log10ugpg", "n66_log10ugpg",
  "pet_log10ugpg", "n6_log10ugpg", "pmma_log10ugpg", "pc_log10ugpg", "pu_log10ugpg",
  "abs_log10ugpg", "ps_log10ugpg",

  # Histology covariates
  "calcifications_number", "pvf_percent_coverage", "Microcalc_presence", "Massive_PVF",

  # Manuscript module scores
  "ECM_Remodeling_score", "Vascular_Interface_score", "Macrophage_Remodeling_score",
  "NK_Cytotoxic_Regulatory_score", "Trophoblast_Remodeling_score",

  # Optional sex QC scores
  "XIST_score", "Y_score"
)

keep_cols_present <- intersect(keep_cols, colnames(obj_forced@meta.data))
keep_cols_missing <- setdiff(keep_cols, colnames(obj_forced@meta.data))

# ----------------------------
# 3. Create public object copy
# ----------------------------
obj_geo <- obj_forced
obj_geo@meta.data <- obj_geo@meta.data[, keep_cols_present, drop = FALSE]

# ----------------------------
# 4. Rename coded fields to public names
# Only retain these if deidentified/coded.
# ----------------------------
rename_map <- c(
  "Subject.ID" = "donor_id",
  "Pregnancy.ID" = "pregnancy_id",
  "sample_id_std" = "sample_id_standardized",
  "GA" = "gestational_age_weeks",
  "GA_continuous" = "gestational_age_weeks_continuous",
  "birth_weight" = "birth_weight_g",
  "Maternal_age_years" = "maternal_age_years",
  "BMI" = "maternal_bmi"
)
rename_map <- rename_map[names(rename_map) %in% colnames(obj_geo@meta.data)]
for (old_name in names(rename_map)) {
  colnames(obj_geo@meta.data)[colnames(obj_geo@meta.data) == old_name] <- rename_map[[old_name]]
}

# ----------------------------
# 5. Save audit files documenting what happened
# ----------------------------
geo_keep_audit <- tibble(
  column_name_requested = keep_cols,
  present_in_obj_forced = keep_cols %in% colnames(obj_forced@meta.data)
)
write_csv(geo_keep_audit, file.path(outdir, "geo_public_object_keep_columns_audit.csv"))
write_csv(tibble(column_name_missing_from_obj_forced = keep_cols_missing),
          file.path(outdir, "geo_public_object_requested_columns_missing.csv"))
write_csv(tibble(column_name_in_obj_geo = colnames(obj_geo@meta.data)),
          file.path(outdir, "geo_public_object_final_metadata_columns.csv"))

# ----------------------------
# 6. Save public object and MD5 checksum
# ----------------------------
geo_rds_path <- file.path(outdir, "spatial_mnp_publication_processed_seurat_v1.rds")
saveRDS(obj_geo, geo_rds_path, compress = "xz")
md5 <- tools::md5sum(geo_rds_path)
write_csv(tibble(file = geo_rds_path, md5 = as.character(md5)),
          file.path(outdir, "spatial_mnp_publication_processed_seurat_v1_md5.csv"))

cat("\nSaved GEO/publication-facing object:\n", geo_rds_path, "\n\n", sep = "")
cat("Metadata columns retained:", ncol(obj_geo@meta.data), "\n")
cat("Requested keep columns missing:", length(keep_cols_missing), "\n")
cat("MD5:", as.character(md5), "\n")
