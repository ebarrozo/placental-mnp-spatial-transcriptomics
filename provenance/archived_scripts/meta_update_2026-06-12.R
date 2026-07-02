## meta_update.R
# 6.12.26
# server 3, Enrico Barrozo, PhD

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(ggplot2)
mypal3 <-get_palette("igv",40)

# -------------------- USER CONFIG --------------------
setwd("/home/ebarrozo/visiumhd_mnp/results/audra")
# final is located here -> "/home/ebarrozo/visiumhd_mnp/results/audra/metadata_update_outputs/spatial_object_metadata_updated_mnp_qc.rds"))

obj <- readRDS("/home/ebarrozo/visiumhd_mnp/results/audra/mnp_obj_fullcells_maternal_fetal_rule_based.rds")
DefaultAssay(obj) <- "RNA"

Idents(obj) <- "nn.predicted.celltype"
levels(obj)
# [1] "Syncytiotrophoblast"      "Villous Cytotrophoblast" 
# [3] "Extravillous Trophoblast" "Fibroblast"              
# [5] "Perivascular Cell"        "Placental Endothelial"   
# [7] "Hofbauer Cell"            "Macrophage"              
# [9] "NK Cell"                  "Decidual NK Cell"        
#[11] "ILC3"                     "Plasma Cell"             
#[13] "B Cell"                   "T Cell"                  
#[15] "Treg"                     "RCTD_lowconf"  
Idents(obj) <- "platform"
levels(obj)
# [1] "Xenium"   "VisiumHD" "VisiumV1"
Idents(obj) <- "sample_id_std"
levels(obj)
# [1] "HP1_xen"   "HP13a_v1"  "HP13b_v1"  "HP13b_xen" "HP2_hd"    "HP2_xen"  
# [7] "HT1_v1"    "HT1_xen"   "HT2_hd"    "HT2_v1"    "HT2a_xen"  "HT2b_xen" 
#[13] "HT3_hd"    "HT3_xen"   "LP1_hd"    "LP1_xen"   "LP2_xen"   "LT1_xen"  
#[19] "LT2_xen"   "NC1a_v1"   "NC1b_v1"   "NC1c_v1"   "NC2_v1"    "NC3_v1"   
#[25] "NC4_v1"    "NC4_xen"   "ND5_v1"    "ND6_v1"    "ND7_v1"    "SP10_v1"  
#[31] "SP11_v1"   "SP12_v1"   "SP8_v1"    "SP9_v1" 
Idents(obj) <- "study_id"
levels(obj)
#  [1] "HP1"  "HP13" "HP2"  "HT1"  "HT2"  "HT3"  "LP1"  "LP2"  "LT1"  "NC4" 
#[11] "LT2"  "NC1"  "NC2"  "NC3"  "ND5"  "ND6"  "ND7"  "SP10" "SP11" "SP12"
#[21] "SP8"  "SP9" 





# ============================================================
# MNP Spatial Transcriptomics Metadata Update
# Purpose:
#   1. Clean old metadata artifacts
#   2. Join updated sample-level metadata by sample_id_std
#   3. Preserve Seurat/RCTD/origin/QC metadata
#   4. Create analysis-ready metadata columns
#   5. Save updated Seurat object and QA tables
#
# Primary exposure definition:
#   High vs low MNP burden by cumulative abundance threshold = 146.3 ug/g
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(janitor)
  library(tidyr)
})

# -----------------------------
# User inputs
# -----------------------------
meta_file <- "meta_model_ready_sample_id_std.csv"

outdir <- "metadata_update_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Optional: set to TRUE if you want to drop Integration entirely.
# The script archives it first as Integration_prev before dropping.
drop_integration_active <- TRUE

# Fixed exposure threshold
mnp_threshold_ug_g <- 146.3

# -----------------------------
# 0. Backup current object
# -----------------------------
saveRDS(obj, file.path(outdir, "spatial_object_before_metadata_update.rds"))

# -----------------------------
# 1. Capture current metadata
# -----------------------------
obj_meta0 <- obj@meta.data %>%
  rownames_to_column("cell_barcode")

stopifnot("sample_id_std" %in% colnames(obj_meta0))

cat("\nInitial metadata dimensions:\n")
print(dim(obj_meta0))

cat("\nInitial metadata columns:\n")
print(colnames(obj_meta0))

# -----------------------------
# 2. Archive or remove old artifacts
# -----------------------------
# Keep nn.predicted.celltype because it is your forced KNN annotation.
# Do NOT drop:
#   nn.predicted.celltype
#   predicted.celltype
#   predicted.celltype.score
#   predicted.celltype_support
#   celltype_clean
#   origin_* columns
#   nCount_RNA/nFeature_RNA/nCount_SCT/nFeature_SCT
#   sample/chunk/seq_id

old_join_artifacts <- c(
  "sample_id_std_prev",
  "platform_prev",
  "sample_id_std_prev2",
  "platform_prev2",
  "sample_id_prev2",
  "seurat_id_prev2"
)

old_join_artifacts_present <- intersect(old_join_artifacts, colnames(obj_meta0))

cat("\nOld join artifact columns to drop:\n")
print(old_join_artifacts_present)

obj_meta1 <- obj_meta0 %>%
  select(-any_of(old_join_artifacts_present))

# Archive Integration if present.
if ("Integration" %in% colnames(obj_meta1) && !"Integration_prev" %in% colnames(obj_meta1)) {
  obj_meta1 <- obj_meta1 %>%
    mutate(Integration_prev = Integration)
}

if (drop_integration_active && "Integration" %in% colnames(obj_meta1)) {
  obj_meta1 <- obj_meta1 %>%
    select(-Integration)
}

# Optional archive old platform fields before updated join.
if ("Platform" %in% colnames(obj_meta1) && !"Platform_prev" %in% colnames(obj_meta1)) {
  obj_meta1 <- obj_meta1 %>%
    mutate(Platform_prev = Platform)
}

# Remove old active Platform if you do not want ambiguity with platform/platform_clean.
# This preserves Platform_prev.
if ("Platform" %in% colnames(obj_meta1)) {
  obj_meta1 <- obj_meta1 %>%
    select(-Platform)
}

# -----------------------------
# 3. Load and standardize new metadata
# -----------------------------
meta_new_raw <- read_csv(meta_file, show_col_types = FALSE)

cat("\nNew metadata dimensions:\n")
print(dim(meta_new_raw))

cat("\nNew metadata columns before cleaning:\n")
print(colnames(meta_new_raw))

meta_new <- meta_new_raw %>%
  clean_names()

# janitor will convert some names.
# Standardize critical columns explicitly.
# Current expected cleaned names include:
#   sample_id_std
#   platform
#   study_id
#   delivery_class
#   mnp_burden
#   mnp_burden_delivery_class
#   cumulative_mn_ps_ug_per_g
#   cumulative_mn_ps_log10ugpg
#   calcifications_number
#   perivillous_fibrin_pvf_percent_coverage

# Rename awkward janitor conversions if present.
meta_new <- meta_new %>%
  rename(
    cumulative_mnp_ug_per_g =
      any_of(c("cumulative_mn_ps_ug_per_g", "cumulative_mnps_ug_per_g")),
    cumulative_mnp_log10ugpg =
      any_of(c("cumulative_mn_ps_log10ugpg", "cumulative_mnps_log10ugpg")),
    infant_gender =
      any_of(c("infant_gender")),
    birth_weight =
      any_of(c("birth_weight")),
    calcifications_number =
      any_of(c("calcifications_number")),
    pvf_percent_coverage =
      any_of(c("perivillous_fibrin_pvf_percent_coverage"))
  )

# If rename(any_of()) did not find a column, it is harmless.
# Now make clean analysis fields.
meta_new <- meta_new %>%
  mutate(
    sample_id_std = as.character(sample_id_std),
    study_id = as.character(study_id),

    # sample_id should be a human-readable sequencing/sample identifier.
    # Prefer sample_id if present, otherwise SampleID-derived cleaned name.
    sample_id = dplyr::coalesce(
      as.character(.data$sample_id),
      as.character(.data$sample_id_std)
    ),

    # Infer platform from sample_id_std suffix.
    platform_inferred = case_when(
      str_detect(sample_id_std, "_xen$") ~ "Xenium",
      str_detect(sample_id_std, "_hd$")  ~ "VisiumHD",
      str_detect(sample_id_std, "_v1$")  ~ "VisiumV1",
      TRUE ~ NA_character_
    ),

    # Trust suffix when available, but preserve original CSV platform too.
    platform_csv = as.character(platform),
    platform_clean = coalesce(platform_inferred, platform_csv),

    platform = platform_clean,

    Delivery_class = case_when(
      str_to_lower(delivery_class) == "term" ~ "Term",
      str_to_lower(delivery_class) == "preterm" ~ "Preterm",
      TRUE ~ as.character(delivery_class)
    ),

    MNPBurden = case_when(
      str_to_lower(mnp_burden) == "high" ~ "High",
      str_to_lower(mnp_burden) == "low" ~ "Low",
      TRUE ~ as.character(mnp_burden)
    ),

    # Enforce locked threshold if MNP burden missing but cumulative burden is present.
    MNPBurden_thresholded = case_when(
      !is.na(cumulative_mnp_ug_per_g) & cumulative_mnp_ug_per_g >= mnp_threshold_ug_g ~ "High",
      !is.na(cumulative_mnp_ug_per_g) & cumulative_mnp_ug_per_g <  mnp_threshold_ug_g ~ "Low",
      TRUE ~ NA_character_
    ),

    MNPBurden_final = coalesce(MNPBurden, MNPBurden_thresholded),

    MNPBurden_Delivery_class = case_when(
      !is.na(MNPBurden_final) & !is.na(Delivery_class) ~ paste(MNPBurden_final, Delivery_class, sep = "_"),
      TRUE ~ NA_character_
    ),

    high_mnp_binary = case_when(
      MNPBurden_final == "High" ~ 1,
      MNPBurden_final == "Low" ~ 0,
      TRUE ~ NA_real_
    ),

    mnp_threshold_ug_g = mnp_threshold_ug_g,

    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden_final, levels = c("Low", "High")),
    MNPBurden_Delivery_class = factor(
      MNPBurden_Delivery_class,
      levels = c("Low_Term", "High_Term", "Low_Preterm", "High_Preterm")
    ),
    platform = factor(platform, levels = c("VisiumV1", "Xenium", "VisiumHD")),
    platform_clean = factor(platform_clean, levels = c("VisiumV1", "Xenium", "VisiumHD")),

    GA_continuous = suppressWarnings(as.numeric(ga_continuous)),
    Maternal_age_years = suppressWarnings(as.numeric(maternal_age_years)),
    BMI = suppressWarnings(as.numeric(bmi)),
    SDI_score_from_zip = suppressWarnings(as.numeric(sdi_score_from_zip)),
    SVI_score = suppressWarnings(as.numeric(svi_score)),
    Prior_PTB = suppressWarnings(as.numeric(prior_ptb)),
    Py_GC_MS = as.logical(py_gc_ms),

    Microcalc_presence = suppressWarnings(as.numeric(microcalc_presence)),
    Massive_PVF = suppressWarnings(as.numeric(massive_pvf)),
    calcifications_number = suppressWarnings(as.numeric(calcifications_number)),
    pvf_percent_coverage = suppressWarnings(as.numeric(pvf_percent_coverage))
  )

# -----------------------------
# 4. Restrict new metadata to approved downstream columns
# -----------------------------
# These are the ONLY new columns allowed to overwrite/add.
# This avoids accidental deletion of RCTD/origin/QC metadata.
approved_update_cols <- c(
  # join / core sample identifiers
  "sample_id_std",
  "sample_id",
  "study_id",

  # platform fields
  "platform",
  "platform_csv",
  "platform_inferred",
  "platform_clean",

  # optional sample provenance from new table
  "xen_id_seurat",
  "vis_v1_seurat",
  "vis_hd_seurat",
  "xenium_id_seq",
  "visium_v1_seq",
  "visium_hd_seq",
  "xenium_platform",
  "visium_platform",
  "multiple_spatial_ro_is",
  "novel_spatial",
  "dataset",
  "xenium500_region",
  "xenium5k_run_id",
  "xenium5k_run_date",
  "xenium5k_slide_id",
  "xenium5k_block",
  "xenium5k_slide",
  "xenium5k_run",

  # clinical / delivery
  "GA_continuous",
  "Delivery_class",
  "Ethnicity_Race",
  "Preeclampsia_bin",
  "infant_gender",
  "birth_weight",
  "Prior_PTB",
  "Maternal_age_years",
  "BMI",
  "SDI_score_from_zip",
  "SVI_score",

  # exposure group
  "Py_GC_MS",
  "MNPBurden",
  "MNPBurden_thresholded",
  "MNPBurden_final",
  "MNPBurden_Delivery_class",
  "high_mnp_binary",
  "mnp_threshold_ug_g",

  # cumulative exposure
  "cumulative_mnp_ug_per_g",
  "cumulative_mnp_log10ugpg",

  # older readable polymer values, from current QC'd table
  "pe_ug_per_g",
  "sbr_ug_per_g",
  "pp_ug_per_g",
  "pvc_ug_per_g",
  "n66_ug_per_g",
  "pet_ug_per_g",
  "n6_ug_per_g",
  "pmma_ug_per_g",
  "pc_ug_per_g",
  "pu_ug_per_g",
  "abs_ug_per_g",
  "ps_ug_per_g",
  "pe_log10ugpg",
  "sbr_log10ugpg",
  "pp_log10ugpg",
  "pvc_log10ugpg",
  "n66_log10ugpg",
  "pet_log10ugpg",
  "n6_log10ugpg",
  "pmma_log10ugpg",
  "pc_log10ugpg",
  "pu_log10ugpg",
  "abs_log10ugpg",
  "ps_log10ugpg",

  # primary downstream polymer models
  "log2_mnp_abs_ug_g_plus1",
  "log2_mnp_n6_ug_g_plus1",
  "log2_mnp_n66_ug_g_plus1",
  "log2_mnp_pc_ug_g_plus1",
  "log2_mnp_pe_ug_g_plus1",
  "log2_mnp_pet_ug_g_plus1",
  "log2_mnp_pmma_ug_g_plus1",
  "log2_mnp_pp_ug_g_plus1",
  "log2_mnp_ps_ug_g_plus1",
  "log2_mnp_pu_ug_g_plus1",
  "log2_mnp_pvc_ug_g_plus1",
  "log2_mnp_sbr_ug_g_plus1",

  # histology
  "calcifications_number",
  "pvf_percent_coverage",
  "Microcalc_presence",
  "Massive_PVF"
)

meta_new_subset <- meta_new %>%
  select(any_of(approved_update_cols))

# Confirm uniqueness of join key.
if (anyDuplicated(meta_new_subset$sample_id_std)) {
  dupes <- meta_new_subset %>%
    count(sample_id_std) %>%
    filter(n > 1)

  print(dupes)
  stop("sample_id_std is duplicated in new metadata. Fix before joining.")
}

# -----------------------------
# 5. Audit columns that will be overwritten
# -----------------------------
cols_to_replace <- intersect(
  colnames(obj_meta1),
  setdiff(colnames(meta_new_subset), "sample_id_std")
)

cols_preserved <- setdiff(colnames(obj_meta1), cols_to_replace)

cat("\nColumns that WILL be overwritten/replaced:\n")
print(cols_to_replace)

cat("\nColumns that WILL be preserved:\n")
print(cols_preserved)

write_csv(
  tibble(column = cols_to_replace),
  file.path(outdir, "metadata_columns_overwritten.csv")
)

write_csv(
  tibble(column = cols_preserved),
  file.path(outdir, "metadata_columns_preserved.csv")
)

# -----------------------------
# 6. Join new metadata
# -----------------------------
obj_meta_clean <- obj_meta1 %>%
  select(-all_of(cols_to_replace))

obj_meta_joined <- obj_meta_clean %>%
  left_join(meta_new_subset, by = "sample_id_std")

# Row-level integrity checks.
stopifnot(nrow(obj_meta_joined) == nrow(obj_meta1))
stopifnot(all(obj_meta_joined$cell_barcode == obj_meta1$cell_barcode))

# -----------------------------
# 7. Optional compatibility aliases
# -----------------------------
# These help older plotting scripts while creating cleaner downstream names.
obj_meta_joined <- obj_meta_joined %>%
  mutate(
    # Keep legacy-style aliases populated from new QC'd columns.
    Cumulative_MNPs_ug_per_g = cumulative_mnp_ug_per_g,
    Cumulative_MNPs_log10ugpg = cumulative_mnp_log10ugpg,

    # Keep GA alias but downstream should use GA_continuous.
    GA = GA_continuous,

    # Legacy binary aliases if desired.
    Prior_PTB_TRUE = case_when(
      !is.na(Prior_PTB) & Prior_PTB == 1 ~ TRUE,
      !is.na(Prior_PTB) & Prior_PTB == 0 ~ FALSE,
      TRUE ~ NA
    )
  )

# -----------------------------
# 8. Restore object metadata
# -----------------------------
obj@meta.data <- obj_meta_joined %>%
  column_to_rownames("cell_barcode")

# -----------------------------
# 9. QA summaries
# -----------------------------
meta_cell <- obj@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  as_tibble()

# Cell/spot-level counts
qa_cell_counts <- meta_cell %>%
  count(platform_clean, Delivery_class, MNPBurden, name = "n_cells_or_spots") %>%
  arrange(platform_clean, Delivery_class, MNPBurden)

# Sample-level design
qa_sample_design <- meta_cell %>%
  distinct(sample_id_std, study_id, platform_clean, Delivery_class, MNPBurden,
           cumulative_mnp_ug_per_g, cumulative_mnp_log10ugpg) %>%
  count(platform_clean, Delivery_class, MNPBurden, name = "n_sample_rows") %>%
  arrange(platform_clean, Delivery_class, MNPBurden)

# Subject-level MNP analytical design
qa_subject_design <- meta_cell %>%
  distinct(study_id, Delivery_class, MNPBurden,
           cumulative_mnp_ug_per_g, cumulative_mnp_log10ugpg) %>%
  filter(!is.na(MNPBurden)) %>%
  count(Delivery_class, MNPBurden, name = "n_subjects") %>%
  arrange(Delivery_class, MNPBurden)

# Subject-level list
qa_subject_list <- meta_cell %>%
  distinct(study_id, Delivery_class, MNPBurden,
           cumulative_mnp_ug_per_g, cumulative_mnp_log10ugpg,
           GA_continuous, platform_clean) %>%
  arrange(study_id, platform_clean)

# Platform mismatch report
qa_platform_mismatch <- meta_new %>%
  filter(!is.na(platform_inferred), !is.na(platform_csv), platform_inferred != platform_csv) %>%
  select(sample_id_std, platform_csv, platform_inferred, platform_clean)

# Missingness at sample level
qa_missingness <- meta_cell %>%
  distinct(sample_id_std, study_id, platform_clean, Delivery_class, MNPBurden,
           cumulative_mnp_ug_per_g, cumulative_mnp_log10ugpg,
           GA_continuous, Maternal_age_years, BMI,
           SDI_score_from_zip, SVI_score,
           Microcalc_presence, Massive_PVF) %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "n_missing") %>%
  arrange(desc(n_missing), column)

cat("\nCell/spot-level counts:\n")
print(qa_cell_counts, n = Inf)

cat("\nSample-level design:\n")
print(qa_sample_design, n = Inf)

cat("\nSubject-level MNP analytical design:\n")
print(qa_subject_design, n = Inf)

cat("\nPlatform mismatch report:\n")
print(qa_platform_mismatch, n = Inf)

cat("\nSample-level missingness:\n")
print(qa_missingness, n = Inf)

# Save QA outputs.
write_csv(qa_cell_counts, file.path(outdir, "qa_cell_or_spot_counts_by_platform_delivery_mnp.csv"))
write_csv(qa_sample_design, file.path(outdir, "qa_sample_level_design.csv"))
write_csv(qa_subject_design, file.path(outdir, "qa_subject_level_mnp_design.csv"))
write_csv(qa_subject_list, file.path(outdir, "qa_subject_sample_list.csv"))
write_csv(qa_platform_mismatch, file.path(outdir, "qa_platform_mismatch_report.csv"))
write_csv(qa_missingness, file.path(outdir, "qa_sample_level_missingness.csv"))

# -----------------------------
# 10. Save updated object
# -----------------------------
saveRDS(obj, file.path(outdir, "spatial_object_metadata_updated_mnp_qc.rds"))

cat("\nDone. Updated object saved to:\n")


# explicit analysis flags
obj@meta.data <- obj@meta.data %>%
  dplyr::mutate(
    mnp_model_sample = !is.na(MNPBurden) & !is.na(cumulative_mnp_ug_per_g),
    mnp_model_subject = !is.na(MNPBurden) & !is.na(cumulative_mnp_ug_per_g),
    histology_model_sample = !is.na(Microcalc_presence) | !is.na(Massive_PVF)
  )

obj@meta.data %>%
  tibble::as_tibble() %>%
  dplyr::distinct(sample_id_std, study_id, platform_clean, Delivery_class,
                  MNPBurden, mnp_model_sample) %>%
  dplyr::count(mnp_model_sample, platform_clean, Delivery_class, MNPBurden) 
  





  obj_all <- obj 

  ## mnp-bearing samples
  obj_mnp <- subset(obj, subset = mnp_model_sample == TRUE)


obj_mnp@meta.data %>%
  tibble::as_tibble() %>%
  dplyr::distinct(sample_id_std, study_id, platform_clean, Delivery_class, MNPBurden,
                  cumulative_mnp_ug_per_g, cumulative_mnp_log10ugpg) %>%
  dplyr::arrange(Delivery_class, MNPBurden, study_id, platform_clean)  

obj_mnp@meta.data %>%
  tibble::as_tibble() %>%
  dplyr::distinct(study_id, Delivery_class, MNPBurden,
                  cumulative_mnp_ug_per_g, cumulative_mnp_log10ugpg) %>%
  dplyr::count(Delivery_class, MNPBurden)  




## module scores

library(Seurat)
library(dplyr)
library(tibble)
library(purrr)
library(broom)
library(sandwich)
library(lmtest)

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

# Check which genes are present
module_gene_check <- purrr::imap_dfr(module_list, function(genes, module) {
  tibble(
    module = module,
    gene = genes,
    present = genes %in% rownames(obj)
  )
})

print(module_gene_check %>% count(module, present))

# Keep only genes present in the object
module_list_present <- purrr::map(module_list, ~ intersect(.x, rownames(obj)))

module_list_present <- module_list_present[lengths(module_list_present) >= 2]

obj <- AddModuleScore(
  object = obj,
  features = module_list_present,
  name = names(module_list_present),
  assay = DefaultAssay(obj)
)

# AddModuleScore appends numbers; rename cleanly
new_score_cols <- tail(colnames(obj@meta.data), length(module_list_present))
clean_score_cols <- paste0(names(module_list_present), "_score")

colnames(obj@meta.data)[match(new_score_cols, colnames(obj@meta.data))] <- clean_score_cols

print(clean_score_cols)




## now remake MNP obj subset
obj_mnp <- subset(obj, subset = mnp_model_sample == TRUE)


# sample-level aggregation
module_score_cols <- c(
  "ECM_Remodeling_score",
  "Vascular_Interface_score",
  "Macrophage_Remodeling_score",
  "NK_Cytotoxic_Regulatory_score",
  "Trophoblast_Remodeling_score"
)

sample_scores <- obj_mnp@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  as_tibble() %>%
  group_by(
    sample_id_std,
    study_id,
    platform_clean,
    Delivery_class,
    MNPBurden,
    MNPBurden_Delivery_class,
    cumulative_mnp_ug_per_g,
    cumulative_mnp_log10ugpg,
    GA_continuous
  ) %>%
  summarise(
    n_cells_or_spots = n(),
    across(all_of(score_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

sample_scores %>%
  arrange(Delivery_class, MNPBurden, study_id, platform_clean) %>%
  print(n = Inf)

## subject-level aggregation, primary
subject_scores <- sample_scores %>%
  group_by(
    study_id,
    Delivery_class,
    MNPBurden,
    MNPBurden_Delivery_class,
    cumulative_mnp_ug_per_g,
    cumulative_mnp_log10ugpg,
    GA_continuous
  ) %>%
  summarise(
    n_samples = n(),
    platforms = paste(sort(unique(as.character(platform_clean))), collapse = ";"),
    across(all_of(score_cols), ~ weighted.mean(.x, w = n_cells_or_spots, na.rm = TRUE)),
    .groups = "drop"
  )

subject_scores %>%
  arrange(Delivery_class, MNPBurden, study_id) %>%
  print(n = Inf)



## primary screen
run_subject_model <- function(score_col, df = subject_scores) {
  f <- as.formula(paste0(score_col, " ~ MNPBurden * Delivery_class"))
  fit <- lm(f, data = df)
  broom::tidy(fit, conf.int = TRUE) %>%
    mutate(
      score = score_col,
      model = "subject_level_MNPBurden_x_Delivery"
    )
}

subject_model_results <- purrr::map_dfr(score_cols, run_subject_model)

subject_model_results %>%
  arrange(score, term) %>%
  print(n = Inf)


## sensitivity model

run_sample_model <- function(score_col, df = sample_scores) {
  f <- as.formula(paste0(score_col, " ~ MNPBurden * Delivery_class + platform_clean"))
  fit <- lm(f, data = df)

  broom::tidy(fit, conf.int = TRUE) %>%
    mutate(
      score = score_col,
      model = "sample_level_platform_adjusted"
    )
}

sample_model_results <- purrr::map_dfr(score_cols, run_sample_model)

sample_model_results %>%
  arrange(score, term) %>%
  print(n = Inf)

## cluster-robust SE by subject
run_sample_model_clustered <- function(score_col, df = sample_scores) {
  f <- as.formula(paste0(score_col, " ~ MNPBurden * Delivery_class + platform_clean"))
  fit <- lm(f, data = df)

  ct <- lmtest::coeftest(
    fit,
    vcov = sandwich::vcovCL(fit, cluster = df$study_id, type = "HC1")
  )

  broom::tidy(ct) %>%
    mutate(
      score = score_col,
      model = "sample_level_platform_adjusted_clustered_by_subject"
    )
}

sample_clustered_results <- purrr::map_dfr(score_cols, run_sample_model_clustered)

sample_clustered_results %>%
  arrange(score, term) %>%
  print(n = Inf)

## continuous MNP sensitivity
run_continuous_subject_model <- function(score_col, df = subject_scores) {
  f <- as.formula(paste0(score_col, " ~ cumulative_mnp_log10ugpg * Delivery_class"))
  fit <- lm(f, data = df)

  broom::tidy(fit, conf.int = TRUE) %>%
    mutate(
      score = score_col,
      model = "subject_level_continuous_MNP_x_Delivery"
    )
}

continuous_subject_results <- purrr::map_dfr(score_cols, run_continuous_subject_model)

continuous_subject_results %>%
  arrange(score, term) %>%
  print(n = Inf)

# save outputs
dir.create("module_score_outputs", showWarnings = FALSE)

write.csv(sample_scores, "module_score_outputs/sample_level_module_scores.csv", row.names = FALSE)
write.csv(subject_scores, "module_score_outputs/subject_level_module_scores.csv", row.names = FALSE)
write.csv(subject_model_results, "module_score_outputs/subject_level_mnpburden_interaction_models.csv", row.names = FALSE)
write.csv(sample_model_results, "module_score_outputs/sample_level_platform_adjusted_models.csv", row.names = FALSE)
write.csv(sample_clustered_results, "module_score_outputs/sample_level_platform_adjusted_clustered_models.csv", row.names = FALSE)
write.csv(continuous_subject_results, "module_score_outputs/subject_level_continuous_mnp_interaction_models.csv", row.names = FALSE)

#saveRDS(obj, "metadata_update_outputs/spatial_object_metadata_updated_with_module_scores.rds") 


# obj <- readRDS("/home/ebarrozo/visiumhd_mnp/results/audra/metadata_update_outputs/spatial_object_metadata_updated_with_module_scores.rds")


library(dplyr)
library(tibble)
library(purrr)
library(broom)
library(readr)

module_score_cols <- c(
  "ECM_Remodeling_score",
  "Vascular_Interface_score",
  "Macrophage_Remodeling_score",
  "NK_Cytotoxic_Regulatory_score",
  "Trophoblast_Remodeling_score"
)

sample_scores_z <- sample_scores %>%
  group_by(platform_clean) %>%
  mutate(
    across(
      all_of(module_score_cols),
      ~ as.numeric(scale(.x)),
      .names = "{.col}_platform_z"
    )
  ) %>%
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
    GA_continuous
  ) %>%
  summarise(
    n_samples = n(),
    platforms = paste(sort(unique(as.character(platform_clean))), collapse = ";"),
    across(all_of(module_z_cols), ~ weighted.mean(.x, w = n_cells_or_spots, na.rm = TRUE)),
    .groups = "drop"
  )

run_subject_model_z <- function(score_col, df = subject_scores_z) {
  f <- as.formula(paste0(score_col, " ~ MNPBurden * Delivery_class"))
  fit <- lm(f, data = df)

  broom::tidy(fit, conf.int = TRUE) %>%
    mutate(
      score = score_col,
      model = "subject_level_platform_z_MNPBurden_x_Delivery"
    )
}

subject_model_results_z <- purrr::map_dfr(module_z_cols, run_subject_model_z)

run_continuous_subject_model_z <- function(score_col, df = subject_scores_z) {
  f <- as.formula(paste0(score_col, " ~ cumulative_mnp_log10ugpg * Delivery_class"))
  fit <- lm(f, data = df)

  broom::tidy(fit, conf.int = TRUE) %>%
    mutate(
      score = score_col,
      model = "subject_level_platform_z_continuous_MNP_x_Delivery"
    )
}

continuous_subject_results_z <- purrr::map_dfr(module_z_cols, run_continuous_subject_model_z)

write_csv(sample_scores_z, "module_score_outputs/sample_level_module_scores_platform_z.csv")
write_csv(subject_scores_z, "module_score_outputs/subject_level_module_scores_platform_z.csv")
write_csv(subject_model_results_z, "module_score_outputs/subject_level_platform_z_mnpburden_interaction_models.csv")
write_csv(continuous_subject_results_z, "module_score_outputs/subject_level_platform_z_continuous_mnp_interaction_models.csv")










library(emmeans)
library(dplyr)
library(purrr)
library(readr)
library(tibble)

module_z_cols <- c(
  "ECM_Remodeling_score_platform_z",
  "Vascular_Interface_score_platform_z",
  "Macrophage_Remodeling_score_platform_z",
  "NK_Cytotoxic_Regulatory_score_platform_z",
  "Trophoblast_Remodeling_score_platform_z"
)

run_within_delivery_contrasts <- function(score_col, df = subject_scores_z) {
  f <- as.formula(paste0(score_col, " ~ MNPBurden * Delivery_class"))
  fit <- lm(f, data = df)

  emm <- emmeans(fit, ~ MNPBurden | Delivery_class)

  contrast(emm, method = "revpairwise") %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(
      score = score_col,
      model = "within_delivery_high_vs_low_subject_platform_z"
    )
}

within_delivery_contrasts <- purrr::map_dfr(
  module_z_cols,
  run_within_delivery_contrasts
) %>%
  mutate(
    q.value = p.adjust(p.value, method = "BH")
  ) %>%
  arrange(score, Delivery_class)

write_csv(
  within_delivery_contrasts,
  "module_score_outputs/within_delivery_high_vs_low_platform_z_contrasts.csv"
)

within_delivery_contrasts %>%
  print(n = Inf)


# ============================================================
# Plot: Subject-level platform-z module scores
# Purpose:
#   High vs low MNP burden within delivery class
#   Exact p-value annotation from emmeans contrasts
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggbeeswarm)
})

# -----------------------------
# Inputs
# -----------------------------
subject_file <- "module_score_outputs/subject_level_module_scores_platform_z.csv"
contrast_file <- "module_score_outputs/within_delivery_high_vs_low_platform_z_contrasts.csv"

outdir <- "module_score_outputs/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load data
# -----------------------------
subject_scores_z <- read_csv(subject_file, show_col_types = FALSE)
within_delivery_contrasts <- read_csv(contrast_file, show_col_types = FALSE)

# -----------------------------
# Define module columns and labels
# -----------------------------
module_z_cols <- c(
  "ECM_Remodeling_score_platform_z",
  "Vascular_Interface_score_platform_z",
  "Macrophage_Remodeling_score_platform_z",
  "NK_Cytotoxic_Regulatory_score_platform_z",
  "Trophoblast_Remodeling_score_platform_z"
)

module_labels <- c(
  ECM_Remodeling_score_platform_z = "ECM\nremodeling",
  Vascular_Interface_score_platform_z = "Vascular\ninterface",
  Macrophage_Remodeling_score_platform_z = "Macrophage\nremodeling",
  NK_Cytotoxic_Regulatory_score_platform_z = "NK cytotoxic /\nregulatory",
  Trophoblast_Remodeling_score_platform_z = "Trophoblast\nremodeling"
)

module_order <- names(module_labels)

# -----------------------------
# Long-format subject scores
# -----------------------------
plot_df <- subject_scores_z %>%
  select(
    study_id,
    Delivery_class,
    MNPBurden,
    n_samples,
    platforms,
    all_of(module_z_cols)
  ) %>%
  pivot_longer(
    cols = all_of(module_z_cols),
    names_to = "module",
    values_to = "module_score_platform_z"
  ) %>%
  mutate(
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High")),
    module = factor(module, levels = module_order, labels = module_labels[module_order])
  )

# -----------------------------
# P-value annotation table
# -----------------------------
p_annot <- within_delivery_contrasts %>%
  filter(score %in% module_z_cols) %>%
  mutate(
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    module = factor(score, levels = module_order, labels = module_labels[module_order]),
    p_label = paste0("p = ", signif(p.value, 2))
  ) %>%
  group_by(Delivery_class, module) %>%
  summarise(
    p_label = first(p_label),
    estimate = first(estimate),
    p.value = first(p.value),
    .groups = "drop"
  )

# Put p labels slightly above each module's local max in each facet
p_positions <- plot_df %>%
  group_by(Delivery_class, module) %>%
  summarise(
    y.position = max(module_score_platform_z, na.rm = TRUE) + 0.35,
    .groups = "drop"
  )

p_annot <- p_annot %>%
  left_join(p_positions, by = c("Delivery_class", "module"))

# Optional: only annotate nominal p < 0.10 to avoid clutter
p_annot_nominal <- p_annot %>%
  filter(p.value < 0.10)

# If you want every exact p-value, use p_annot instead of p_annot_nominal below.
annotation_df <- p_annot_nominal

# -----------------------------
# Plot
# -----------------------------
p_module <- ggplot(
  plot_df,
  aes(
    x = module,
    y = module_score_platform_z,
    color = MNPBurden,
    fill = MNPBurden
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    color = "grey75"
  ) +
  geom_boxplot(
    width = 0.52,
    alpha = 0.18,
    outlier.shape = NA,
    linewidth = 0.45,
    position = position_dodge(width = 0.68)
  ) +
  ggbeeswarm::geom_quasirandom(
    aes(shape = MNPBurden),
    dodge.width = 0.68,
    width = 0.08,
    size = 2.8,
    alpha = 0.92,
    stroke = 0.45
  ) +
  geom_text(
    data = annotation_df,
    aes(
      x = module,
      y = y.position,
      label = p_label
    ),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "plain",
    color = "black"
  ) +
  facet_wrap(
    ~ Delivery_class,
    nrow = 1
  ) +
  scale_color_manual(
    values = c(
      Low = "#64748B",
      High = "#D97706"
    )
  ) +
  scale_fill_manual(
    values = c(
      Low = "#64748B",
      High = "#D97706"
    )
  ) +
  scale_shape_manual(
    values = c(
      Low = 21,
      High = 24
    )
  ) +
  labs(
    x = NULL,
    y = "Platform-normalized module score",
    color = "MNP burden",
    fill = "MNP burden",
    shape = "MNP burden"
  ) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      color = "black",
      size = 9.5,
      lineheight = 0.9
    ),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(color = "black", size = 11),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12,
      color = "black"
    ),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    panel.spacing.x = unit(1.0, "lines"),
    plot.margin = margin(8, 10, 8, 10)
  )

p_module

# -----------------------------
# Save
# -----------------------------

ggsave(
  filename = file.path(outdir, "Fig_module_scores_platform_z_subject_level.png"),
  plot = p_module,
  width = 9,
  height = 3.8,
  units = "in",
  dpi = 600
)










# 
# All spatial samples
obj_all@meta.data %>%
  as_tibble() %>%
  distinct(sample_id_std, study_id, platform_clean, Delivery_class, MNPBurden, mnp_model_sample) %>%
  count(platform_clean, Delivery_class, MNPBurden, mnp_model_sample)

# MNP-bearing subset only
obj_mnp@meta.data %>%
  as_tibble() %>%
  distinct(sample_id_std, study_id, platform_clean, Delivery_class, MNPBurden) %>%
  count(platform_clean, Delivery_class, MNPBurden)

# Subject-level MNP subset
obj_mnp@meta.data %>%
  as_tibble() %>%
  distinct(study_id, Delivery_class, MNPBurden) %>%
  count(Delivery_class, MNPBurden)



# ============================================================
# Fig. 2 composition summaries
# Purpose:
#   1. Summarize predicted cell-type composition across platforms/subjects
#   2. Compare high vs low MNP burden within delivery class
#   3. Correlate cell-type composition with module scores
#
# Main annotation column:
#   nn.predicted.celltype
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
  library(broom)
})

# -----------------------------
# Inputs / outputs
# -----------------------------
outdir <- "fig2_composition_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

annotation_col <- "nn.predicted.celltype"

stopifnot(annotation_col %in% colnames(obj@meta.data))
stopifnot("sample_id_std" %in% colnames(obj@meta.data))
stopifnot("study_id" %in% colnames(obj@meta.data))
stopifnot("platform_clean" %in% colnames(obj@meta.data))
stopifnot("Delivery_class" %in% colnames(obj@meta.data))
stopifnot("MNPBurden" %in% colnames(obj@meta.data))
stopifnot("mnp_model_sample" %in% colnames(obj@meta.data))

# -----------------------------
# Pull metadata
# -----------------------------
meta_all <- obj@meta.data %>%
  rownames_to_column("cell_barcode") %>%
  as_tibble() %>%
  mutate(
    celltype_main = .data[[annotation_col]],
    celltype_main = as.character(celltype_main),
    celltype_main = if_else(is.na(celltype_main) | celltype_main == "", "Unassigned", celltype_main),
    platform_clean = factor(platform_clean, levels = c("VisiumV1", "Xenium", "VisiumHD")),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

# -----------------------------
# Basic annotation QA
# -----------------------------
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
  ungroup() %>%
  arrange(platform_clean, Delivery_class, MNPBurden, desc(prop))

write_csv(annotation_counts_mnp, file.path(outdir, "annotation_counts_mnp_by_platform_delivery_burden.csv"))

# -----------------------------
# Sample-level composition
# -----------------------------
sample_comp_long <- meta_all %>%
  group_by(
    sample_id_std,
    study_id,
    platform_clean,
    Delivery_class,
    MNPBurden,
    mnp_model_sample,
    celltype_main
  ) %>%
  summarise(
    n_celltype = n(),
    .groups = "drop"
  ) %>%
  group_by(sample_id_std) %>%
  mutate(
    n_total = sum(n_celltype),
    prop = n_celltype / n_total,
    pct = prop * 100
  ) %>%
  ungroup()

# Fill absent cell types with zeros per sample.
all_sample_keys <- meta_all %>%
  distinct(sample_id_std, study_id, platform_clean, Delivery_class, MNPBurden, mnp_model_sample)

all_celltypes <- sample_comp_long %>%
  distinct(celltype_main)

sample_comp_long <- all_sample_keys %>%
  tidyr::crossing(all_celltypes) %>%
  left_join(
    sample_comp_long,
    by = c(
      "sample_id_std", "study_id", "platform_clean",
      "Delivery_class", "MNPBurden", "mnp_model_sample",
      "celltype_main"
    )
  ) %>%
  group_by(sample_id_std) %>%
  mutate(
    n_celltype = replace_na(n_celltype, 0L),
    n_total = if_else(is.na(n_total), sum(n_celltype), n_total),
    prop = replace_na(prop, 0),
    pct = replace_na(pct, 0)
  ) %>%
  ungroup()

write_csv(sample_comp_long, file.path(outdir, "sample_level_celltype_composition_long.csv"))

sample_comp_wide <- sample_comp_long %>%
  select(sample_id_std, study_id, platform_clean, Delivery_class, MNPBurden, mnp_model_sample,
         celltype_main, prop) %>%
  mutate(celltype_main = make.names(celltype_main)) %>%
  pivot_wider(
    names_from = celltype_main,
    values_from = prop,
    names_prefix = "prop_",
    values_fill = 0
  )

write_csv(sample_comp_wide, file.path(outdir, "sample_level_celltype_composition.csv"))

# -----------------------------
# Subject-level composition
# -----------------------------
# Weighted by sample total so large spatial runs contribute proportionally.
# Also creates unweighted option for sensitivity.
subject_comp_long_weighted <- sample_comp_long %>%
  filter(mnp_model_sample == TRUE) %>%
  group_by(study_id, Delivery_class, MNPBurden, celltype_main) %>%
  summarise(
    n_celltype = sum(n_celltype),
    n_total = sum(n_total),
    prop = n_celltype / n_total,
    pct = prop * 100,
    n_samples = n_distinct(sample_id_std),
    platforms = paste(sort(unique(as.character(platform_clean))), collapse = ";"),
    .groups = "drop"
  )

subject_comp_long_unweighted <- sample_comp_long %>%
  filter(mnp_model_sample == TRUE) %>%
  group_by(study_id, Delivery_class, MNPBurden, celltype_main) %>%
  summarise(
    prop = mean(prop, na.rm = TRUE),
    pct = prop * 100,
    n_samples = n_distinct(sample_id_std),
    platforms = paste(sort(unique(as.character(platform_clean))), collapse = ";"),
    .groups = "drop"
  )

write_csv(
  subject_comp_long_weighted,
  file.path(outdir, "subject_level_celltype_composition_weighted_long.csv")
)

write_csv(
  subject_comp_long_unweighted,
  file.path(outdir, "subject_level_celltype_composition_unweighted_long.csv")
)

subject_comp_wide <- subject_comp_long_weighted %>%
  select(study_id, Delivery_class, MNPBurden, celltype_main, prop) %>%
  mutate(celltype_main = make.names(celltype_main)) %>%
  pivot_wider(
    names_from = celltype_main,
    values_from = prop,
    names_prefix = "prop_",
    values_fill = 0
  )

write_csv(subject_comp_wide, file.path(outdir, "subject_level_celltype_composition.csv"))

# -----------------------------
# Identify top cell types for plotting
# -----------------------------
top_celltypes <- subject_comp_long_weighted %>%
  group_by(celltype_main) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_prop)) %>%
  slice_head(n = 12) %>%
  pull(celltype_main)

cat("\nTop cell types in MNP subset:\n")
print(top_celltypes)

# -----------------------------
# Plot 1: all-platform composition by sample
# -----------------------------
plot_sample_comp <- sample_comp_long %>%
  filter(celltype_main %in% top_celltypes) %>%
  mutate(
    sample_label = paste0(sample_id_std, "\n", platform_clean),
    sample_label = fct_reorder(sample_label, as.numeric(platform_clean))
  )

p_platform <- ggplot(
  plot_sample_comp,
  aes(x = sample_label, y = pct, fill = celltype_main)
) +
  geom_col(width = 0.85) +
  facet_grid(
    mnp_model_sample ~ Delivery_class,
    scales = "free_x",
    space = "free_x",
    labeller = labeller(
      mnp_model_sample = c(`TRUE` = "Matched MNP subset", `FALSE` = "Spatial annotation / validation")
    )
  ) +
  labs(
    x = NULL,
    y = "Predicted cell-type composition (%)",
    fill = "Predicted cell type"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.8, "lines")
  )

ggsave(
  file.path(outdir, "Fig2_celltype_composition_by_platform.pdf"),
  p_platform,
  width = 11,
  height = 5.5,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "Fig2_celltype_composition_by_platform.png"),
  p_platform,
  width = 11,
  height = 5.5,
  dpi = 600
)

# -----------------------------
# Plot 2: MNP subset subject-level composition
# -----------------------------
plot_subject_comp <- subject_comp_long_weighted %>%
  filter(celltype_main %in% top_celltypes) %>%
  mutate(
    subject_label = paste0(study_id, "\n", Delivery_class, " / ", MNPBurden),
    subject_label = fct_reorder(subject_label, as.numeric(Delivery_class) * 10 + as.numeric(MNPBurden))
  )

p_subject <- ggplot(
  plot_subject_comp,
  aes(x = subject_label, y = pct, fill = celltype_main)
) +
  geom_col(width = 0.85) +
  facet_grid(
    ~ Delivery_class + MNPBurden,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    x = NULL,
    y = "Predicted cell-type composition (%)",
    fill = "Predicted cell type"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.spacing.x = unit(0.8, "lines")
  )

ggsave(
  file.path(outdir, "Fig2_mnp_subset_celltype_composition.pdf"),
  p_subject,
  width = 10.5,
  height = 4.8,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "Fig2_mnp_subset_celltype_composition.png"),
  p_subject,
  width = 10.5,
  height = 4.8,
  dpi = 600
)

# -----------------------------
# Within-delivery high vs low contrasts for cell-type proportions
# -----------------------------
# Use logit-transformed proportions with small offset.
# This avoids raw proportion boundary issues.
subject_comp_model <- subject_comp_long_weighted %>%
  mutate(
    prop_offset = (n_celltype + 0.5) / (n_total + 1),
    logit_prop = qlogis(prop_offset),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

run_celltype_contrast <- function(ct, df = subject_comp_model) {
  dat <- df %>%
    filter(celltype_main == ct)

  if (n_distinct(dat$MNPBurden) < 2 || n_distinct(dat$Delivery_class) < 2) {
    return(tibble())
  }

  fit <- lm(logit_prop ~ MNPBurden * Delivery_class, data = dat)

  emm <- emmeans(fit, ~ MNPBurden | Delivery_class)

  contrast(
    emm,
    method = list("High - Low" = c(-1, 1))
  ) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(
      celltype_main = ct,
      model = "subject_level_logit_prop_high_vs_low_within_delivery"
    )
}

celltype_composition_contrasts <- purrr::map_dfr(
  sort(unique(subject_comp_model$celltype_main)),
  run_celltype_contrast
) %>%
  group_by(Delivery_class) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(Delivery_class, p.value)

write_csv(
  celltype_composition_contrasts,
  file.path(outdir, "within_delivery_celltype_composition_contrasts.csv")
)

cat("\nTop cell-type composition contrasts:\n")
celltype_composition_contrasts %>%
  arrange(p.value) %>%
  print(n = 50)

# -----------------------------
# Composition-module correlations
# -----------------------------
# Requires subject_level_module_scores_platform_z.csv from prior step.
module_file <- "module_score_outputs/subject_level_module_scores_platform_z.csv"

if (file.exists(module_file)) {

  subject_scores_z <- read_csv(module_file, show_col_types = FALSE)

  module_z_cols <- c(
    "ECM_Remodeling_score_platform_z",
    "Vascular_Interface_score_platform_z",
    "Macrophage_Remodeling_score_platform_z",
    "NK_Cytotoxic_Regulatory_score_platform_z",
    "Trophoblast_Remodeling_score_platform_z"
  )

  module_z_cols <- intersect(module_z_cols, colnames(subject_scores_z))

  comp_for_corr <- subject_comp_wide %>%
    left_join(
      subject_scores_z %>%
        select(study_id, all_of(module_z_cols)),
      by = "study_id"
    )

  prop_cols <- grep("^prop_", colnames(comp_for_corr), value = TRUE)

  run_corr <- function(prop_col, module_col, df = comp_for_corr) {
    dat <- df %>%
      select(all_of(c(prop_col, module_col))) %>%
      drop_na()

    if (nrow(dat) < 5) {
      return(tibble())
    }

    ct <- suppressWarnings(cor.test(dat[[prop_col]], dat[[module_col]], method = "spearman", exact = FALSE))

    tibble(
      celltype_prop = prop_col,
      module_score = module_col,
      rho = unname(ct$estimate),
      p.value = ct$p.value,
      n = nrow(dat)
    )
  }

  module_celltype_correlations <- tidyr::crossing(
    prop_col = prop_cols,
    module_col = module_z_cols
  ) %>%
    mutate(result = map2(prop_col, module_col, run_corr)) %>%
    select(result) %>%
    unnest(result) %>%
    mutate(q.value = p.adjust(p.value, method = "BH")) %>%
    arrange(p.value)

  write_csv(
    module_celltype_correlations,
    file.path(outdir, "module_celltype_composition_correlations.csv")
  )

  cat("\nTop module-composition correlations:\n")
  module_celltype_correlations %>%
    arrange(p.value) %>%
    print(n = 50)
}

# -----------------------------
# Save RDS of key composition objects
# -----------------------------
saveRDS(
  list(
    sample_comp_long = sample_comp_long,
    sample_comp_wide = sample_comp_wide,
    subject_comp_long_weighted = subject_comp_long_weighted,
    subject_comp_long_unweighted = subject_comp_long_unweighted,
    subject_comp_wide = subject_comp_wide,
    celltype_composition_contrasts = celltype_composition_contrasts
  ),
  file.path(outdir, "fig2_composition_objects.rds")
)

cat("\nDone. Fig. 2 composition outputs written to:\n")
cat(outdir, "\n")














# ============================================================
# Fig. 2 composition plots with UMAP-matched IGV palette
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(scales)
  library(scCustomize)
  library(ggpubr)
})

outdir <- "fig2_composition_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# UMAP-matched palette
# -----------------------------
mypal3 <- get_palette("igv", 40)

cell_levels <- levels(factor(obj$nn.predicted.celltype))
cell_levels <- cell_levels[!is.na(cell_levels)]

col_map <- setNames(mypal3[seq_along(cell_levels)], cell_levels)

# Ensure all top cell types are covered
missing_cols <- setdiff(top_celltypes, names(col_map))
if (length(missing_cols) > 0) {
  warning("These top_celltypes are missing from col_map: ", paste(missing_cols, collapse = ", "))
}

# -----------------------------
# Plot 1: all samples/platforms
# -----------------------------
plot_sample_comp <- sample_comp_long %>%
  filter(celltype_main %in% top_celltypes) %>%
  mutate(
    celltype_main = factor(celltype_main, levels = cell_levels),
    sample_label = paste0(sample_id_std, "\n", platform_clean),
    sample_label = fct_reorder(sample_label, as.numeric(platform_clean))
  )

p_platform <- ggplot(
  plot_sample_comp,
  aes(x = sample_label, y = pct, fill = celltype_main)
) +
  geom_col(width = 0.85, color = NA) +
  facet_grid(
    mnp_model_sample ~ Delivery_class,
    scales = "free_x",
    space = "free_x",
    labeller = labeller(
      mnp_model_sample = c(
        `TRUE` = "Matched MNP subset",
        `FALSE` = "Spatial annotation / validation"
      )
    )
  ) +
  scale_fill_manual(
    values = col_map,
    breaks = cell_levels,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = NULL,
    y = "Predicted cell-type composition",
    fill = "Predicted cell type"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7,
      color = "black"
    ),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(color = "black", size = 10),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.8, "lines"),
    plot.margin = margin(6, 8, 6, 8)
  )

p_platform

ggsave(
  filename = file.path(outdir, "Fig2_celltype_composition_by_platform_UMAP_palette.png"),
  plot = p_platform,
  width = 11,
  height = 5.5,
  units = "in",
  dpi = 600
)

# -----------------------------
# Plot 2: MNP-bearing subject-level subset
# -----------------------------
plot_subject_comp <- subject_comp_long_weighted %>%
  filter(celltype_main %in% top_celltypes) %>%
  mutate(
    celltype_main = factor(celltype_main, levels = cell_levels),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High")),
    subject_label = paste0(study_id, "\n", Delivery_class, " / ", MNPBurden),
    subject_label = fct_reorder(
      subject_label,
      as.numeric(Delivery_class) * 10 + as.numeric(MNPBurden)
    )
  )

p_subject <- ggplot(
  plot_subject_comp,
  aes(x = subject_label, y = pct, fill = celltype_main)
) +
  geom_col(width = 0.85, color = NA) +
  facet_grid(
    ~ Delivery_class + MNPBurden,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = col_map,
    breaks = cell_levels,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = NULL,
    y = "Predicted cell-type composition",
    fill = "Predicted cell type"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8,
      color = "black"
    ),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(color = "black", size = 10),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.spacing.x = unit(1.0, "lines"),
    plot.margin = margin(6, 8, 6, 8)
  )

p_subject

ggsave(
  filename = file.path(outdir, "Fig2_mnp_subset_celltype_composition_UMAP_palette.png"),
  plot = p_subject,
  width = 10.5,
  height = 4.8,
  units = "in",
  dpi = 600
)


immune_celltypes <- c(
  "Macrophage",
  "Hofbauer Cell",
  "NK Cell",
  "T Cell",
  "B Cell",
  "Placental Endothelial",
  "Perivascular Cell"
)

plot_immune_comp <- subject_comp_long_weighted %>%
  filter(celltype_main %in% immune_celltypes) %>%
  mutate(
    celltype_main = factor(celltype_main, levels = immune_celltypes),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

p_immune <- ggplot(
  plot_immune_comp,
  aes(x = MNPBurden, y = pct, color = MNPBurden, fill = MNPBurden)
) +
  geom_boxplot(
    width = 0.48,
    alpha = 0.18,
    outlier.shape = NA,
    linewidth = 0.4
  ) +
  ggbeeswarm::geom_quasirandom(
    width = 0.08,
    size = 2.8,
    alpha = 0.95,
    stroke = 0.4
  ) +
  facet_grid(
    Delivery_class ~ celltype_main,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = c(Low = "#64748B", High = "#D97706")
  ) +
  scale_fill_manual(
    values = c(Low = "#64748B", High = "#D97706")
  ) +
  labs(
    x = NULL,
    y = "Predicted composition (%)",
    color = "MNP burden",
    fill = "MNP burden"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9, color = "black"),
    legend.position = "top",
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.8, "lines")
  )

p_immune

ggsave(
  filename = file.path(outdir, "Fig2_targeted_immune_celltype_composition.png"),
  plot = p_immune,
  width = 10.5,
  height = 4.8,
  units = "in",
  dpi = 600
)

celltype_composition_contrasts %>%
  arrange(p.value) %>%
  print(n = 30)






composition_celltypes_focus <- c(
  "Placental Endothelial",
  "Macrophage",
  "NK Cell",
  "ILC3",
  "Plasma Cell",
  "Villous Cytotrophoblast",
  "Extravillous Trophoblast"
)

plot_focus_comp <- subject_comp_long_weighted %>%
  filter(celltype_main %in% composition_celltypes_focus) %>%
  mutate(
    celltype_main = factor(celltype_main, levels = composition_celltypes_focus),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

# Pull exact p-values for selected cell types
p_comp_annot <- celltype_composition_contrasts %>%
  filter(celltype_main %in% composition_celltypes_focus) %>%
  mutate(
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    celltype_main = factor(celltype_main, levels = composition_celltypes_focus),
    p_label = paste0("p = ", signif(p.value, 2))
  ) %>%
  filter(p.value < 0.10)

p_pos <- plot_focus_comp %>%
  group_by(Delivery_class, celltype_main) %>%
  summarise(
    y_max = max(pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y.position = pmax(y_max * 1.12 + 0.3, 0.9)
  )

p_comp_annot <- p_comp_annot %>%
  left_join(p_pos, by = c("Delivery_class", "celltype_main"))

p_focus_comp <- ggplot(
  plot_focus_comp,
  aes(x = MNPBurden, y = pct, color = MNPBurden, fill = MNPBurden)
) +
  geom_boxplot(
    width = 0.48,
    alpha = 0.18,
    outlier.shape = NA,
    linewidth = 0.4
  ) +
  ggbeeswarm::geom_quasirandom(
    width = 0.08,
    size = 2.8,
    alpha = 0.95,
    stroke = 0.4
  ) +
  geom_text(
    data = p_comp_annot,
    aes(
      x = 1.5,
      y = y.position,
      label = p_label
    ),
    inherit.aes = FALSE,
    size = 3.1,
    color = "black",
    hjust = 0.5,
    vjust = -0.9
  ) +
  facet_grid(
    Delivery_class ~ celltype_main,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = c(Low = "#64748B", High = "#D97706")
  ) +
  scale_fill_manual(
    values = c(Low = "#64748B", High = "#D97706")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.18))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Predicted composition (%)",
    color = "MNP burden",
    fill = "MNP burden"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 10),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9, color = "black"),
    legend.position = "top",
    panel.spacing.x = unit(0.75, "lines"),
    panel.spacing.y = unit(0.75, "lines"),
    plot.margin = margin(10, 8, 8, 8)
  )

p_focus_comp

ggsave(
  filename = file.path(outdir, "Fig2_focused_celltype_composition_MNP_contrasts.png"),
  plot = p_focus_comp,
  width = 12,
  height = 4.8,
  units = "in",
  dpi = 600
)












# ============================================================
# Polymer-specific module models
# Purpose:
#   Test whether polymer-specific MNP burden associates with
#   platform-normalized module scores within delivery class.
#
# Uses:
#   subject_level_module_scores_platform_z.csv
#
# Outputs:
#   polymer_module_outputs/
#     polymer_module_subject_model_terms.csv
#     polymer_module_within_delivery_slopes.csv
#     polymer_module_within_delivery_slopes_heatmap.png
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(broom)
  library(emmeans)
  library(purrr)
  library(ggplot2)
})


## remake z scores with polymer-specific columns

polymer_cols <- grep(
  "^log2_mnp_.*_ug_g_plus1$",
  colnames(obj_mnp@meta.data),
  value = TRUE
)

print(polymer_cols)


module_score_cols <- c(
  "ECM_Remodeling_score",
  "Vascular_Interface_score",
  "Macrophage_Remodeling_score",
  "NK_Cytotoxic_Regulatory_score",
  "Trophoblast_Remodeling_score"
)

module_score_cols <- intersect(module_score_cols, colnames(obj_mnp@meta.data))

polymer_cols <- grep(
  "^log2_mnp_.*_ug_g_plus1$",
  colnames(obj_mnp@meta.data),
  value = TRUE
)

stopifnot(length(module_score_cols) > 0)
stopifnot(length(polymer_cols) > 0)

sample_scores <- obj_mnp@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  tibble::as_tibble() %>%
  group_by(
    sample_id_std,
    study_id,
    platform_clean,
    Delivery_class,
    MNPBurden,
    MNPBurden_Delivery_class,
    cumulative_mnp_ug_per_g,
    cumulative_mnp_log10ugpg,
    GA_continuous
  ) %>%
  summarise(
    n_cells_or_spots = n(),
    across(all_of(module_score_cols), ~ mean(.x, na.rm = TRUE)),
    across(all_of(polymer_cols), ~ first(.x)),
    .groups = "drop"
  )

sample_scores_z <- sample_scores %>%
  group_by(platform_clean) %>%
  mutate(
    across(
      all_of(module_score_cols),
      ~ as.numeric(scale(.x)),
      .names = "{.col}_platform_z"
    )
  ) %>%
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
    GA_continuous
  ) %>%
  summarise(
    n_samples = n(),
    platforms = paste(sort(unique(as.character(platform_clean))), collapse = ";"),
    across(all_of(module_z_cols), ~ weighted.mean(.x, w = n_cells_or_spots, na.rm = TRUE)),
    across(all_of(polymer_cols), ~ first(.x)),
    .groups = "drop"
  )

readr::write_csv(
  subject_scores_z,
  "module_score_outputs/subject_level_module_scores_platform_z_with_polymers.csv"
)

readr::write_csv(
  sample_scores_z,
  "module_score_outputs/sample_level_module_scores_platform_z_with_polymers.csv"
)


# -----------------------------
# Inputs / outputs
# -----------------------------
subject_file <- "module_score_outputs/subject_level_module_scores_platform_z_with_polymers.csv"

outdir <- "polymer_module_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

subject_scores_z <- read_csv(subject_file, show_col_types = FALSE)



# -----------------------------
# Define module and polymer columns
# -----------------------------
module_z_cols <- c(
  "ECM_Remodeling_score_platform_z",
  "Vascular_Interface_score_platform_z",
  "Macrophage_Remodeling_score_platform_z",
  "NK_Cytotoxic_Regulatory_score_platform_z",
  "Trophoblast_Remodeling_score_platform_z"
)

module_z_cols <- intersect(module_z_cols, colnames(subject_scores_z))

polymer_cols <- c(
  "log2_mnp_pet_ug_g_plus1",
  "log2_mnp_pu_ug_g_plus1",
  "log2_mnp_pvc_ug_g_plus1",
  "log2_mnp_sbr_ug_g_plus1",
  "log2_mnp_abs_ug_g_plus1",
  "log2_mnp_pc_ug_g_plus1",
  "log2_mnp_pp_ug_g_plus1",
  "log2_mnp_pe_ug_g_plus1",
  "log2_mnp_n6_ug_g_plus1",
  "log2_mnp_n66_ug_g_plus1",
  "log2_mnp_pmma_ug_g_plus1",
  "log2_mnp_ps_ug_g_plus1"
)

polymer_cols <- intersect(polymer_cols, colnames(subject_scores_z))

polymer_labels <- c(
  log2_mnp_pet_ug_g_plus1  = "PET",
  log2_mnp_pu_ug_g_plus1   = "PU",
  log2_mnp_pvc_ug_g_plus1  = "PVC",
  log2_mnp_sbr_ug_g_plus1  = "SBR",
  log2_mnp_abs_ug_g_plus1  = "ABS",
  log2_mnp_pc_ug_g_plus1   = "PC",
  log2_mnp_pp_ug_g_plus1   = "PP",
  log2_mnp_pe_ug_g_plus1   = "PE",
  log2_mnp_n6_ug_g_plus1   = "N6",
  log2_mnp_n66_ug_g_plus1  = "N66",
  log2_mnp_pmma_ug_g_plus1 = "PMMA",
  log2_mnp_ps_ug_g_plus1   = "PS"
)

module_labels <- c(
  ECM_Remodeling_score_platform_z = "ECM remodeling",
  Vascular_Interface_score_platform_z = "Vascular interface",
  Macrophage_Remodeling_score_platform_z = "Macrophage remodeling",
  NK_Cytotoxic_Regulatory_score_platform_z = "NK cytotoxic/regulatory",
  Trophoblast_Remodeling_score_platform_z = "Trophoblast remodeling"
)

# -----------------------------
# Data prep
# -----------------------------
model_df <- subject_scores_z %>%
  mutate(
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

# Quick missingness check
polymer_missingness <- model_df %>%
  summarise(across(all_of(polymer_cols), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "polymer_col", values_to = "n_missing") %>%
  mutate(polymer = polymer_labels[polymer_col]) %>%
  arrange(desc(n_missing), polymer)

print(polymer_missingness, n = Inf)

write_csv(
  polymer_missingness,
  file.path(outdir, "polymer_missingness_subject_level.csv")
)

# -----------------------------
# Fit model terms
# -----------------------------
run_polymer_module_model <- function(module_col, polymer_col, df = model_df) {

  dat <- df %>%
    select(
      study_id,
      Delivery_class,
      all_of(module_col),
      all_of(polymer_col)
    ) %>%
    drop_na()

  # Need enough observations and variation.
  if (nrow(dat) < 6 ||
      n_distinct(dat$Delivery_class) < 2 ||
      sd(dat[[polymer_col]], na.rm = TRUE) == 0) {
    return(tibble())
  }

  f <- as.formula(
    paste0(module_col, " ~ ", polymer_col, " * Delivery_class")
  )

  fit <- lm(f, data = dat)

  broom::tidy(fit, conf.int = TRUE) %>%
    mutate(
      module_col = module_col,
      module = module_labels[module_col],
      polymer_col = polymer_col,
      polymer = polymer_labels[polymer_col],
      n = nrow(dat),
      model = "subject_level_polymer_x_delivery"
    )
}

polymer_module_model_terms <- crossing(
  module_col = module_z_cols,
  polymer_col = polymer_cols
) %>%
  mutate(result = map2(module_col, polymer_col, run_polymer_module_model)) %>%
  select(result) %>%
  unnest(result) %>%
  mutate(
    q.value = p.adjust(p.value, method = "BH")
  ) %>%
  arrange(p.value)

write_csv(
  polymer_module_model_terms,
  file.path(outdir, "polymer_module_subject_model_terms.csv")
)

cat("\nTop polymer × module model terms:\n")
polymer_module_model_terms %>%
  arrange(p.value) %>%
  print(n = 50)

# -----------------------------
# Within-delivery polymer slopes
# -----------------------------
# This extracts polymer slope separately in Term and Preterm.
run_polymer_within_delivery_slopes <- function(module_col, polymer_col, df = model_df) {

  dat <- df %>%
    select(
      study_id,
      Delivery_class,
      all_of(module_col),
      all_of(polymer_col)
    ) %>%
    drop_na()

  if (nrow(dat) < 6 ||
      n_distinct(dat$Delivery_class) < 2 ||
      sd(dat[[polymer_col]], na.rm = TRUE) == 0) {
    return(tibble())
  }

  f <- as.formula(
    paste0(module_col, " ~ ", polymer_col, " * Delivery_class")
  )

  fit <- lm(f, data = dat)

  emtrends(
    fit,
    specs = ~ Delivery_class,
    var = polymer_col
  ) %>%
    as.data.frame() %>%
    as_tibble() %>%
    rename(
      slope = !!paste0(polymer_col, ".trend"),
      SE = SE,
      df = df,
      lower.CL = lower.CL,
      upper.CL = upper.CL
    ) %>%
    mutate(
      module_col = module_col,
      module = module_labels[module_col],
      polymer_col = polymer_col,
      polymer = polymer_labels[polymer_col],
      n = nrow(dat),
      model = "within_delivery_polymer_slope_subject_level"
    )
}

polymer_module_slopes <- crossing(
  module_col = module_z_cols,
  polymer_col = polymer_cols
) %>%
  mutate(result = map2(module_col, polymer_col, run_polymer_within_delivery_slopes)) %>%
  select(result) %>%
  unnest(result) %>%
  mutate(
    p.value = 2 * pt(abs(slope / SE), df = df, lower.tail = FALSE),
    q.value = p.adjust(p.value, method = "BH"),
    direction = case_when(
      slope > 0 ~ "Positive",
      slope < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    nominal = p.value < 0.05
  ) %>%
  arrange(p.value)

write_csv(
  polymer_module_slopes,
  file.path(outdir, "polymer_module_within_delivery_slopes.csv")
)

cat("\nTop within-delivery polymer slopes:\n")
polymer_module_slopes %>%
  arrange(p.value) %>%
  print(n = 60)

# -----------------------------
# Heatmap for within-delivery slopes
# -----------------------------
polymer_order <- c("PET", "PU", "PVC", "SBR", "ABS", "PC", "PP", "PE", "N6", "N66", "PMMA", "PS")
module_order <- c(
  "ECM remodeling",
  "Vascular interface",
  "Macrophage remodeling",
  "NK cytotoxic/regulatory",
  "Trophoblast remodeling"
)

heatmap_df <- polymer_module_slopes %>%
  mutate(
    polymer = factor(polymer, levels = polymer_order),
    module = factor(module, levels = module_order),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    sig_label = case_when(
      p.value < 0.05 ~ "●",
      TRUE ~ ""
    )
  )

p_poly_heatmap <- ggplot(
  heatmap_df,
  aes(x = polymer, y = module, fill = slope)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = sig_label),
    size = 3.2,
    color = "black"
  ) +
  facet_wrap(~ Delivery_class, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = "Slope"
  ) +
  labs(
    x = "Polymer burden",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 9),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11, color = "black"),
    legend.position = "right"
  )

p_poly_heatmap

ggsave(
  filename = file.path(outdir, "polymer_module_within_delivery_slopes_heatmap.png"),
  plot = p_poly_heatmap,
  width = 9.5,
  height = 3.8,
  units = "in",
  dpi = 600
)

# -----------------------------
# Focused polymer-module plot
# -----------------------------
# Optional: show key polymer relationships with subject points.
# This is useful for QC, not necessarily main figure yet.
focus_polymers <- c(
  "log2_mnp_pet_ug_g_plus1",
  "log2_mnp_pu_ug_g_plus1",
  "log2_mnp_pvc_ug_g_plus1",
  "log2_mnp_sbr_ug_g_plus1",
  "log2_mnp_abs_ug_g_plus1"
)

focus_modules <- c(
  "Macrophage_Remodeling_score_platform_z",
  "NK_Cytotoxic_Regulatory_score_platform_z",
  "Vascular_Interface_score_platform_z"
)

focus_long <- model_df %>%
  select(study_id, Delivery_class, MNPBurden, all_of(focus_polymers), all_of(focus_modules)) %>%
  pivot_longer(
    cols = all_of(focus_polymers),
    names_to = "polymer_col",
    values_to = "polymer_value"
  ) %>%
  pivot_longer(
    cols = all_of(focus_modules),
    names_to = "module_col",
    values_to = "module_score"
  ) %>%
  mutate(
    polymer = factor(polymer_labels[polymer_col],
                     levels = c("PET", "PU", "PVC", "SBR", "ABS")),
    module = factor(module_labels[module_col],
                    levels = c("Macrophage remodeling", "NK cytotoxic/regulatory", "Vascular interface")),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm"))
  )

p_poly_scatter <- ggplot(
  focus_long,
  aes(x = polymer_value, y = module_score, color = Delivery_class)
) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  facet_grid(module ~ polymer, scales = "free_x") +
  scale_color_manual(values = c(Term = "#0F766E", Preterm = "#7C3AED")) +
  labs(
    x = "Polymer burden, log2(µg/g + 1)",
    y = "Platform-normalized module score",
    color = "Delivery class"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title = element_text(color = "black", size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8.5, color = "black"),
    legend.position = "top",
    panel.spacing.x = unit(0.65, "lines"),
    panel.spacing.y = unit(0.65, "lines"),
    plot.margin = margin(6, 8, 6, 8)
  )

p_poly_scatter

ggsave(
  filename = file.path(outdir, "polymer_module_focused_scatter_QC.png"),
  plot = p_poly_scatter,
  width = 10.5,
  height = 5.2,
  units = "in",
  dpi = 600
)

cat("\nDone. Polymer-specific module outputs written to:\n")
cat(outdir, "\n")









# ============================================================
# Leave-one-subject-out robustness
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(emmeans)
  library(broom)
  library(readr)
})

outdir <- "robustness_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Inputs
# -----------------------------
subject_scores_z <- read_csv(
  "module_score_outputs/subject_level_module_scores_platform_z_with_polymers.csv",
  show_col_types = FALSE
) %>%
  mutate(
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

primary_modules <- c(
  "Macrophage_Remodeling_score_platform_z",
  "NK_Cytotoxic_Regulatory_score_platform_z"
)

polymer_focus <- c(
  "log2_mnp_pc_ug_g_plus1",
  "log2_mnp_pmma_ug_g_plus1"
)

polymer_module_focus <- c(
  "Trophoblast_Remodeling_score_platform_z",
  "Vascular_Interface_score_platform_z",
  "Macrophage_Remodeling_score_platform_z",
  "NK_Cytotoxic_Regulatory_score_platform_z"
)

module_labels <- c(
  Macrophage_Remodeling_score_platform_z = "Macrophage remodeling",
  NK_Cytotoxic_Regulatory_score_platform_z = "NK cytotoxic/regulatory",
  Trophoblast_Remodeling_score_platform_z = "Trophoblast remodeling",
  Vascular_Interface_score_platform_z = "Vascular interface"
)

polymer_labels <- c(
  log2_mnp_pc_ug_g_plus1 = "PC",
  log2_mnp_pmma_ug_g_plus1 = "PMMA"
)

# -----------------------------
# A. LOO for primary high-vs-low MNP module findings
# -----------------------------
run_loo_mnp_contrast <- function(drop_subject, module_col, df = subject_scores_z) {

  dat <- df %>%
    filter(study_id != drop_subject)

  # Need each delivery/burden cell represented after dropping
  design_check <- dat %>%
    count(Delivery_class, MNPBurden)

  if (nrow(design_check) < 4 || any(design_check$n < 1)) {
    return(tibble(
      dropped_subject = drop_subject,
      module_col = module_col,
      module = module_labels[module_col],
      Delivery_class = NA_character_,
      estimate = NA_real_,
      SE = NA_real_,
      df = NA_real_,
      p.value = NA_real_,
      note = "insufficient design after drop"
    ))
  }

  f <- as.formula(paste0(module_col, " ~ MNPBurden * Delivery_class"))
  fit <- lm(f, data = dat)

  emm <- emmeans(fit, ~ MNPBurden | Delivery_class)

  contrast(
    emm,
    method = list("High - Low" = c(-1, 1))
  ) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(
      dropped_subject = drop_subject,
      module_col = module_col,
      module = module_labels[module_col],
      note = "ok"
    )
}

loo_mnp_contrasts <- crossing(
  dropped_subject = unique(subject_scores_z$study_id),
  module_col = primary_modules
) %>%
  mutate(result = map2(dropped_subject, module_col, run_loo_mnp_contrast)) %>%
  select(result) %>%
  unnest(result) %>%
  arrange(module, Delivery_class, dropped_subject)

write_csv(
  loo_mnp_contrasts,
  file.path(outdir, "loo_primary_mnpburden_within_delivery_contrasts.csv")
)

# Summary for term primary findings
loo_mnp_summary <- loo_mnp_contrasts %>%
  filter(Delivery_class == "Term", note == "ok") %>%
  group_by(module) %>%
  summarise(
    n_runs = n(),
    min_estimate = min(estimate, na.rm = TRUE),
    max_estimate = max(estimate, na.rm = TRUE),
    median_estimate = median(estimate, na.rm = TRUE),
    n_positive = sum(estimate > 0, na.rm = TRUE),
    n_nominal_p05 = sum(p.value < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  loo_mnp_summary,
  file.path(outdir, "loo_primary_mnpburden_summary.csv")
)

print(loo_mnp_summary)

# -----------------------------
# B. LOO for polymer slopes
# -----------------------------
run_loo_polymer_slope <- function(drop_subject, module_col, polymer_col, df = subject_scores_z) {

  dat <- df %>%
    filter(study_id != drop_subject) %>%
    select(study_id, Delivery_class, all_of(module_col), all_of(polymer_col)) %>%
    drop_na()

  if (nrow(dat) < 6 ||
      n_distinct(dat$Delivery_class) < 2 ||
      sd(dat[[polymer_col]], na.rm = TRUE) == 0) {
    return(tibble(
      dropped_subject = drop_subject,
      module_col = module_col,
      module = module_labels[module_col],
      polymer_col = polymer_col,
      polymer = polymer_labels[polymer_col],
      Delivery_class = NA_character_,
      slope = NA_real_,
      SE = NA_real_,
      df = NA_real_,
      p.value = NA_real_,
      note = "insufficient design after drop"
    ))
  }

  f <- as.formula(paste0(module_col, " ~ ", polymer_col, " * Delivery_class"))
  fit <- lm(f, data = dat)

  emtrends(
    fit,
    specs = ~ Delivery_class,
    var = polymer_col
  ) %>%
    as.data.frame() %>%
    as_tibble() %>%
    rename(slope = !!paste0(polymer_col, ".trend")) %>%
    mutate(
      p.value = 2 * pt(abs(slope / SE), df = df, lower.tail = FALSE),
      dropped_subject = drop_subject,
      module_col = module_col,
      module = module_labels[module_col],
      polymer_col = polymer_col,
      polymer = polymer_labels[polymer_col],
      note = "ok"
    )
}

loo_polymer_slopes <- crossing(
  dropped_subject = unique(subject_scores_z$study_id),
  module_col = polymer_module_focus,
  polymer_col = polymer_focus
) %>%
  mutate(result = pmap(
    list(dropped_subject, module_col, polymer_col),
    run_loo_polymer_slope
  )) %>%
  select(result) %>%
  unnest(result) %>%
  arrange(polymer, module, Delivery_class, dropped_subject)

write_csv(
  loo_polymer_slopes,
  file.path(outdir, "loo_polymer_module_slopes.csv")
)

loo_polymer_summary <- loo_polymer_slopes %>%
  filter(note == "ok") %>%
  group_by(polymer, module, Delivery_class) %>%
  summarise(
    n_runs = n(),
    min_slope = min(slope, na.rm = TRUE),
    max_slope = max(slope, na.rm = TRUE),
    median_slope = median(slope, na.rm = TRUE),
    n_positive = sum(slope > 0, na.rm = TRUE),
    n_negative = sum(slope < 0, na.rm = TRUE),
    n_nominal_p05 = sum(p.value < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  loo_polymer_summary,
  file.path(outdir, "loo_polymer_module_slopes_summary.csv")
)

print(loo_polymer_summary, n = Inf)

# ============================================================
# Optional Fig. 3e robustness inset: LOO estimates
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
})

loo_file <- "robustness_outputs/loo_primary_mnpburden_within_delivery_contrasts.csv"
outdir <- "robustness_outputs/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

loo_df <- read_csv(loo_file, show_col_types = FALSE) %>%
  filter(
    Delivery_class == "Term",
    note == "ok",
    module %in% c("Macrophage remodeling", "NK cytotoxic/regulatory")
  ) %>%
  mutate(
    module = factor(
      module,
      levels = c("Macrophage remodeling", "NK cytotoxic/regulatory")
    )
  )

loo_summary <- loo_df %>%
  group_by(module) %>%
  summarise(
    min_estimate = min(estimate, na.rm = TRUE),
    max_estimate = max(estimate, na.rm = TRUE),
    median_estimate = median(estimate, na.rm = TRUE),
    n_positive = sum(estimate > 0, na.rm = TRUE),
    n_runs = n(),
    label = paste0(n_positive, "/", n_runs, " positive"),
    .groups = "drop"
  )

p_loo <- ggplot(loo_df, aes(x = estimate, y = module)) +
  geom_vline(
    xintercept = 0,
    color = "grey65",
    linewidth = 0.35
  ) +
  geom_segment(
    data = loo_summary,
    aes(
      x = min_estimate,
      xend = max_estimate,
      y = module,
      yend = module
    ),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "grey45"
  ) +
  geom_point(
    color = "#D97706",
    size = 2.8,
    alpha = 0.9
  ) +
  geom_point(
    data = loo_summary,
    aes(x = median_estimate, y = module),
    inherit.aes = FALSE,
    shape = 23,
    size = 3.6,
    fill = "white",
    color = "black",
    stroke = 0.7
  ) +
  geom_text(
    data = loo_summary,
    aes(
      x = max_estimate + 0.10,
      y = module,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.1
  ) +
  labs(
    x = "LOO term High − Low estimate",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 7) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    plot.margin = margin(6, 34, 6, 6)
  )

p_loo

ggsave(
  filename = file.path(outdir, "Fig3_optional_LOO_robustness_inset.png"),
  plot = p_loo,
  width = 4.2,
  height = 1.45,
  units = "in",
  dpi = 600
)



suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

outdir <- "robustness_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

subject_scores_z <- read_csv(
  "module_score_outputs/subject_level_module_scores_platform_z_with_polymers.csv",
  show_col_types = FALSE
) %>%
  mutate(
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

primary_modules <- c(
  "Macrophage_Remodeling_score_platform_z",
  "NK_Cytotoxic_Regulatory_score_platform_z"
)

module_labels <- c(
  Macrophage_Remodeling_score_platform_z = "Macrophage remodeling",
  NK_Cytotoxic_Regulatory_score_platform_z = "NK cytotoxic/regulatory"
)

# Direct within-term High - Low difference.
get_term_high_low_direct <- function(df, module_col, exposure_col = "MNPBurden") {
  df %>%
    filter(Delivery_class == "Term") %>%
    group_by(.data[[exposure_col]]) %>%
    summarise(mean_score = mean(.data[[module_col]], na.rm = TRUE), .groups = "drop") %>%
    summarise(
      estimate = mean_score[.data[[exposure_col]] == "High"] -
        mean_score[.data[[exposure_col]] == "Low"]
    ) %>%
    pull(estimate)
}

# Generate all unique label permutations within delivery class.
make_exact_perms <- function(df) {
  term_df <- df %>% filter(Delivery_class == "Term")
  preterm_df <- df %>% filter(Delivery_class == "Preterm")

  term_ids <- term_df$study_id
  preterm_ids <- preterm_df$study_id

  term_high_sets <- combn(term_ids, sum(term_df$MNPBurden == "High"), simplify = FALSE)
  preterm_high_sets <- combn(preterm_ids, sum(preterm_df$MNPBurden == "High"), simplify = FALSE)

  tidyr::crossing(
    term_idx = seq_along(term_high_sets),
    preterm_idx = seq_along(preterm_high_sets)
  ) %>%
    mutate(
      term_high = map(term_idx, ~ term_high_sets[[.x]]),
      preterm_high = map(preterm_idx, ~ preterm_high_sets[[.x]])
    )
}

perm_grid <- make_exact_perms(subject_scores_z)

cat("Number of exact within-delivery permutations:", nrow(perm_grid), "\n")
# Should be 60.

exact_perm_results <- map_dfr(primary_modules, function(module_col) {

  obs <- get_term_high_low_direct(subject_scores_z, module_col)

  perm_est <- map2_dbl(
    perm_grid$term_high,
    perm_grid$preterm_high,
    function(term_high, preterm_high) {

      df_perm <- subject_scores_z %>%
        mutate(
          MNPBurden_perm = case_when(
            Delivery_class == "Term" & study_id %in% term_high ~ "High",
            Delivery_class == "Term" & !study_id %in% term_high ~ "Low",
            Delivery_class == "Preterm" & study_id %in% preterm_high ~ "High",
            Delivery_class == "Preterm" & !study_id %in% preterm_high ~ "Low",
            TRUE ~ NA_character_
          ),
          MNPBurden_perm = factor(MNPBurden_perm, levels = c("Low", "High"))
        )

      get_term_high_low_direct(df_perm, module_col, exposure_col = "MNPBurden_perm")
    }
  )

  tibble(
    module_col = module_col,
    module = module_labels[module_col],
    observed_estimate = obs,
    perm_estimate = perm_est,
    perm_p_two_sided = mean(abs(perm_est) >= abs(obs)),
    perm_p_one_sided_positive = mean(perm_est >= obs)
  )
})

exact_perm_summary <- exact_perm_results %>%
  group_by(module_col, module, observed_estimate) %>%
  summarise(
    n_exact_permutations = n(),
    perm_p_two_sided = first(perm_p_two_sided),
    perm_p_one_sided_positive = first(perm_p_one_sided_positive),
    perm_mean = mean(perm_estimate),
    perm_sd = sd(perm_estimate),
    perm_q025 = quantile(perm_estimate, 0.025),
    perm_q975 = quantile(perm_estimate, 0.975),
    .groups = "drop"
  )

print(exact_perm_summary)

write_csv(
  exact_perm_results,
  file.path(outdir, "exact_permutation_primary_mnpburden_module_results_long.csv")
)

write_csv(
  exact_perm_summary,
  file.path(outdir, "exact_permutation_primary_mnpburden_module_summary.csv")
)

# Plot exact permutation distributions
observed_perm <- exact_perm_summary %>%
  select(module, observed_estimate)

p_perm <- ggplot(
  exact_perm_results,
  aes(x = perm_estimate)
) +
  geom_histogram(
    bins = 20,
    fill = "grey80",
    color = "white"
  ) +
  geom_vline(
    data = observed_perm,
    aes(xintercept = observed_estimate),
    color = "#D97706",
    linewidth = 0.9
  ) +
  facet_wrap(~ module, scales = "free") +
  labs(
    x = "Exact permuted term High − Low estimate",
    y = "Permutation count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black")
  )

ggsave(
  filename = file.path(outdir, "exact_permutation_primary_mnpburden_module_histograms.png"),
  plot = p_perm,
  width = 6.5,
  height = 3.2,
  units = "in",
  dpi = 600
)


# ============================================================
# Focused Xenium kNN neighborhood analysis for MNP spatial study
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(tidyverse)
  library(readr)
  library(FNN)
})

BASE <- "/media/jochum00/Aagaard_Raid3/ebarrozo/seurat_all_ref_merged/with_annotations"

outdir <- "xenium_neighborhood_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Xenium sample manifest
# ------------------------------------------------------------

xenium_manifest <- tribble(
  ~sample_id_std, ~rds_file,                                  ~fov,
  "HP13b_xen",    "seurat_HP13b_xen_withPredCols.rds",         "S23b",
  "HP2_xen",      "seurat_HP2_xen_withPredCols.rds",           "S31542",
  "HT1_xen",      "seurat_HT1_xen_withPredCols.rds",           "S11231",
  "HT2a_xen",     "seurat_HT2a_xen_withPredCols.rds",          "S13164a",
  "HT2b_xen",     "seurat_HT2b_withPredCols.rds",              "S13164b",
  "HT3_xen",      "seurat_HT3_xen_withPredCols.rds",           "S3065",
  "HP1_xen",      "seurat_HP1_withPredCols.rds",               "S9380",
  "LP2_xen",      "seurat_LP2_withPredCols.rds",               "S10963",
  "LT2_xen",      "seurat_LT2_withPredCols.rds",               "S30578",
  "NC4_xen",      "seurat_NC4_xen_withPredCols.rds",           "S26",
  "LT1_xen",      "seurat_LT1_xen_withPredCols.rds",           "S31175",
  "LP1_xen",      "seurat_LP1_xen_withPredCols.rds",           "S12129"
) %>%
  mutate(rds_path = file.path(BASE, rds_file))

# ------------------------------------------------------------
# 2. Sample-level metadata from merged object
# ------------------------------------------------------------

sample_meta_for_neighborhood <- obj@meta.data %>%
  as_tibble() %>%
  distinct(
    sample_id_std,
    study_id,
    platform_clean,
    Delivery_class,
    MNPBurden,
    MNPBurden_Delivery_class,
    cumulative_mnp_ug_per_g,
    cumulative_mnp_log10ugpg,
    mnp_model_sample
  ) %>%
  filter(platform_clean == "Xenium") %>%
  mutate(
    sample_id_std = as.character(sample_id_std),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

write_csv(
  sample_meta_for_neighborhood,
  file.path(outdir, "sample_meta_for_neighborhood.csv")
)

cat("\nSample-level metadata available:\n")
sample_meta_for_neighborhood %>%
  count(sample_id_std, study_id, Delivery_class, MNPBurden, mnp_model_sample) %>%
  arrange(Delivery_class, MNPBurden, sample_id_std) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

get_coords_by_fov_safe <- function(sample_obj, fov) {

  coords <- tryCatch(
    Seurat::GetTissueCoordinates(sample_obj, image = fov),
    error = function(e) NULL
  )

  if (!is.null(coords) && nrow(coords) > 0) {
    df <- as.data.frame(coords)

    if (!"cell" %in% names(df)) {
      df$cell <- rownames(df)
    }

    if (all(c("x", "y") %in% names(df))) {
      return(df %>% transmute(cell = as.character(cell), x = as.numeric(x), y = as.numeric(y)))
    }

    if (all(c("imagecol", "imagerow") %in% names(df))) {
      return(df %>% transmute(cell = as.character(cell), x = as.numeric(imagecol), y = as.numeric(imagerow)))
    }

    if (all(c("x_centroid", "y_centroid") %in% names(df))) {
      return(df %>% transmute(cell = as.character(cell), x = as.numeric(x_centroid), y = as.numeric(y_centroid)))
    }

    num_cols <- names(df)[vapply(df, is.numeric, TRUE)]
    if (length(num_cols) >= 2) {
      return(
        df %>%
          transmute(
            cell = as.character(cell),
            x = as.numeric(.data[[num_cols[1]]]),
            y = as.numeric(.data[[num_cols[2]]])
          )
      )
    }
  }

  img <- sample_obj@images[[fov]]
  bnd <- try(img@boundaries[["centroids"]], silent = TRUE)

  if (!inherits(bnd, "try-error") && !is.null(bnd)) {

    cell_vec <- try(slot(bnd, "cells"), silent = TRUE)
    coord_df <- NULL

    if ("coords" %in% slotNames(bnd)) {
      coord_df <- as.data.frame(slot(bnd, "coords"))
    } else if ("data" %in% slotNames(bnd)) {
      tmp <- as.data.frame(slot(bnd, "data"))
      num_cols <- names(tmp)[vapply(tmp, is.numeric, TRUE)]
      if (length(num_cols) >= 2) {
        coord_df <- tmp[, num_cols[1:2], drop = FALSE]
      }
    }

    if (!inherits(cell_vec, "try-error") &&
        !is.null(coord_df) &&
        nrow(coord_df) == length(cell_vec)) {

      names(coord_df)[1:2] <- c("x", "y")

      return(
        tibble(
          cell = as.character(cell_vec),
          x = as.numeric(coord_df$x),
          y = as.numeric(coord_df$y)
        )
      )
    }
  }

  stop("Could not extract coordinates for FOV: ", fov)
}

collapse_celltype <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "RCTD_lowconf"

  case_when(
    x %in% c("NK Cell", "Decidual NK Cell") ~ "NK family",
    x %in% c("Macrophage", "Hofbauer Cell") ~ "Macrophage/Hofbauer",
    x %in% c("Placental Endothelial") ~ "Placental Endothelial",
    x %in% c("Perivascular Cell") ~ "Perivascular Cell",
    x %in% c("Villous Cytotrophoblast") ~ "Villous Cytotrophoblast",
    x %in% c("Extravillous Trophoblast") ~ "Extravillous Trophoblast",
    x %in% c("Syncytiotrophoblast") ~ "Syncytiotrophoblast",
    x %in% c("ILC3") ~ "ILC3",
    x %in% c("Plasma Cell") ~ "Plasma Cell",
    x %in% c("T Cell", "Treg") ~ "T/Treg",
    x %in% c("B Cell") ~ "B Cell",
    x %in% c("Fibroblast") ~ "Fibroblast",
    TRUE ~ "Other/lowconf"
  )
}

target_groups <- c(
  "Macrophage/Hofbauer",
  "NK family",
  "Placental Endothelial"
)

neighbor_groups <- c(
  "Placental Endothelial",
  "Perivascular Cell",
  "Macrophage/Hofbauer",
  "NK family",
  "ILC3",
  "Plasma Cell",
  "Villous Cytotrophoblast",
  "Extravillous Trophoblast",
  "Syncytiotrophoblast",
  "Fibroblast"
)

# ------------------------------------------------------------
# 4. Extract cells from one sample object
# ------------------------------------------------------------

extract_one_sample_cells <- function(sample_id_std, rds_path, fov) {

  message("Loading sample: ", sample_id_std)

  if (!file.exists(rds_path)) {
    warning("Missing RDS path: ", rds_path)
    return(tibble())
  }

  sample_obj <- readRDS(rds_path)

  if (!fov %in% names(sample_obj@images)) {
    warning(
      "FOV not found for ", sample_id_std,
      ". Requested: ", fov,
      ". Available: ", paste(names(sample_obj@images), collapse = ", ")
    )
    return(tibble())
  }

  coords <- get_coords_by_fov_safe(sample_obj, fov)

  meta <- sample_obj@meta.data %>%
    rownames_to_column("cell") %>%
    as_tibble() %>%
    transmute(
      cell = as.character(cell),
      predicted.celltype = as.character(predicted.celltype),
      predicted.celltype.score = suppressWarnings(as.numeric(predicted.celltype.score))
    )

  coords %>%
    left_join(meta, by = "cell") %>%
    mutate(
      sample_id_std = sample_id_std,
      fov = fov,
      celltype_raw = predicted.celltype,
      celltype_collapsed = collapse_celltype(predicted.celltype)
    ) %>%
    select(
      sample_id_std,
      fov,
      cell,
      x,
      y,
      predicted.celltype,
      predicted.celltype.score,
      celltype_raw,
      celltype_collapsed
    )
}

xenium_cells_all <- pmap_dfr(
  xenium_manifest %>% select(sample_id_std, rds_path, fov),
  extract_one_sample_cells
) %>%
  as_tibble() %>%
  left_join(sample_meta_for_neighborhood, by = "sample_id_std")

write_csv(
  xenium_cells_all,
  file.path(outdir, "xenium_cells_all_with_sample_metadata.csv")
)

cat("\nExtracted cells by sample:\n")
xenium_cells_all %>%
  count(sample_id_std, study_id, Delivery_class, MNPBurden, mnp_model_sample, name = "n_cells") %>%
  arrange(Delivery_class, MNPBurden, sample_id_std) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 5. Keep MNP-model Xenium samples and focused cell types
# ------------------------------------------------------------

xenium_cells_mnp <- xenium_cells_all %>%
  filter(
    mnp_model_sample == TRUE,
    !is.na(Delivery_class),
    !is.na(MNPBurden),
    !is.na(x),
    !is.na(y),
    !is.na(celltype_collapsed),
    celltype_collapsed != "Other/lowconf"
  )

write_csv(
  xenium_cells_mnp,
  file.path(outdir, "xenium_cells_mnp_neighborhood_input.csv")
)

cat("\nMNP neighborhood input cells by sample:\n")
xenium_cells_mnp %>%
  count(sample_id_std, study_id, Delivery_class, MNPBurden, name = "n_cells") %>%
  arrange(Delivery_class, MNPBurden, sample_id_std) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 6. kNN neighborhood calculation by sample
# ------------------------------------------------------------

calculate_sample_neighborhoods <- function(df_sample, k = 30) {

  sid <- unique(df_sample$sample_id_std)
  message("Calculating kNN neighborhoods: ", sid)

  df_sample <- df_sample %>%
    mutate(row_index = row_number())

  all_xy <- as.matrix(df_sample[, c("x", "y")])
  all_group <- df_sample$celltype_collapsed

  target_df <- df_sample %>%
    filter(celltype_collapsed %in% target_groups)

  if (nrow(target_df) == 0) {
    warning("No target cells for sample: ", sid)
    return(tibble())
  }

  query_xy <- as.matrix(target_df[, c("x", "y")])

  # k + 1 because the query cell itself may appear as its nearest neighbor.
  nn <- FNN::get.knnx(
    data = all_xy,
    query = query_xy,
    k = min(k + 1, nrow(df_sample))
  )

  nn_idx_raw <- nn$nn.index
  self_idx <- target_df$row_index

  # Remove self if present and retain k neighbors.
  nn_idx <- t(vapply(seq_len(nrow(nn_idx_raw)), function(i) {
    vals <- nn_idx_raw[i, ]
    vals <- vals[vals != self_idx[i]]
    vals <- vals[seq_len(min(k, length(vals)))]
    if (length(vals) < k) vals <- c(vals, rep(NA_integer_, k - length(vals)))
    vals
  }, integer(k)))

  # For each neighbor group, calculate fraction among k nearest neighbors.
  neigh_frac <- map_dfc(neighbor_groups, function(g) {
    vals <- rowMeans(matrix(all_group[nn_idx], nrow = nrow(nn_idx)) == g, na.rm = TRUE)
    tibble(!!paste0("neighbor_frac__", make.names(g)) := vals)
  })

  target_df %>%
    select(
      sample_id_std,
      study_id,
      Delivery_class,
      MNPBurden,
      MNPBurden_Delivery_class,
      cumulative_mnp_ug_per_g,
      cumulative_mnp_log10ugpg,
      cell,
      target_celltype = celltype_collapsed
    ) %>%
    bind_cols(neigh_frac)
}

k_use <- 30

cell_level_neighborhoods <- xenium_cells_mnp %>%
  group_by(sample_id_std) %>%
  group_split() %>%
  map_dfr(calculate_sample_neighborhoods, k = k_use)

write_csv(
  cell_level_neighborhoods,
  file.path(outdir, paste0("cell_level_target_neighborhoods_k", k_use, ".csv"))
)

cat("\nTarget-cell neighborhood rows:\n")
cell_level_neighborhoods %>%
  count(sample_id_std, Delivery_class, MNPBurden, target_celltype, name = "n_target_cells") %>%
  arrange(Delivery_class, MNPBurden, sample_id_std, target_celltype) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 7. Aggregate to sample level
# ------------------------------------------------------------

neighbor_frac_cols <- grep("^neighbor_frac__", colnames(cell_level_neighborhoods), value = TRUE)

sample_level_neighborhoods <- cell_level_neighborhoods %>%
  group_by(
    sample_id_std,
    study_id,
    Delivery_class,
    MNPBurden,
    MNPBurden_Delivery_class,
    cumulative_mnp_ug_per_g,
    cumulative_mnp_log10ugpg,
    target_celltype
  ) %>%
  summarise(
    n_target_cells = n(),
    across(all_of(neighbor_frac_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

write_csv(
  sample_level_neighborhoods,
  file.path(outdir, paste0("sample_level_target_neighborhoods_k", k_use, ".csv"))
)

# ------------------------------------------------------------
# 8. Aggregate to subject level
# ------------------------------------------------------------

subject_level_neighborhoods <- sample_level_neighborhoods %>%
  mutate(weight_n_target_cells = n_target_cells) %>%
  group_by(
    study_id,
    Delivery_class,
    MNPBurden,
    MNPBurden_Delivery_class,
    target_celltype
  ) %>%
  summarise(
    n_samples = n(),
    total_target_cells = sum(weight_n_target_cells, na.rm = TRUE),
    across(
      all_of(neighbor_frac_cols),
      ~ weighted.mean(.x, w = weight_n_target_cells, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  rename(n_target_cells = total_target_cells)

write_csv(
  subject_level_neighborhoods,
  file.path(outdir, paste0("subject_level_target_neighborhoods_k", k_use, ".csv"))
)

cat("\nSubject-level neighborhood summary:\n")
subject_level_neighborhoods %>%
  count(study_id, Delivery_class, MNPBurden, target_celltype) %>%
  arrange(Delivery_class, MNPBurden, study_id, target_celltype) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 9. Long format and within-delivery contrasts
# ------------------------------------------------------------

subject_neighborhood_long <- subject_level_neighborhoods %>%
  pivot_longer(
    cols = all_of(neighbor_frac_cols),
    names_to = "neighbor_metric",
    values_to = "neighbor_fraction"
  ) %>%
  mutate(
    neighbor_celltype = neighbor_metric %>%
      str_replace("^neighbor_frac__", "") %>%
      str_replace_all("\\.", " "),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

write_csv(
  subject_neighborhood_long,
  file.path(outdir, paste0("subject_level_target_neighborhoods_long_k", k_use, ".csv"))
)

# Simple within-delivery high-low contrasts.
# This is intentionally descriptive/small-n, matching the rest of this section.
within_delivery_neighborhood_contrasts <- subject_neighborhood_long %>%
  group_by(target_celltype, neighbor_celltype, Delivery_class) %>%
  filter(n_distinct(MNPBurden) == 2) %>%
  summarise(
    n_subjects = n_distinct(study_id),
    mean_low = mean(neighbor_fraction[MNPBurden == "Low"], na.rm = TRUE),
    mean_high = mean(neighbor_fraction[MNPBurden == "High"], na.rm = TRUE),
    estimate_high_minus_low = mean_high - mean_low,
    p_value = tryCatch(
      t.test(neighbor_fraction ~ MNPBurden)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)

write_csv(
  within_delivery_neighborhood_contrasts,
  file.path(outdir, paste0("within_delivery_neighborhood_contrasts_k", k_use, ".csv"))
)

cat("\nTop neighborhood contrasts:\n")
within_delivery_neighborhood_contrasts %>%
  arrange(p_value) %>%
  slice_head(n = 40) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 10. Focused heatmap
# ------------------------------------------------------------

plot_df <- within_delivery_neighborhood_contrasts %>%
  filter(
    target_celltype %in% target_groups,
    neighbor_celltype %in% neighbor_groups
  ) %>%
  mutate(
    target_celltype = factor(target_celltype, levels = target_groups),
    neighbor_celltype = factor(neighbor_celltype, levels = neighbor_groups),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    sig_label = case_when(
      p_value < 0.05 ~ "●",
      p_value < 0.10 ~ "·",
      TRUE ~ ""
    )
  )

max_abs <- max(abs(plot_df$estimate_high_minus_low), na.rm = TRUE)

p_neigh_heatmap <- ggplot(
  plot_df,
  aes(x = neighbor_celltype, y = target_celltype, fill = estimate_high_minus_low)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3.2, color = "black") +
  facet_wrap(~ Delivery_class, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = "High − Low\nneighbor fraction"
  ) +
  labs(
    x = "Neighbor cell type",
    y = "Target cell type"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black", size = 11),
    legend.position = "right"
  )

p_neigh_heatmap

ggsave(
  filename = file.path(outdir, paste0("focused_xenium_neighborhood_high_low_heatmap_k", k_use, ".png")),
  plot = p_neigh_heatmap,
  width = 7.4,
  height = 3.4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------
# 11. Focused box/dot plot for the most biologically relevant axes
# ------------------------------------------------------------

focused_pairs <- tribble(
  ~target_celltype,          ~neighbor_celltype,
  "Macrophage/Hofbauer",    "Placental Endothelial",
  "Macrophage/Hofbauer",    "Perivascular Cell",
  "Macrophage/Hofbauer",    "NK family",
  "NK family",              "Placental Endothelial",
  "NK family",              "Macrophage/Hofbauer",
  "Placental Endothelial",  "Macrophage/Hofbauer",
  "Placental Endothelial",  "NK family"
)

focused_plot_df <- subject_neighborhood_long %>%
  inner_join(focused_pairs, by = c("target_celltype", "neighbor_celltype")) %>%
  mutate(
    pair_label = paste0(target_celltype, " → ", neighbor_celltype),
    pair_label = factor(pair_label, levels = paste0(focused_pairs$target_celltype, " → ", focused_pairs$neighbor_celltype)),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High")),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm"))
  )

p_focused <- ggplot(
  focused_plot_df,
  aes(x = MNPBurden, y = neighbor_fraction)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
  geom_point(
    aes(group = study_id),
    position = position_jitter(width = 0.08, height = 0),
    size = 2.0,
    alpha = 0.9
  ) +
  facet_grid(pair_label ~ Delivery_class, scales = "free_y") +
  labs(
    x = "MNP burden",
    y = paste0("Fraction of ", k_use, " nearest neighbors")
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

p_focused

ggsave(
  filename = file.path(outdir, paste0("focused_xenium_neighborhood_dotplots_k", k_use, ".png")),
  plot = p_focused,
  width = 7.2,
  height = 7.8,
  units = "in",
  dpi = 600
)

cat("\nDone. Outputs written to: ", outdir, "\n")













# ============================================================
# Forced-annotation Xenium kNN neighborhood analysis
# Uses sample-specific Xenium objects for centroids/native cell IDs
# Uses mnp_obj_w_nn_predicted_celltype.rds for forced cell labels
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(tidyverse)
  library(readr)
  library(FNN)
})


# ============================================================
# Force/rescue annotations inside image-bearing Xenium objects
# Purpose:
#   Preserve native Xenium cell IDs + images + centroids while adding
#   nn.predicted.celltype for neighborhood analysis.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(tidyverse)
  library(RcppAnnoy)
})

BASE <- "/media/jochum00/Aagaard_Raid3/ebarrozo/seurat_all_ref_merged/with_annotations"

OUT_DIR <- file.path(BASE, "xenium_objects_forced_annotation_native_ids")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

xenium_manifest <- tribble(
  ~sample_id_std, ~rds_file,                                  ~fov,
  "HP13b_xen",    "seurat_HP13b_xen_withPredCols.rds",         "S23b",
  "HP2_xen",      "seurat_HP2_xen_withPredCols.rds",           "S31542",
  "HT1_xen",      "seurat_HT1_xen_withPredCols.rds",           "S11231",
  "HT2a_xen",     "seurat_HT2a_xen_withPredCols.rds",          "S13164a",
  "HT2b_xen",     "seurat_HT2b_withPredCols.rds",              "S13164b",
  "HT3_xen",      "seurat_HT3_xen_withPredCols.rds",           "S3065",
  "HP1_xen",      "seurat_HP1_withPredCols.rds",               "S9380",
  "LP2_xen",      "seurat_LP2_withPredCols.rds",               "S10963",
  "LT2_xen",      "seurat_LT2_withPredCols.rds",               "S30578",
  "NC4_xen",      "seurat_NC4_xen_withPredCols.rds",           "S26",
  "LT1_xen",      "seurat_LT1_xen_withPredCols.rds",           "S31175",
  "LP1_xen",      "seurat_LP1_xen_withPredCols.rds",           "S12129"
) %>%
  mutate(
    rds_path = file.path(BASE, rds_file),
    forced_rds_path = file.path(
      OUT_DIR,
      paste0(sample_id_std, "_with_native_forced_nn_predicted_celltype.rds")
    )
  )

is_lowconf_label <- function(x) {
  x <- as.character(x)
  is.na(x) |
    x == "" |
    x %in% c("RCTD_lowconf", "RCTD-lowconf") |
    grepl("^UMAP_", x)
}

prepare_pca_for_sample <- function(obj, npcs = 30, nfeatures = 2000) {
  if ("pca" %in% Reductions(obj)) return(obj)

  assay_use <- DefaultAssay(obj)
  message("  Default assay: ", assay_use)

  obj <- NormalizeData(obj, verbose = FALSE)

  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = FALSE
  )

  obj <- ScaleData(
    obj,
    features = VariableFeatures(obj),
    verbose = FALSE
  )

  obj <- RunPCA(
    obj,
    features = VariableFeatures(obj),
    npcs = npcs,
    verbose = FALSE
  )

  obj
}

rescue_lowconf_one_object <- function(
  obj,
  celltype_col = "predicted.celltype",
  output_col = "nn.predicted.celltype",
  support_col = "nn.predicted.celltype_support",
  source_col = "nn.predicted.celltype_source",
  reduction = "pca",
  dims = 1:30,
  K = 75,
  n_trees = 50,
  support_min = 0.40,
  assign_all_remaining = TRUE,
  chunk = 20000,
  min_hi = 200
) {
  stopifnot(celltype_col %in% colnames(obj@meta.data))
  stopifnot(reduction %in% Reductions(obj))

  E <- Embeddings(obj, reduction)
  dims <- dims[dims <= ncol(E)]
  E <- E[, dims, drop = FALSE]

  lab0 <- as.character(obj@meta.data[[celltype_col]])
  bad0 <- is_lowconf_label(lab0)

  new_lab <- lab0
  support <- rep(NA_real_, length(lab0))
  source <- ifelse(bad0, "lowconf_initial", "highconf_original")

  hi_idx <- which(!bad0)
  lo_idx <- which(bad0)

  message("  high-confidence seeds: ", length(hi_idx))
  message("  low-confidence/NA cells to rescue: ", length(lo_idx))

  if (length(lo_idx) == 0) {
    obj@meta.data[[output_col]] <- new_lab
    obj@meta.data[[support_col]] <- 1
    obj@meta.data[[source_col]] <- source
    return(obj)
  }

  if (length(hi_idx) < min_hi) {
    warning("  Too few high-confidence seeds. Copying original labels only.")
    obj@meta.data[[output_col]] <- new_lab
    obj@meta.data[[support_col]] <- ifelse(!bad0, 1, NA_real_)
    obj@meta.data[[source_col]] <- source
    return(obj)
  }

  hi_labels <- new_lab[hi_idx]

  ann <- new(AnnoyEuclidean, ncol(E))

  for (i in seq_along(hi_idx)) {
    ann$addItem(i - 1L, as.numeric(E[hi_idx[i], ]))
  }

  ann$build(n_trees)

  seq_starts <- seq(1, length(lo_idx), by = chunk)

  for (s in seq_starts) {
    sel <- lo_idx[s:min(s + chunk - 1, length(lo_idx))]

    message(
      "  rescue chunk ",
      s,
      " to ",
      min(s + chunk - 1, length(lo_idx)),
      " of ",
      length(lo_idx)
    )

    for (j in sel) {
      v <- as.numeric(E[j, ])

      nn <- tryCatch(
        ann$getNNsByVector(v, K, -1, TRUE),
        error = function(e) ann$getNNsByVectorList(v, K, -1, TRUE)
      )

      nn_i <- nn$item
      nn_d <- nn$distance

      labs <- hi_labels[nn_i + 1L]
      w <- 1 / (nn_d + 1e-6)

      tab <- tapply(w, labs, sum)
      best <- names(tab)[which.max(tab)]
      sup <- max(tab) / sum(tab)

      support[j] <- sup

      if (!is.na(sup) && (sup >= support_min || assign_all_remaining)) {
        new_lab[j] <- best
        source[j] <- ifelse(sup >= support_min, "rescued_support_pass", "forced_low_support")
      }
    }
  }

  support[!bad0 & is.na(support)] <- 1

  obj@meta.data[[paste0(celltype_col, "_orig")]] <- lab0
  obj@meta.data[[output_col]] <- new_lab
  obj@meta.data[[support_col]] <- support
  obj@meta.data[[source_col]] <- source

  obj
}

force_one_xenium_object <- function(sample_id_std, rds_path, fov, forced_rds_path) {
  message("\n==============================")
  message("Processing ", sample_id_std)
  message("==============================")

  obj_s <- readRDS(rds_path)

  if (!fov %in% names(obj_s@images)) {
    warning("FOV not found for ", sample_id_std, ": ", fov)
  }

  if (!"predicted.celltype" %in% colnames(obj_s@meta.data)) {
    stop("Missing predicted.celltype in ", sample_id_std)
  }

  obj_s$predicted.celltype[obj_s$predicted.celltype %in% c("UMAP_1", "UMAP_2")] <- "RCTD_lowconf"

  before_tab <- tibble(
    sample_id_std = sample_id_std,
    n_cells = ncol(obj_s),
    n_lowconf_before = sum(is_lowconf_label(obj_s$predicted.celltype)),
    frac_lowconf_before = mean(is_lowconf_label(obj_s$predicted.celltype))
  )

  obj_s <- prepare_pca_for_sample(obj_s, npcs = 30, nfeatures = 2000)

  obj_s <- rescue_lowconf_one_object(
    obj_s,
    celltype_col = "predicted.celltype",
    output_col = "nn.predicted.celltype",
    support_col = "nn.predicted.celltype_support",
    source_col = "nn.predicted.celltype_source",
    reduction = "pca",
    dims = 1:30,
    K = 75,
    support_min = 0.40,
    assign_all_remaining = TRUE,
    min_hi = 200
  )

  after_tab <- tibble(
    sample_id_std = sample_id_std,
    n_lowconf_after = sum(is_lowconf_label(obj_s$nn.predicted.celltype)),
    frac_lowconf_after = mean(is_lowconf_label(obj_s$nn.predicted.celltype)),
    median_support = median(obj_s$nn.predicted.celltype_support, na.rm = TRUE),
    frac_low_support_lt_040 = mean(obj_s$nn.predicted.celltype_support < 0.40, na.rm = TRUE)
  )

  qc <- bind_cols(before_tab, after_tab %>% select(-sample_id_std))

  saveRDS(obj_s, forced_rds_path)

  qc
}

force_annotation_qc <- pmap_dfr(
  xenium_manifest %>% select(sample_id_std, rds_path, fov, forced_rds_path),
  force_one_xenium_object
)

write_csv(
  force_annotation_qc,
  file.path(OUT_DIR, "native_xenium_forced_annotation_qc.csv")
)

force_annotation_qc %>%
  arrange(desc(frac_low_support_lt_040)) %>%
  print(n = Inf)

# ============================================================
# Native-ID forced-annotation Xenium kNN neighborhood analysis
#
# Correct workflow:
#   1. Force/rescue annotations inside image-bearing Xenium objects.
#   2. Preserve native Xenium cell IDs, images, and centroids.
#   3. Use current metadata-updated merged obj only for sample-level metadata.
#   4. Do NOT join merged-object numeric cell IDs back to centroids.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(tidyverse)
  library(readr)
  library(FNN)
  library(RcppAnnoy)
})

BASE <- "/media/jochum00/Aagaard_Raid3/ebarrozo/seurat_all_ref_merged/with_annotations"

FORCED_NATIVE_DIR <- file.path(BASE, "xenium_objects_forced_annotation_native_ids")
dir.create(FORCED_NATIVE_DIR, showWarnings = FALSE, recursive = TRUE)

outdir <- "xenium_neighborhood_outputs_forced_native"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# obj must already be loaded and metadata-updated.
# It should include: sample_id_std, platform_clean, study_id, Delivery_class,
# MNPBurden, MNPBurden_Delivery_class, mnp_model_sample, and MNP burden columns.

required_meta_cols <- c(
  "sample_id_std",
  "platform_clean",
  "study_id",
  "Delivery_class",
  "MNPBurden",
  "MNPBurden_Delivery_class",
  "mnp_model_sample"
)

missing_meta_cols <- setdiff(required_meta_cols, colnames(obj@meta.data))

if (length(missing_meta_cols) > 0) {
  stop(
    "Current merged obj is missing required metadata columns: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 1. Xenium manifest
# ------------------------------------------------------------

xenium_manifest <- tribble(
  ~sample_id_std, ~rds_file,                                  ~fov,
  "HP13b_xen",    "seurat_HP13b_xen_withPredCols.rds",         "S23b",
  "HP2_xen",      "seurat_HP2_xen_withPredCols.rds",           "S31542",
  "HT1_xen",      "seurat_HT1_xen_withPredCols.rds",           "S11231",
  "HT2a_xen",     "seurat_HT2a_xen_withPredCols.rds",          "S13164a",
  "HT2b_xen",     "seurat_HT2b_withPredCols.rds",              "S13164b",
  "HT3_xen",      "seurat_HT3_xen_withPredCols.rds",           "S3065",
  "HP1_xen",      "seurat_HP1_withPredCols.rds",               "S9380",
  "LP2_xen",      "seurat_LP2_withPredCols.rds",               "S10963",
  "LT2_xen",      "seurat_LT2_withPredCols.rds",               "S30578",
  "NC4_xen",      "seurat_NC4_xen_withPredCols.rds",           "S26",
  "LT1_xen",      "seurat_LT1_xen_withPredCols.rds",           "S31175",
  "LP1_xen",      "seurat_LP1_xen_withPredCols.rds",           "S12129"
) %>%
  mutate(
    original_rds_path = file.path(BASE, rds_file),
    forced_rds_path = file.path(
      FORCED_NATIVE_DIR,
      paste0(sample_id_std, "_with_native_forced_nn_predicted_celltype.rds")
    )
  )

# ------------------------------------------------------------
# 2. Helpers for native-object forced annotation
# ------------------------------------------------------------

is_lowconf_label <- function(x) {
  x <- as.character(x)
  is.na(x) |
    x == "" |
    x %in% c("RCTD_lowconf", "RCTD-lowconf") |
    grepl("^UMAP_", x)
}

prepare_pca_for_sample <- function(obj_s, npcs = 30, nfeatures = 2000) {
  if ("pca" %in% Reductions(obj_s)) return(obj_s)

  message("  Default assay: ", DefaultAssay(obj_s))

  obj_s <- NormalizeData(obj_s, verbose = FALSE)

  obj_s <- FindVariableFeatures(
    obj_s,
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = FALSE
  )

  obj_s <- ScaleData(
    obj_s,
    features = VariableFeatures(obj_s),
    verbose = FALSE
  )

  obj_s <- RunPCA(
    obj_s,
    features = VariableFeatures(obj_s),
    npcs = npcs,
    verbose = FALSE
  )

  obj_s
}

rescue_lowconf_one_object <- function(
  obj_s,
  celltype_col = "predicted.celltype",
  output_col = "nn.predicted.celltype",
  support_col = "nn.predicted.celltype_support",
  source_col = "nn.predicted.celltype_source",
  reduction = "pca",
  dims = 1:30,
  K = 75,
  n_trees = 50,
  support_min = 0.40,
  assign_all_remaining = TRUE,
  chunk = 20000,
  min_hi = 200
) {
  stopifnot(celltype_col %in% colnames(obj_s@meta.data))
  stopifnot(reduction %in% Reductions(obj_s))

  E <- Embeddings(obj_s, reduction)
  dims <- dims[dims <= ncol(E)]
  E <- E[, dims, drop = FALSE]

  lab0 <- as.character(obj_s@meta.data[[celltype_col]])
  bad0 <- is_lowconf_label(lab0)

  new_lab <- lab0
  support <- rep(NA_real_, length(lab0))
  source <- ifelse(bad0, "lowconf_initial", "highconf_original")

  hi_idx <- which(!bad0)
  lo_idx <- which(bad0)

  message("  high-confidence seeds: ", length(hi_idx))
  message("  low-confidence/NA cells to rescue: ", length(lo_idx))

  if (length(lo_idx) == 0) {
    obj_s@meta.data[[output_col]] <- new_lab
    obj_s@meta.data[[support_col]] <- 1
    obj_s@meta.data[[source_col]] <- source
    return(obj_s)
  }

  if (length(hi_idx) < min_hi) {
    warning("  Too few high-confidence seeds. Copying original labels only.")
    obj_s@meta.data[[output_col]] <- new_lab
    obj_s@meta.data[[support_col]] <- ifelse(!bad0, 1, NA_real_)
    obj_s@meta.data[[source_col]] <- source
    return(obj_s)
  }

  hi_labels <- new_lab[hi_idx]

  ann <- new(AnnoyEuclidean, ncol(E))

  for (i in seq_along(hi_idx)) {
    ann$addItem(i - 1L, as.numeric(E[hi_idx[i], ]))
  }

  ann$build(n_trees)

  seq_starts <- seq(1, length(lo_idx), by = chunk)

  for (s in seq_starts) {
    sel <- lo_idx[s:min(s + chunk - 1, length(lo_idx))]

    message(
      "  rescue chunk ",
      s,
      " to ",
      min(s + chunk - 1, length(lo_idx)),
      " of ",
      length(lo_idx)
    )

    for (j in sel) {
      v <- as.numeric(E[j, ])

      nn <- tryCatch(
        ann$getNNsByVector(v, K, -1, TRUE),
        error = function(e) ann$getNNsByVectorList(v, K, -1, TRUE)
      )

      nn_i <- nn$item
      nn_d <- nn$distance

      labs <- hi_labels[nn_i + 1L]
      w <- 1 / (nn_d + 1e-6)

      tab <- tapply(w, labs, sum)
      best <- names(tab)[which.max(tab)]
      sup <- max(tab) / sum(tab)

      support[j] <- sup

      if (!is.na(sup) && (sup >= support_min || assign_all_remaining)) {
        new_lab[j] <- best
        source[j] <- ifelse(
          sup >= support_min,
          "rescued_support_pass",
          "forced_low_support"
        )
      }
    }
  }

  support[!bad0 & is.na(support)] <- 1

  obj_s@meta.data[[paste0(celltype_col, "_orig")]] <- lab0
  obj_s@meta.data[[output_col]] <- new_lab
  obj_s@meta.data[[support_col]] <- support
  obj_s@meta.data[[source_col]] <- source

  obj_s
}

force_one_xenium_object <- function(sample_id_std, original_rds_path, fov, forced_rds_path) {
  message("\n==============================")
  message("Processing ", sample_id_std)
  message("==============================")

  if (!file.exists(original_rds_path)) {
    warning("Missing original RDS: ", original_rds_path)
    return(tibble(
      sample_id_std = sample_id_std,
      status = "missing_original_rds"
    ))
  }

  obj_s <- readRDS(original_rds_path)

  if (!fov %in% names(obj_s@images)) {
    warning("FOV not found for ", sample_id_std, ": ", fov)
  }

  if (!"predicted.celltype" %in% colnames(obj_s@meta.data)) {
    stop("Missing predicted.celltype in ", sample_id_std)
  }

  obj_s$predicted.celltype[obj_s$predicted.celltype %in% c("UMAP_1", "UMAP_2")] <- "RCTD_lowconf"

  before_tab <- tibble(
    sample_id_std = sample_id_std,
    status = "ok",
    n_cells = ncol(obj_s),
    n_lowconf_before = sum(is_lowconf_label(obj_s$predicted.celltype)),
    frac_lowconf_before = mean(is_lowconf_label(obj_s$predicted.celltype))
  )

  obj_s <- prepare_pca_for_sample(obj_s, npcs = 30, nfeatures = 2000)

  obj_s <- rescue_lowconf_one_object(
    obj_s,
    celltype_col = "predicted.celltype",
    output_col = "nn.predicted.celltype",
    support_col = "nn.predicted.celltype_support",
    source_col = "nn.predicted.celltype_source",
    reduction = "pca",
    dims = 1:30,
    K = 75,
    support_min = 0.40,
    assign_all_remaining = TRUE,
    min_hi = 200
  )

  after_tab <- tibble(
    n_lowconf_after = sum(is_lowconf_label(obj_s$nn.predicted.celltype)),
    frac_lowconf_after = mean(is_lowconf_label(obj_s$nn.predicted.celltype)),
    median_support = median(obj_s$nn.predicted.celltype_support, na.rm = TRUE),
    frac_low_support_lt_040 = mean(obj_s$nn.predicted.celltype_support < 0.40, na.rm = TRUE)
  )

  saveRDS(obj_s, forced_rds_path)

  bind_cols(before_tab, after_tab)
}

# ------------------------------------------------------------
# 3. Run forced annotation if needed
# ------------------------------------------------------------

run_force_annotation <- TRUE

if (run_force_annotation) {
  force_annotation_qc <- pmap_dfr(
    xenium_manifest %>%
      select(sample_id_std, original_rds_path, fov, forced_rds_path),
    force_one_xenium_object
  )

  write_csv(
    force_annotation_qc,
    file.path(FORCED_NATIVE_DIR, "native_xenium_forced_annotation_qc.csv")
  )

  force_annotation_qc %>%
    arrange(desc(frac_low_support_lt_040)) %>%
    as_tibble() %>%
    print(n = Inf)
}

# ------------------------------------------------------------
# 4. Sample-level metadata from current merged object
# ------------------------------------------------------------

meta_cols_available <- colnames(obj@meta.data)

optional_mnp_cols <- intersect(
  c(
    "cumulative_mnp_ug_per_g",
    "cumulative_mnp_log10ugpg",
    "Cumulative_MNPs_ug_per_g",
    "Cumulative_MNPs_log10ugpg"
  ),
  meta_cols_available
)

sample_meta_for_neighborhood <- obj@meta.data %>%
  as_tibble() %>%
  distinct(
    sample_id_std,
    study_id,
    platform_clean,
    Delivery_class,
    MNPBurden,
    MNPBurden_Delivery_class,
    mnp_model_sample,
    across(all_of(optional_mnp_cols))
  ) %>%
  filter(platform_clean == "Xenium") %>%
  mutate(
    sample_id_std = as.character(sample_id_std),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

# Standardize optional MNP names for downstream select().
if (!"cumulative_mnp_ug_per_g" %in% colnames(sample_meta_for_neighborhood)) {
  if ("Cumulative_MNPs_ug_per_g" %in% colnames(sample_meta_for_neighborhood)) {
    sample_meta_for_neighborhood <- sample_meta_for_neighborhood %>%
      mutate(cumulative_mnp_ug_per_g = Cumulative_MNPs_ug_per_g)
  } else {
    sample_meta_for_neighborhood <- sample_meta_for_neighborhood %>%
      mutate(cumulative_mnp_ug_per_g = NA_real_)
  }
}

if (!"cumulative_mnp_log10ugpg" %in% colnames(sample_meta_for_neighborhood)) {
  if ("Cumulative_MNPs_log10ugpg" %in% colnames(sample_meta_for_neighborhood)) {
    sample_meta_for_neighborhood <- sample_meta_for_neighborhood %>%
      mutate(cumulative_mnp_log10ugpg = Cumulative_MNPs_log10ugpg)
  } else {
    sample_meta_for_neighborhood <- sample_meta_for_neighborhood %>%
      mutate(cumulative_mnp_log10ugpg = NA_real_)
  }
}

write_csv(
  sample_meta_for_neighborhood,
  file.path(outdir, "sample_meta_for_neighborhood_forced_native.csv")
)

cat("\nSample metadata available:\n")
sample_meta_for_neighborhood %>%
  count(sample_id_std, study_id, Delivery_class, MNPBurden, mnp_model_sample) %>%
  arrange(Delivery_class, MNPBurden, sample_id_std) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 5. Coordinate extraction helpers
# ------------------------------------------------------------

get_coords_by_fov_safe <- function(sample_obj, fov) {

  coords <- tryCatch(
    Seurat::GetTissueCoordinates(sample_obj, image = fov),
    error = function(e) NULL
  )

  if (!is.null(coords) && nrow(coords) > 0) {
    df <- as.data.frame(coords)

    if (!"cell" %in% names(df)) {
      df$cell <- rownames(df)
    }

    if (all(c("x", "y") %in% names(df))) {
      return(df %>% transmute(cell = as.character(cell), x = as.numeric(x), y = as.numeric(y)))
    }

    if (all(c("imagecol", "imagerow") %in% names(df))) {
      return(df %>% transmute(cell = as.character(cell), x = as.numeric(imagecol), y = as.numeric(imagerow)))
    }

    if (all(c("x_centroid", "y_centroid") %in% names(df))) {
      return(df %>% transmute(cell = as.character(cell), x = as.numeric(x_centroid), y = as.numeric(y_centroid)))
    }

    num_cols <- names(df)[vapply(df, is.numeric, TRUE)]

    if (length(num_cols) >= 2) {
      return(
        df %>%
          transmute(
            cell = as.character(cell),
            x = as.numeric(.data[[num_cols[1]]]),
            y = as.numeric(.data[[num_cols[2]]])
          )
      )
    }
  }

  img <- sample_obj@images[[fov]]
  bnd <- try(img@boundaries[["centroids"]], silent = TRUE)

  if (!inherits(bnd, "try-error") && !is.null(bnd)) {

    cell_vec <- try(slot(bnd, "cells"), silent = TRUE)
    coord_df <- NULL

    if ("coords" %in% slotNames(bnd)) {
      coord_df <- as.data.frame(slot(bnd, "coords"))
    } else if ("data" %in% slotNames(bnd)) {
      tmp <- as.data.frame(slot(bnd, "data"))
      num_cols <- names(tmp)[vapply(tmp, is.numeric, TRUE)]

      if (length(num_cols) >= 2) {
        coord_df <- tmp[, num_cols[1:2], drop = FALSE]
      }
    }

    if (!inherits(cell_vec, "try-error") &&
        !is.null(coord_df) &&
        nrow(coord_df) == length(cell_vec)) {

      names(coord_df)[1:2] <- c("x", "y")

      return(
        tibble(
          cell = as.character(cell_vec),
          x = as.numeric(coord_df$x),
          y = as.numeric(coord_df$y)
        )
      )
    }
  }

  stop("Could not extract coordinates for FOV: ", fov)
}

collapse_celltype <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "RCTD_lowconf"

  case_when(
    x %in% c("NK Cell", "Decidual NK Cell") ~ "NK family",
    x %in% c("Macrophage", "Hofbauer Cell") ~ "Macrophage/Hofbauer",
    x %in% c("Placental Endothelial") ~ "Placental Endothelial",
    x %in% c("Perivascular Cell") ~ "Perivascular Cell",
    x %in% c("Villous Cytotrophoblast") ~ "Villous Cytotrophoblast",
    x %in% c("Extravillous Trophoblast") ~ "Extravillous Trophoblast",
    x %in% c("Syncytiotrophoblast") ~ "Syncytiotrophoblast",
    x %in% c("ILC3") ~ "ILC3",
    x %in% c("Plasma Cell") ~ "Plasma Cell",
    x %in% c("T Cell", "Treg") ~ "T/Treg",
    x %in% c("B Cell") ~ "B Cell",
    x %in% c("Fibroblast") ~ "Fibroblast",
    TRUE ~ "Other/lowconf"
  )
}

target_groups <- c(
  "Macrophage/Hofbauer",
  "NK family",
  "Placental Endothelial"
)

neighbor_groups <- c(
  "Placental Endothelial",
  "Perivascular Cell",
  "Macrophage/Hofbauer",
  "NK family",
  "ILC3",
  "Plasma Cell",
  "Villous Cytotrophoblast",
  "Extravillous Trophoblast",
  "Syncytiotrophoblast",
  "Fibroblast"
)

# ------------------------------------------------------------
# 6. Extract cells from forced native objects
# ------------------------------------------------------------

extract_one_sample_cells_forced_native <- function(sample_id_std, forced_rds_path, fov) {

  message("Loading forced native sample: ", sample_id_std)

  if (!file.exists(forced_rds_path)) {
    warning("Missing forced native RDS: ", forced_rds_path)
    return(tibble())
  }

  sample_obj <- readRDS(forced_rds_path)

  if (!fov %in% names(sample_obj@images)) {
    warning(
      "FOV not found for ", sample_id_std,
      ". Requested: ", fov,
      ". Available: ", paste(names(sample_obj@images), collapse = ", ")
    )
    return(tibble())
  }

  if (!"nn.predicted.celltype" %in% colnames(sample_obj@meta.data)) {
    stop("Missing nn.predicted.celltype in forced native object: ", sample_id_std)
  }

  coords <- get_coords_by_fov_safe(sample_obj, fov)

  meta_native <- sample_obj@meta.data %>%
    rownames_to_column("cell") %>%
    as_tibble() %>%
    transmute(
      cell = as.character(cell),
      highconf_predicted_celltype = as.character(predicted.celltype),
      highconf_predicted_score = suppressWarnings(as.numeric(predicted.celltype.score)),
      nn.predicted.celltype = as.character(nn.predicted.celltype),
      nn.predicted.celltype_support = suppressWarnings(as.numeric(nn.predicted.celltype_support)),
      nn.predicted.celltype_source = as.character(nn.predicted.celltype_source)
    )

  coords %>%
    left_join(meta_native, by = "cell") %>%
    mutate(
      sample_id_std = sample_id_std,
      fov = fov,
      celltype_raw = nn.predicted.celltype,
      celltype_collapsed = collapse_celltype(nn.predicted.celltype),
      highconf_collapsed = collapse_celltype(highconf_predicted_celltype)
    ) %>%
    select(
      sample_id_std,
      fov,
      cell,
      x,
      y,
      highconf_predicted_celltype,
      highconf_predicted_score,
      highconf_collapsed,
      nn.predicted.celltype,
      nn.predicted.celltype_support,
      nn.predicted.celltype_source,
      celltype_raw,
      celltype_collapsed
    )
}

xenium_cells_all_forced <- pmap_dfr(
  xenium_manifest %>% select(sample_id_std, forced_rds_path, fov),
  extract_one_sample_cells_forced_native
) %>%
  as_tibble() %>%
  left_join(sample_meta_for_neighborhood, by = "sample_id_std")

write_csv(
  xenium_cells_all_forced,
  file.path(outdir, "xenium_cells_all_with_forced_native_labels_and_sample_metadata.csv")
)

cat("\nForced native cells by sample:\n")
xenium_cells_all_forced %>%
  count(sample_id_std, study_id, Delivery_class, MNPBurden, mnp_model_sample, name = "n_cells") %>%
  arrange(Delivery_class, MNPBurden, sample_id_std) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 7. Keep MNP-model samples
# ------------------------------------------------------------

support_threshold <- NA_real_
# support_threshold <- 0.40

xenium_cells_mnp_forced <- xenium_cells_all_forced %>%
  filter(
    mnp_model_sample == TRUE,
    !is.na(Delivery_class),
    !is.na(MNPBurden),
    !is.na(x),
    !is.na(y),
    !is.na(celltype_collapsed),
    celltype_collapsed != "Other/lowconf"
  )

if (!is.na(support_threshold)) {
  xenium_cells_mnp_forced <- xenium_cells_mnp_forced %>%
    filter(is.na(nn.predicted.celltype_support) | nn.predicted.celltype_support >= support_threshold)
}

write_csv(
  xenium_cells_mnp_forced,
  file.path(outdir, "xenium_cells_mnp_neighborhood_input_forced_native.csv")
)

cat("\nMNP forced-native neighborhood input cells by sample:\n")
xenium_cells_mnp_forced %>%
  count(sample_id_std, study_id, Delivery_class, MNPBurden, name = "n_cells") %>%
  arrange(Delivery_class, MNPBurden, sample_id_std) %>%
  as_tibble() %>%
  print(n = Inf)

cat("\nTarget-group counts by sample:\n")
xenium_cells_mnp_forced %>%
  filter(celltype_collapsed %in% target_groups) %>%
  count(sample_id_std, study_id, Delivery_class, MNPBurden, celltype_collapsed, name = "n_target_cells") %>%
  arrange(Delivery_class, MNPBurden, sample_id_std, celltype_collapsed) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 8. kNN neighborhood calculation
# ------------------------------------------------------------

calculate_sample_neighborhoods <- function(df_sample, k = 30) {

  sid <- unique(df_sample$sample_id_std)
  message("Calculating kNN neighborhoods: ", sid)

  df_sample <- df_sample %>%
    mutate(row_index = row_number())

  all_xy <- as.matrix(df_sample[, c("x", "y")])
  all_group <- df_sample$celltype_collapsed

  target_df <- df_sample %>%
    filter(celltype_collapsed %in% target_groups)

  if (nrow(target_df) == 0) {
    warning("No target cells for sample: ", sid)
    return(tibble())
  }

  query_xy <- as.matrix(target_df[, c("x", "y")])

  nn <- FNN::get.knnx(
    data = all_xy,
    query = query_xy,
    k = min(k + 1, nrow(df_sample))
  )

  nn_idx_raw <- nn$nn.index
  self_idx <- target_df$row_index

  nn_idx <- t(vapply(seq_len(nrow(nn_idx_raw)), function(i) {
    vals <- nn_idx_raw[i, ]
    vals <- vals[vals != self_idx[i]]
    vals <- vals[seq_len(min(k, length(vals)))]

    if (length(vals) < k) {
      vals <- c(vals, rep(NA_integer_, k - length(vals)))
    }

    vals
  }, integer(k)))

  neigh_frac <- map_dfc(neighbor_groups, function(g) {
    vals <- rowMeans(
      matrix(all_group[nn_idx], nrow = nrow(nn_idx)) == g,
      na.rm = TRUE
    )

    tibble(!!paste0("neighbor_frac__", make.names(g)) := vals)
  })

  target_df %>%
    select(
      sample_id_std,
      study_id,
      Delivery_class,
      MNPBurden,
      MNPBurden_Delivery_class,
      cumulative_mnp_ug_per_g,
      cumulative_mnp_log10ugpg,
      cell,
      nn.predicted.celltype_support,
      target_celltype = celltype_collapsed
    ) %>%
    bind_cols(neigh_frac)
}

k_use <- 30

cell_level_neighborhoods <- xenium_cells_mnp_forced %>%
  group_by(sample_id_std) %>%
  group_split() %>%
  map_dfr(calculate_sample_neighborhoods, k = k_use)

write_csv(
  cell_level_neighborhoods,
  file.path(outdir, paste0("cell_level_target_neighborhoods_forced_native_k", k_use, ".csv"))
)

# ------------------------------------------------------------
# 9. Aggregate sample level
# ------------------------------------------------------------

neighbor_frac_cols <- grep("^neighbor_frac__", colnames(cell_level_neighborhoods), value = TRUE)

sample_level_neighborhoods <- cell_level_neighborhoods %>%
  group_by(
    sample_id_std,
    study_id,
    Delivery_class,
    MNPBurden,
    MNPBurden_Delivery_class,
    cumulative_mnp_ug_per_g,
    cumulative_mnp_log10ugpg,
    target_celltype
  ) %>%
  summarise(
    n_target_cells = n(),
    across(all_of(neighbor_frac_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

write_csv(
  sample_level_neighborhoods,
  file.path(outdir, paste0("sample_level_target_neighborhoods_forced_native_k", k_use, ".csv"))
)

# ------------------------------------------------------------
# 10. Aggregate subject level
# ------------------------------------------------------------

subject_level_neighborhoods <- sample_level_neighborhoods %>%
  mutate(weight_n_target_cells = n_target_cells) %>%
  group_by(
    study_id,
    Delivery_class,
    MNPBurden,
    MNPBurden_Delivery_class,
    target_celltype
  ) %>%
  summarise(
    n_samples = n(),
    total_target_cells = sum(weight_n_target_cells, na.rm = TRUE),
    across(
      all_of(neighbor_frac_cols),
      ~ weighted.mean(.x, w = weight_n_target_cells, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  rename(n_target_cells = total_target_cells)

write_csv(
  subject_level_neighborhoods,
  file.path(outdir, paste0("subject_level_target_neighborhoods_forced_native_k", k_use, ".csv"))
)

# ------------------------------------------------------------
# 11. Long format and high-low contrasts
# ------------------------------------------------------------

subject_neighborhood_long <- subject_level_neighborhoods %>%
  pivot_longer(
    cols = all_of(neighbor_frac_cols),
    names_to = "neighbor_metric",
    values_to = "neighbor_fraction"
  ) %>%
  mutate(
    neighbor_celltype = neighbor_metric %>%
      str_replace("^neighbor_frac__", "") %>%
      str_replace_all("\\.", " "),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High"))
  )

write_csv(
  subject_neighborhood_long,
  file.path(outdir, paste0("subject_level_target_neighborhoods_long_forced_native_k", k_use, ".csv"))
)

within_delivery_neighborhood_contrasts <- subject_neighborhood_long %>%
  group_by(target_celltype, neighbor_celltype, Delivery_class) %>%
  filter(n_distinct(MNPBurden) == 2) %>%
  summarise(
    n_subjects = n_distinct(study_id),
    mean_low = mean(neighbor_fraction[MNPBurden == "Low"], na.rm = TRUE),
    mean_high = mean(neighbor_fraction[MNPBurden == "High"], na.rm = TRUE),
    estimate_high_minus_low = mean_high - mean_low,
    p_value = tryCatch(
      t.test(neighbor_fraction ~ MNPBurden)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write_csv(
  within_delivery_neighborhood_contrasts,
  file.path(outdir, paste0("within_delivery_neighborhood_contrasts_forced_native_k", k_use, ".csv"))
)

cat("\nTop forced-native neighborhood contrasts:\n")
within_delivery_neighborhood_contrasts %>%
  arrange(p_value) %>%
  slice_head(n = 50) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 12. Heatmap
# ------------------------------------------------------------

plot_df <- within_delivery_neighborhood_contrasts %>%
  filter(
    target_celltype %in% target_groups,
    neighbor_celltype %in% neighbor_groups
  ) %>%
  mutate(
    target_celltype = factor(target_celltype, levels = target_groups),
    neighbor_celltype = factor(neighbor_celltype, levels = neighbor_groups),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    sig_label = case_when(
      p_value < 0.05 ~ "●",
      p_value < 0.10 ~ "·",
      TRUE ~ ""
    )
  )

max_abs <- max(abs(plot_df$estimate_high_minus_low), na.rm = TRUE)

p_neigh_heatmap <- ggplot(
  plot_df,
  aes(x = neighbor_celltype, y = target_celltype, fill = estimate_high_minus_low)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3.2, color = "black") +
  facet_wrap(~ Delivery_class, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = "High − Low\nneighbor fraction"
  ) +
  labs(
    x = "Neighbor cell type",
    y = "Target cell type"
  ) +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black", size = 11),
    legend.position = "right"
  )

p_neigh_heatmap

ggsave(
  filename = file.path(outdir, paste0("focused_xenium_neighborhood_high_low_heatmap_forced_native_k", k_use, ".png")),
  plot = p_neigh_heatmap,
  width = 7.4,
  height = 3.4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------
# 13. Focused dotplots
# ------------------------------------------------------------

focused_pairs <- tribble(
  ~target_celltype,          ~neighbor_celltype,
  "Macrophage/Hofbauer",    "Placental Endothelial",
  "Macrophage/Hofbauer",    "Perivascular Cell",
  "Macrophage/Hofbauer",    "NK family",
  "NK family",              "Placental Endothelial",
  "NK family",              "Macrophage/Hofbauer",
  "Placental Endothelial",  "Macrophage/Hofbauer",
  "Placental Endothelial",  "NK family"
)

focused_plot_df <- subject_neighborhood_long %>%
  inner_join(focused_pairs, by = c("target_celltype", "neighbor_celltype")) %>%
  mutate(
    pair_label = paste0(target_celltype, " → ", neighbor_celltype),
    pair_label = factor(
      pair_label,
      levels = paste0(focused_pairs$target_celltype, " → ", focused_pairs$neighbor_celltype)
    ),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High")),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm"))
  )

p_focused <- ggplot(
  focused_plot_df,
  aes(x = MNPBurden, y = neighbor_fraction)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
  geom_point(
    aes(group = study_id),
    position = position_jitter(width = 0.08, height = 0),
    size = 2.0,
    alpha = 0.9
  ) +
  facet_grid(pair_label ~ Delivery_class, scales = "free_y") +
  labs(
    x = "MNP burden",
    y = paste0("Fraction of ", k_use, " nearest neighbors")
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

p_focused

ggsave(
  filename = file.path(outdir, paste0("focused_xenium_neighborhood_dotplots_forced_native_k", k_use, ".png")),
  plot = p_focused,
  width = 7.2,
  height = 7.8,
  units = "in",
  dpi = 600
)

cat("\nDone. Outputs written to: ", outdir, "\n")







### 

# ============================================================
# VisiumHD local module co-enrichment analysis
# MNP spatial transcriptomics manuscript
#
# Goal:
#   Test whether immune module scores locally co-enrich with
#   vascular-interface or trophoblast-remodeling module scores
#   in VisiumHD bins, then compare high vs low MNP burden
#   within delivery class.
#
# Unit:
#   VisiumHD bin/spot, not single cell.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(tidyverse)
  library(readr)
  library(FNN)
  library(broom)
})

BASE <- "/media/jochum00/Aagaard_Raid3/ebarrozo/seurat_all_ref_merged/with_annotations"

outdir <- "visiumhd_local_module_coenrichment_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. VisiumHD manifest
# ------------------------------------------------------------
# Update fov names if your HD objects use slightly different image names.
# These fov names came from the prior remap script comments.

visiumhd_manifest <- tribble(
  ~sample_id_std, ~rds_file,                         ~fov,
  "HT2_hd",       "seurat_HT2_hd_withPredCols.rds",   "HT2_hd.008um",
  "HT3_hd",       "seurat_HT3_withPredCols.rds",      "HT3.008um",
  "HP2_hd",       "seurat_HP2_hd_withPredCols.rds",   "HP2.008um",
  "LP1_hd",       "seurat_LP1_hd_withPredCols.rds",   "LP1_hd.008um"
) %>%
  mutate(rds_path = file.path(BASE, rds_file))

# ------------------------------------------------------------
# 2. Required metadata/module columns
# ------------------------------------------------------------

module_score_cols <- c(
  "ECM_Remodeling_score",
  "Vascular_Interface_score",
  "Macrophage_Remodeling_score",
  "NK_Cytotoxic_Regulatory_score",
  "Trophoblast_Remodeling_score"
)

# Prefer platform-normalized module scores if they exist.
module_z_cols <- paste0(module_score_cols, "_platform_z")

available_z <- intersect(module_z_cols, colnames(obj@meta.data))
available_raw <- intersect(module_score_cols, colnames(obj@meta.data))

if (length(available_z) == length(module_score_cols)) {
  module_use_cols <- module_z_cols
  names(module_use_cols) <- module_score_cols
  message("Using platform-normalized module scores.")
} else if (length(available_raw) == length(module_score_cols)) {
  module_use_cols <- module_score_cols
  names(module_use_cols) <- module_score_cols
  message("Using raw module scores. Consider platform-z scoring later.")
} else {
  stop(
    "Missing module score columns. Available raw: ",
    paste(available_raw, collapse = ", "),
    " | Available z: ",
    paste(available_z, collapse = ", ")
  )
}

required_meta_cols <- c(
  "sample_id_std",
  "study_id",
  "platform_clean",
  "Delivery_class",
  "MNPBurden",
  "MNPBurden_Delivery_class",
  "mnp_model_sample"
)

missing_meta <- setdiff(required_meta_cols, colnames(obj@meta.data))

if (length(missing_meta) > 0) {
  stop("Missing metadata columns in obj: ", paste(missing_meta, collapse = ", "))
}

# ------------------------------------------------------------
# 3. Coordinate helper
# ------------------------------------------------------------

get_coords_by_fov_safe <- function(sample_obj, fov) {

  coords <- tryCatch(
    Seurat::GetTissueCoordinates(sample_obj, image = fov),
    error = function(e) NULL
  )

  if (!is.null(coords) && nrow(coords) > 0) {
    df <- as.data.frame(coords)

    if (!"cell" %in% names(df)) {
      df$cell <- rownames(df)
    }

    if (all(c("x", "y") %in% names(df))) {
      return(df %>% transmute(cell = as.character(cell), x = as.numeric(x), y = as.numeric(y)))
    }

    if (all(c("imagecol", "imagerow") %in% names(df))) {
      return(df %>% transmute(cell = as.character(cell), x = as.numeric(imagecol), y = as.numeric(imagerow)))
    }

    if (all(c("pxl_col_in_fullres", "pxl_row_in_fullres") %in% names(df))) {
      return(df %>% transmute(cell = as.character(cell), x = as.numeric(pxl_col_in_fullres), y = as.numeric(pxl_row_in_fullres)))
    }

    num_cols <- names(df)[vapply(df, is.numeric, TRUE)]

    if (length(num_cols) >= 2) {
      return(
        df %>%
          transmute(
            cell = as.character(cell),
            x = as.numeric(.data[[num_cols[1]]]),
            y = as.numeric(.data[[num_cols[2]]])
          )
      )
    }
  }

  img <- sample_obj@images[[fov]]

  if ("coordinates" %in% slotNames(img)) {
    coords <- try(slot(img, "coordinates"), silent = TRUE)

    if (!inherits(coords, "try-error") && !is.null(coords) && nrow(coords) > 0) {
      df <- as.data.frame(coords)
      df$cell <- rownames(df)

      if (all(c("imagecol", "imagerow") %in% names(df))) {
        return(df %>% transmute(cell = as.character(cell), x = as.numeric(imagecol), y = as.numeric(imagerow)))
      }

      if (all(c("pxl_col_in_fullres", "pxl_row_in_fullres") %in% names(df))) {
        return(df %>% transmute(cell = as.character(cell), x = as.numeric(pxl_col_in_fullres), y = as.numeric(pxl_row_in_fullres)))
      }

      num_cols <- names(df)[vapply(df, is.numeric, TRUE)]

      if (length(num_cols) >= 2) {
        return(
          df %>%
            transmute(
              cell = as.character(cell),
              x = as.numeric(.data[[num_cols[1]]]),
              y = as.numeric(.data[[num_cols[2]]])
            )
        )
      }
    }
  }

  stop("Could not extract coordinates for FOV: ", fov)
}

# ------------------------------------------------------------
# 4. Extract HD coordinates and join current merged metadata
# ------------------------------------------------------------
# The merge-rowname problem may also affect HD. First try sample_id_std + cell.
# If that fails, we fall back to sample object metadata for scores if present.
# ------------------------------------------------------------
# 4. Compute module scores directly in each VisiumHD object
# ------------------------------------------------------------

module_gene_sets <- list(
  ECM_Remodeling = c(
    "COL1A1","COL1A2","COL3A1","COL4A1","COL4A2",
    "FN1","SPARC","TIMP3","MMP14","SERPINE1"
  ),
  Vascular_Interface = c(
    "PECAM1","VWF","KDR","FLT1","ENG","ESAM","NRP1","PGF"
  ),
  Macrophage_Remodeling = c(
    "CD68","CD163","VSIG4","FCGR1A","FCGR3A",
    "C1QA","C1QB","C1QC","MRC1","MSR1"
  ),
  NK_Cytotoxic_Regulatory = c(
    "NKG7","GNLY","PRF1","GZMB","KLRD1","KLRB1","XCL1","XCL2"
  ),
  Trophoblast_Remodeling = c(
    "KRT7","GATA3","HLA-G","FLT1","ENG","TFPI","LGALS13","ERVFRD-1"
  )
)

score_one_module <- function(obj_s, genes, score_name, assay_use = NULL) {
  if (is.null(assay_use)) {
    assay_use <- if ("SCT" %in% Assays(obj_s)) "SCT" else DefaultAssay(obj_s)
  }

  genes_present <- intersect(genes, rownames(obj_s[[assay_use]]))

  if (length(genes_present) < 2) {
    warning(score_name, ": fewer than 2 genes present in assay ", assay_use)
    obj_s@meta.data[[paste0(score_name, "_score")]] <- NA_real_
    return(obj_s)
  }

  old_cols <- colnames(obj_s@meta.data)

  obj_s <- AddModuleScore(
    object = obj_s,
    features = list(genes_present),
    assay = assay_use,
    name = paste0(score_name, "_score"),
    search = FALSE
  )

  new_col <- setdiff(colnames(obj_s@meta.data), old_cols)

  # AddModuleScore appends "1"; standardize the name.
  if (length(new_col) == 1) {
    obj_s@meta.data[[paste0(score_name, "_score")]] <- obj_s@meta.data[[new_col]]
    if (new_col != paste0(score_name, "_score")) {
      obj_s@meta.data[[new_col]] <- NULL
    }
  }

  obj_s
}

add_all_module_scores_hd <- function(obj_s) {
  assay_use <- if ("SCT" %in% Assays(obj_s)) "SCT" else DefaultAssay(obj_s)

  message("  Scoring modules using assay: ", assay_use)

  for (nm in names(module_gene_sets)) {
    obj_s <- score_one_module(
      obj_s,
      genes = module_gene_sets[[nm]],
      score_name = nm,
      assay_use = assay_use
    )
  }

  # Within-sample z-score for local co-enrichment.
  # This makes slopes/correlations more interpretable across HD samples.
  score_cols <- paste0(names(module_gene_sets), "_score")

  for (sc in score_cols) {
    zc <- paste0(sc, "_sample_z")
    obj_s@meta.data[[zc]] <- as.numeric(scale(obj_s@meta.data[[sc]]))
  }

  obj_s
}
extract_one_hd_sample <- function(sample_id_std, rds_path, fov) {

  message("Loading VisiumHD sample: ", sample_id_std)

  if (!file.exists(rds_path)) {
    warning("Missing RDS: ", rds_path)
    return(tibble())
  }

  sample_obj <- readRDS(rds_path)

  if (!fov %in% names(sample_obj@images)) {
    warning(
      "Requested FOV not found for ", sample_id_std,
      ": ", fov,
      ". Available images: ", paste(names(sample_obj@images), collapse = ", ")
    )
    return(tibble())
  }

  coords <- get_coords_by_fov_safe(sample_obj, fov)

  sample_obj <- add_all_module_scores_hd(sample_obj)

  score_cols <- c(
    "ECM_Remodeling_score_sample_z",
    "Vascular_Interface_score_sample_z",
    "Macrophage_Remodeling_score_sample_z",
    "NK_Cytotoxic_Regulatory_score_sample_z",
    "Trophoblast_Remodeling_score_sample_z"
  )

  sample_meta <- sample_obj@meta.data %>%
    rownames_to_column("cell") %>%
    as_tibble() %>%
    mutate(cell = as.character(cell)) %>%
    select(cell, all_of(score_cols))

  sample_level_meta <- obj@meta.data %>%
    as_tibble() %>%
    distinct(
      sample_id_std,
      study_id,
      platform_clean,
      Delivery_class,
      MNPBurden,
      MNPBurden_Delivery_class,
      mnp_model_sample
    ) %>%
    filter(sample_id_std == !!sample_id_std)

  out <- coords %>%
    left_join(sample_meta, by = "cell") %>%
    mutate(sample_id_std = sample_id_std) %>%
    left_join(sample_level_meta, by = "sample_id_std") %>%
    rename(
      ECM_Remodeling_score_use = ECM_Remodeling_score_sample_z,
      Vascular_Interface_score_use = Vascular_Interface_score_sample_z,
      Macrophage_Remodeling_score_use = Macrophage_Remodeling_score_sample_z,
      NK_Cytotoxic_Regulatory_score_use = NK_Cytotoxic_Regulatory_score_sample_z,
      Trophoblast_Remodeling_score_use = Trophoblast_Remodeling_score_sample_z
    )

  out
}

visiumhd_bins_all <- pmap_dfr(
  visiumhd_manifest %>% select(sample_id_std, rds_path, fov),
  extract_one_hd_sample
) %>%
  as_tibble()

write_csv(
  visiumhd_bins_all,
  file.path(outdir, "visiumhd_bins_all_with_module_scores.csv")
)


# ------------------------------------------------------------
# 5. Keep MNP-model VisiumHD samples
# ------------------------------------------------------------

visiumhd_bins_mnp <- visiumhd_bins_all %>%
  filter(
    mnp_model_sample == TRUE,
    platform_clean == "VisiumHD",
    !is.na(Delivery_class),
    !is.na(MNPBurden),
    !is.na(x),
    !is.na(y)
  ) %>%
  drop_na(
    Macrophage_Remodeling_score_use,
    NK_Cytotoxic_Regulatory_score_use,
    Vascular_Interface_score_use,
    Trophoblast_Remodeling_score_use
  )

write_csv(
  visiumhd_bins_mnp,
  file.path(outdir, "visiumhd_bins_mnp_input.csv")
)

cat("\nMNP VisiumHD input bins:\n")
visiumhd_bins_mnp %>%
  count(sample_id_std, study_id, Delivery_class, MNPBurden, name = "n_bins") %>%
  arrange(Delivery_class, MNPBurden, sample_id_std) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 6. Local mean of neighboring bins
# ------------------------------------------------------------

score_cols_use <- c(
  "ECM_Remodeling_score_use",
  "Vascular_Interface_score_use",
  "Macrophage_Remodeling_score_use",
  "NK_Cytotoxic_Regulatory_score_use",
  "Trophoblast_Remodeling_score_use"
)

calculate_local_neighbor_scores <- function(df_sample, k = 30) {

  sid <- unique(df_sample$sample_id_std)
  message("Calculating local module co-enrichment: ", sid)

  df_sample <- df_sample %>%
    mutate(row_index = row_number())

  xy <- as.matrix(df_sample[, c("x", "y")])

  nn <- FNN::get.knn(
    data = xy,
    k = min(k + 1, nrow(df_sample))
  )

  nn_idx_raw <- nn$nn.index

  # Remove self if included and retain k neighbors.
  nn_idx <- t(vapply(seq_len(nrow(nn_idx_raw)), function(i) {
    vals <- nn_idx_raw[i, ]
    vals <- vals[vals != i]
    vals <- vals[seq_len(min(k, length(vals)))]

    if (length(vals) < k) {
      vals <- c(vals, rep(NA_integer_, k - length(vals)))
    }

    vals
  }, integer(k)))

  local_scores <- map_dfc(score_cols_use, function(sc) {
    m <- matrix(df_sample[[sc]][nn_idx], nrow = nrow(nn_idx))
    tibble(!!paste0("local_", sc) := rowMeans(m, na.rm = TRUE))
  })

  bind_cols(df_sample, local_scores)
}

k_use <- 30

visiumhd_bins_local <- visiumhd_bins_mnp %>%
  group_by(sample_id_std) %>%
  group_split() %>%
  map_dfr(calculate_local_neighbor_scores, k = k_use)

write_csv(
  visiumhd_bins_local,
  file.path(outdir, paste0("visiumhd_bins_local_module_scores_k", k_use, ".csv"))
)

# ------------------------------------------------------------
# 7. Define local co-enrichment pairs
# ------------------------------------------------------------

coenrichment_pairs <- tribble(
  ~target_score,                            ~local_score,                                      ~pair_label,
  "Macrophage_Remodeling_score_use",        "local_Vascular_Interface_score_use",               "Macrophage remodeling vs local vascular interface",
  "NK_Cytotoxic_Regulatory_score_use",      "local_Vascular_Interface_score_use",               "NK cytotoxic/regulatory vs local vascular interface",
  "Macrophage_Remodeling_score_use",        "local_Trophoblast_Remodeling_score_use",           "Macrophage remodeling vs local trophoblast remodeling",
  "NK_Cytotoxic_Regulatory_score_use",      "local_Trophoblast_Remodeling_score_use",           "NK cytotoxic/regulatory vs local trophoblast remodeling",
  "Vascular_Interface_score_use",           "local_Trophoblast_Remodeling_score_use",           "Vascular interface vs local trophoblast remodeling",
  "Trophoblast_Remodeling_score_use",       "local_Vascular_Interface_score_use",               "Trophoblast remodeling vs local vascular interface"
)

# ------------------------------------------------------------
# 8. Sample-level co-enrichment slopes/correlations
# ------------------------------------------------------------

calculate_pair_stats_one_sample <- function(df_sample, target_score, local_score, pair_label) {

  sid <- unique(df_sample$sample_id_std)

  dat <- df_sample %>%
    select(
      sample_id_std,
      study_id,
      Delivery_class,
      MNPBurden,
      MNPBurden_Delivery_class,
      all_of(target_score),
      all_of(local_score)
    ) %>%
    drop_na()

  if (nrow(dat) < 30 ||
      sd(dat[[target_score]], na.rm = TRUE) == 0 ||
      sd(dat[[local_score]], na.rm = TRUE) == 0) {
    return(tibble(
      sample_id_std = sid,
      study_id = unique(df_sample$study_id)[1],
      Delivery_class = unique(df_sample$Delivery_class)[1],
      MNPBurden = unique(df_sample$MNPBurden)[1],
      MNPBurden_Delivery_class = unique(df_sample$MNPBurden_Delivery_class)[1],
      pair_label = pair_label,
      n_bins = nrow(dat),
      pearson_r = NA_real_,
      spearman_rho = NA_real_,
      slope = NA_real_,
      p_value_slope = NA_real_
    ))
  }

  fit <- lm(
    as.formula(paste0(target_score, " ~ ", local_score)),
    data = dat
  )

  slope_row <- broom::tidy(fit) %>%
    filter(term == local_score)

  tibble(
    sample_id_std = sid,
    study_id = unique(df_sample$study_id)[1],
    Delivery_class = unique(df_sample$Delivery_class)[1],
    MNPBurden = unique(df_sample$MNPBurden)[1],
    MNPBurden_Delivery_class = unique(df_sample$MNPBurden_Delivery_class)[1],
    pair_label = pair_label,
    n_bins = nrow(dat),
    pearson_r = suppressWarnings(cor(dat[[target_score]], dat[[local_score]], method = "pearson")),
    spearman_rho = suppressWarnings(cor(dat[[target_score]], dat[[local_score]], method = "spearman")),
    slope = slope_row$estimate[1],
    p_value_slope = slope_row$p.value[1]
  )
}

sample_level_coenrichment <- visiumhd_bins_local %>%
  group_by(sample_id_std) %>%
  group_split() %>%
  map_dfr(function(df_s) {
    pmap_dfr(
      coenrichment_pairs,
      ~ calculate_pair_stats_one_sample(
        df_sample = df_s,
        target_score = ..1,
        local_score = ..2,
        pair_label = ..3
      )
    )
  })

write_csv(
  sample_level_coenrichment,
  file.path(outdir, paste0("sample_level_visiumhd_local_module_coenrichment_k", k_use, ".csv"))
)

cat("\nSample-level local co-enrichment:\n")
sample_level_coenrichment %>%
  arrange(pair_label, Delivery_class, MNPBurden, sample_id_std) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 9. Subject-level co-enrichment
# ------------------------------------------------------------

subject_level_coenrichment <- sample_level_coenrichment %>%
  group_by(
    study_id,
    Delivery_class,
    MNPBurden,
    MNPBurden_Delivery_class,
    pair_label
  ) %>%
  summarise(
    n_samples = n(),
    n_bins = sum(n_bins, na.rm = TRUE),
    pearson_r = weighted.mean(pearson_r, w = n_bins, na.rm = TRUE),
    spearman_rho = weighted.mean(spearman_rho, w = n_bins, na.rm = TRUE),
    slope = weighted.mean(slope, w = n_bins, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  subject_level_coenrichment,
  file.path(outdir, paste0("subject_level_visiumhd_local_module_coenrichment_k", k_use, ".csv"))
)

# ------------------------------------------------------------
# 10. Within-delivery high-low contrasts
# ------------------------------------------------------------

within_delivery_coenrichment_contrasts <- subject_level_coenrichment %>%
  pivot_longer(
    cols = c(pearson_r, spearman_rho, slope),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(pair_label, metric, Delivery_class) %>%
  filter(n_distinct(MNPBurden) == 2) %>%
  summarise(
    n_subjects = n_distinct(study_id),
    mean_low = mean(value[MNPBurden == "Low"], na.rm = TRUE),
    mean_high = mean(value[MNPBurden == "High"], na.rm = TRUE),
    estimate_high_minus_low = mean_high - mean_low,
    p_value = tryCatch(
      t.test(value ~ MNPBurden)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write_csv(
  within_delivery_coenrichment_contrasts,
  file.path(outdir, paste0("within_delivery_visiumhd_local_module_coenrichment_contrasts_k", k_use, ".csv"))
)

cat("\nTop VisiumHD local co-enrichment contrasts:\n")
within_delivery_coenrichment_contrasts %>%
  arrange(p_value) %>%
  as_tibble() %>%
  print(n = Inf)

# ------------------------------------------------------------
# 11. Heatmap for slope high-low differences
# ------------------------------------------------------------

plot_df <- within_delivery_coenrichment_contrasts %>%
  filter(metric %in% c("spearman_rho", "slope")) %>%
  mutate(
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm")),
    metric = recode(
      metric,
      spearman_rho = "Spearman rho",
      slope = "Slope"
    ),
    sig_label = case_when(
      p_value < 0.05 ~ "●",
      p_value < 0.10 ~ "·",
      TRUE ~ ""
    )
  )

max_abs <- max(abs(plot_df$estimate_high_minus_low), na.rm = TRUE)

p_hd_heatmap <- ggplot(
  plot_df,
  aes(x = pair_label, y = metric, fill = estimate_high_minus_low)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3.2, color = "black") +
  facet_wrap(~ Delivery_class, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = "High − Low"
  ) +
  labs(
    x = "Local module co-enrichment pair",
    y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 7.2),
    axis.text.y = element_text(color = "black", size = 8.5),
    axis.title = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black", size = 11),
    legend.position = "right"
  )

p_hd_heatmap

ggsave(
  filename = file.path(outdir, paste0("visiumhd_local_module_coenrichment_high_low_heatmap_k", k_use, ".png")),
  plot = p_hd_heatmap,
  width = 12,
  height = 3.4,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------
# 12. Focused dotplot
# ------------------------------------------------------------

focused_plot_df <- subject_level_coenrichment %>%
  filter(
    pair_label %in% c(
      "Macrophage remodeling vs local vascular interface",
      "NK cytotoxic/regulatory vs local vascular interface",
      "Macrophage remodeling vs local trophoblast remodeling",
      "NK cytotoxic/regulatory vs local trophoblast remodeling"
    )
  ) %>%
  mutate(
    pair_label = factor(
      pair_label,
      levels = c(
        "Macrophage remodeling vs local vascular interface",
        "NK cytotoxic/regulatory vs local vascular interface",
        "Macrophage remodeling vs local trophoblast remodeling",
        "NK cytotoxic/regulatory vs local trophoblast remodeling"
      )
    ),
    MNPBurden = factor(MNPBurden, levels = c("Low", "High")),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm"))
  )

p_hd_focused <- ggplot(
  focused_plot_df,
  aes(x = MNPBurden, y = spearman_rho)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
  geom_point(
    aes(group = study_id),
    position = position_jitter(width = 0.08, height = 0),
    size = 2.0,
    alpha = 0.9
  ) +
  facet_grid(pair_label ~ Delivery_class, scales = "free_y") +
  labs(
    x = "MNP burden",
    y = "Local module co-enrichment, Spearman rho"
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

p_hd_focused

ggsave(
  filename = file.path(outdir, paste0("visiumhd_local_module_coenrichment_dotplots_k", k_use, ".png")),
  plot = p_hd_focused,
  width = 7.8,
  height = 5.8,
  units = "in",
  dpi = 600
)

cat("\nDone. Outputs written to: ", outdir, "\n")




















suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(pheatmap)
})
base_de <- "/home/ebarrozo/visiumhd_mnp/results/seurat_all_v1/seurat_all_ref_merged/de_MNPBurden_Delivery_class_bulk_sensitive_go"

term_dir <- file.path(base_de, "pb_Term_MNPBurden_limma_SCT")
preterm_dir <- file.path(base_de, "pb_Preterm_MNPBurden_limma_SCT")

list.files(term_dir, pattern = "pb_DE", full.names = FALSE)
list.files(preterm_dir, pattern = "pb_DE", full.names = FALSE)

term_de <- readr::read_tsv(
  file.path(term_dir, "pb_DE_ALL.tsv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(Delivery_class = "Term")

preterm_de <- readr::read_tsv(
  file.path(preterm_dir, "pb_DE_ALL.tsv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(Delivery_class = "Preterm")

read_pb_de_dir <- function(de_dir, delivery_class) {
  files <- list.files(
    de_dir,
    pattern = "^pb_DE_.*\\.tsv$",
    full.names = TRUE
  )

  files <- files[
    !grepl("pb_DE_ALL\\.tsv$", basename(files)) &
      !grepl("pb_DE_summary_by_celltype\\.tsv$", basename(files))
  ]

  message("Reading ", length(files), " DE files from: ", de_dir)
  print(basename(files))

  purrr::map_dfr(files, function(fp) {
    readr::read_tsv(fp, show_col_types = FALSE) %>%
      dplyr::mutate(
        source_file = basename(fp),
        Delivery_class = delivery_class
      )
  })
}

term_de <- read_pb_de_dir(term_dir, "Term")

preterm_de <- readr::read_tsv(
  file.path(preterm_dir, "pb_DE_ALL.tsv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(Delivery_class = "Preterm")

# Optional: write the merged Term table so this does not happen again
readr::write_tsv(
  term_de,
  file.path(term_dir, "pb_DE_ALL.tsv")
)

de_all <- dplyr::bind_rows(term_de, preterm_de) %>%
  dplyr::filter(!is.na(gene), !is.na(predicted.celltype), !is.na(logFC)) %>%
  dplyr::mutate(
    contrast = paste(Delivery_class, predicted.celltype, sep = " | "),
    sig_q05 = adj.P.Val < 0.05 & abs(logFC) >= 0.25,
    sig_q10 = adj.P.Val < 0.10 & abs(logFC) >= 0.25
  )

# Choose genes for Fig. 3a.
# Use q05 if you want the title to say q<0.05.
top_genes <- de_all %>%
  dplyr::filter(sig_q05) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    max_abs_logFC = max(abs(logFC), na.rm = TRUE),
    n_sig = sum(sig_q05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_sig), dplyr::desc(max_abs_logFC)) %>%
  dplyr::slice_head(n = 60) %>%
  dplyr::pull(gene)

fig3a_mat <- de_all %>%
  dplyr::filter(gene %in% top_genes) %>%
  dplyr::select(gene, contrast, logFC) %>%
  dplyr::group_by(gene, contrast) %>%
  dplyr::summarise(logFC = logFC[which.max(abs(logFC))], .groups = "drop") %>%
  tidyr::pivot_wider(names_from = contrast, values_from = logFC, values_fill = 0) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# Optional: order columns by delivery class, then cell type
fig3a_mat <- fig3a_mat[, order(colnames(fig3a_mat)), drop = FALSE]

max_abs <- quantile(abs(fig3a_mat), 0.98, na.rm = TRUE)
fig3a_mat_plot <- pmax(pmin(fig3a_mat, max_abs), -max_abs)

png("Fig3a_subject_aware_limma_logFC_heatmap.png", width = 2600, height = 2200, res = 300)
pheatmap::pheatmap(
  fig3a_mat_plot,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks = seq(-max_abs, max_abs, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize_row = 6,
  fontsize_col = 7,
  main = "Subject-aware limma MNP burden effects"
)
dev.off()








exclude_pattern <- "^(MT-|RPL|RPS|DDX3Y|KDM5D|RPS4Y1|EIF1AY|UTY|USP9Y|XIST|MALAT1|H19$|HBB|HBA|HBM|HBD|HBE|HBQ)"

de_all_clean <- de_all %>%
  dplyr::filter(
    predicted.celltype != "RCTD-lowconf",
    !grepl(exclude_pattern, gene)
  )

# prioritize interpretable genes with recurrent/significant effects
top_genes_clean <- de_all_clean %>%
  dplyr::filter(sig_q05) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    max_abs_logFC = max(abs(logFC), na.rm = TRUE),
    n_sig = sum(sig_q05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_sig), dplyr::desc(max_abs_logFC)) %>%
  dplyr::slice_head(n = 40) %>%
  dplyr::pull(gene)

fig3a_mat <- de_all_clean %>%
  dplyr::filter(gene %in% top_genes_clean) %>%
  dplyr::select(gene, Delivery_class, predicted.celltype, logFC) %>%
  dplyr::mutate(
    contrast = paste(Delivery_class, predicted.celltype, sep = " | ")
  ) %>%
  dplyr::group_by(gene, contrast) %>%
  dplyr::summarise(logFC = logFC[which.max(abs(logFC))], .groups = "drop") %>%
  tidyr::pivot_wider(names_from = contrast, values_from = logFC, values_fill = 0) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# manual column order
celltype_order <- c(
  "B Cell",
  "T Cell",
  "NK Cell",
  "Macrophage",
  "Hofbauer Cell",
  "Placental Endothelial",
  "Perivascular Cell",
  "Villous Cytotrophoblast",
  "Extravillous Trophoblast",
  "Syncytiotrophoblast"
)

desired_cols <- c(
  paste("Preterm", celltype_order, sep = " | "),
  paste("Term", celltype_order, sep = " | ")
)

desired_cols <- desired_cols[desired_cols %in% colnames(fig3a_mat)]
fig3a_mat <- fig3a_mat[, desired_cols, drop = FALSE]

max_abs <- stats::quantile(abs(fig3a_mat), 0.98, na.rm = TRUE)
fig3a_mat_plot <- pmax(pmin(fig3a_mat, max_abs), -max_abs)

fig_outdir <- file.path(base_de, "figures")
dir.create(fig_outdir, showWarnings = FALSE, recursive = TRUE)

png(
  file.path(fig_outdir, "Fig3a_subject_aware_limma_logFC_heatmap_clean.png"),
  width = 2400,
  height = 1900,
  res = 300
)

pheatmap::pheatmap(
  fig3a_mat_plot,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks = seq(-max_abs, max_abs, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col = 90,
  main = "MNP burden-associated pseudobulk effects",
  legend_breaks = c(-4, -2, 0, 2, 4),
  legend_labels = c("-4", "-2", "0", "2", "4")
)

dev.off()






library(grid)
library(gtable)

png(
  file.path(fig_outdir, "Fig3a_subject_aware_limma_logFC_heatmap_clean2.png"),
  width = 2400,
  height = 1900,
  res = 300
)

ph <- pheatmap::pheatmap(
  fig3a_mat_plot,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks = seq(-max_abs, max_abs, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col = 90,
  main = "MNP burden-associated pseudobulk effects",
  legend_breaks = c(-4, -2, 0, 2, 4),
  legend_labels = c("-4", "-2", "0", "2", "4"),
  silent = TRUE
)

# Add scale bar / legend title
legend_idx <- which(ph$gtable$layout$name == "legend")

ph$gtable <- gtable::gtable_add_rows(
  ph$gtable,
  heights = grid::unit(0.25, "in"),
  pos = ph$gtable$layout$t[legend_idx] - 1
)

ph$gtable <- gtable::gtable_add_grob(
  ph$gtable,
  grobs = grid::textGrob(
    "limma log2FC\nHigh vs Low MNP",
    gp = grid::gpar(fontsize = 8, fontface = "bold"),
    just = "center"
  ),
  t = ph$gtable$layout$t[legend_idx],
  l = ph$gtable$layout$l[legend_idx],
  b = ph$gtable$layout$t[legend_idx],
  r = ph$gtable$layout$r[legend_idx],
  name = "legend_title"
)

grid::grid.newpage()
grid::grid.draw(ph$gtable)

dev.off()



suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

# Convert matrix to long format
fig3a_long <- fig3a_mat_plot %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "contrast",
    values_to = "logFC"
  ) %>%
  tidyr::separate(
    contrast,
    into = c("Delivery_class", "predicted.celltype"),
    sep = " \\| ",
    remove = FALSE
  )

# Preserve current row/column order
gene_order <- rownames(fig3a_mat_plot)
contrast_order <- colnames(fig3a_mat_plot)

fig3a_long <- fig3a_long %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(gene_order)),
    contrast = factor(contrast, levels = contrast_order),
    Delivery_class = factor(Delivery_class, levels = c("Preterm", "Term"))
  )

p_fig3a <- ggplot(fig3a_long, aes(x = contrast, y = gene, fill = logFC)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    oob = scales::squish,
    name = "limma log2FC\nHigh vs Low MNP"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "MNP burden-associated pseudobulk effects"
  ) +
  theme_classic(base_size = 8) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", size = 6.5),
    axis.text.y = element_text(color = "black", size = 6.5),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25)
  )

p_fig3a


fig_outdir <- file.path(base_de, "figures")
dir.create(fig_outdir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(fig_outdir, "Fig3a_subject_aware_limma_logFC_heatmap_clean_ggplot.png"),
  plot = p_fig3a,
  width = 8.0,
  height = 6.2,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(fig_outdir, "Fig3a_subject_aware_limma_logFC_heatmap_clean_ggplot.pdf"),
  plot = p_fig3a,
  width = 8.0,
  height = 6.2,
  units = "in"
)

# saved here: /home/ebarrozo/visiumhd_mnp/results/seurat_all_v1/seurat_all_ref_merged/de_MNPBurden_Delivery_class_bulk_sensitive_go/figures




# ============================================================
# Supplementary analyses:
# 1) Cross-platform cell-type annotation support
# 2) Maternal-origin prediction by cell type
# 3) NK / T-cell maternal-origin and decidual-like marker support
#
# Input:
#   obj = integrated Seurat object with predicted cell types
#
# Outputs:
#   supplement_annotation_qc/
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(purrr)
})

outdir <- "supplement_annotation_qc"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

stopifnot(exists("obj"))
md <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  tibble::as_tibble()

# ------------------------------------------------------------
# 0. Standardize likely column names
# ------------------------------------------------------------

# Platform
if ("platform_clean" %in% colnames(md)) {
  md$platform_use <- as.character(md$platform_clean)
} else if ("platform" %in% colnames(md)) {
  md$platform_use <- as.character(md$platform)
} else {
  stop("Could not find platform_clean or platform in obj@meta.data.")
}

# Predicted cell type
if ("nn.predicted.celltype" %in% colnames(md)) {
  md$celltype_use <- as.character(md$nn.predicted.celltype)
} else if ("predicted.celltype" %in% colnames(md)) {
  md$celltype_use <- as.character(md$predicted.celltype)
} else {
  stop("Could not find nn.predicted.celltype or predicted.celltype in obj@meta.data.")
}

# Subject / sample IDs
if (!"study_id" %in% colnames(md)) {
  warning("study_id not found. Using sample_id_std as subject surrogate.")
  md$study_id <- md$sample_id_std
}

if (!"sample_id_std" %in% colnames(md)) {
  warning("sample_id_std not found. Using study_id as sample surrogate.")
  md$sample_id_std <- md$study_id
}

# MNP metadata if present
if (!"Delivery_class" %in% colnames(md)) md$Delivery_class <- NA_character_
if (!"MNPBurden" %in% colnames(md)) md$MNPBurden <- NA_character_

# ------------------------------------------------------------
# 1. Identify maternal-origin column
# ------------------------------------------------------------

origin_candidates <- grep(
  "maternal|fetal|origin|XIST|Y_score|Yscore|sex",
  colnames(md),
  value = TRUE,
  ignore.case = TRUE
)

cat("\nCandidate origin columns:\n")
print(origin_candidates)

# Edit this if needed after inspecting candidates.
# Common likely values might be:
#   predicted_origin
#   predicted.maternal.origin
#   maternal_origin
#   origin_prediction
#   predicted_parental_origin
origin_col <- NA_character_

possible_origin_cols <- c(
  "predicted_origin",
  "predicted.maternal.origin",
  "predicted_maternal_origin",
  "maternal_origin",
  "origin_prediction",
  "predicted_parental_origin",
  "predicted.origin",
  "origin", "origin_rule_based"
)

origin_col <- possible_origin_cols[possible_origin_cols %in% colnames(md)][1]

if (is.na(origin_col)) {
  message("No obvious maternal-origin categorical column found.")
  message("A fallback will be created if XIST_score/Y_score columns exist.")
}

# ------------------------------------------------------------
# 2. Create maternal-origin flag
# ------------------------------------------------------------

if (!is.na(origin_col)) {

  md$origin_use <- as.character(md[[origin_col]])

  md$maternal_origin_pred <- case_when(
    str_detect(tolower(md$origin_use), "maternal|mother|mat") ~ TRUE,
    str_detect(tolower(md$origin_use), "fetal|foetal|fetus|placental|fet") ~ FALSE,
    TRUE ~ NA
  )

} else if (all(c("XIST_score", "Y_score") %in% colnames(md))) {

  # Fallback heuristic only. Adjust thresholds if needed.
  # This is not a definitive origin classifier.
  md <- md %>%
    mutate(
      XIST_score_z = as.numeric(scale(XIST_score)),
      Y_score_z = as.numeric(scale(Y_score)),
      maternal_origin_pred = case_when(
        XIST_score_z > 1 & Y_score_z < 0 ~ TRUE,
        Y_score_z > 1 & XIST_score_z < 0 ~ FALSE,
        TRUE ~ NA
      ),
      origin_use = case_when(
        maternal_origin_pred == TRUE ~ "Predicted maternal",
        maternal_origin_pred == FALSE ~ "Predicted fetal/placental",
        TRUE ~ "Unassigned"
      )
    )

  warning(
    "Maternal origin was inferred using XIST/Y-score heuristic. ",
    "Use only as supportive annotation unless validated."
  )

} else {

  md$origin_use <- NA_character_
  md$maternal_origin_pred <- NA

  warning(
    "No maternal-origin column and no XIST_score/Y_score fallback found. ",
    "Maternal-origin outputs will be empty."
  )
}

readr::write_csv(
  tibble(candidate_origin_columns = origin_candidates),
  file.path(outdir, "candidate_maternal_origin_columns.csv")
)

# ------------------------------------------------------------
# 3. Cross-platform predicted cell-type composition
# ------------------------------------------------------------
# ------------------------------------------------------------
# Refine platform labels: VisiumHD, VisiumV1, Xenium500, Xenium5k
# ------------------------------------------------------------

md <- md %>%
  dplyr::mutate(
    sample_id_std = as.character(sample_id_std),
    platform_use = as.character(platform_use)
  )

# First inspect sample IDs by platform
md %>%
  dplyr::distinct(sample_id_std, platform_use) %>%
  dplyr::arrange(platform_use, sample_id_std) %>%
  print(n = Inf)

## meta_model_ready.csv
# Edit this after checking sample_id_std.
# These are placeholders. Use your true Xenium 500 sample_id_std values.
xenium500_ids <- c(
  "HT1_xen",
  "HP13b_xen",
  "NC4_xen",
  "HT2a_xen"
)

md <- md %>%
  dplyr::mutate(
    platform_detail = dplyr::case_when(
      platform_use %in% c("VisiumHD", "Visium_HD", "Visium HD") ~ "VisiumHD",
      platform_use %in% c("VisiumV1", "Visium_v1", "Visium v1", "VisiumV1") ~ "VisiumV1",
      platform_use %in% c("Xenium", "Xenium500", "Xenium5k") &
        sample_id_std %in% xenium500_ids ~ "Xenium500",
      platform_use %in% c("Xenium", "Xenium500", "Xenium5k") &
        !sample_id_std %in% xenium500_ids ~ "Xenium5k",
      TRUE ~ platform_use
    ),
    platform_detail = factor(
      platform_detail,
      levels = c("VisiumV1", "Xenium500", "Xenium5k", "VisiumHD")
    )
  )

md %>%
  dplyr::distinct(sample_id_std, platform_use, platform_detail) %>%
  dplyr::arrange(platform_detail, sample_id_std) %>%
  print(n = Inf)

# ------------------------------------------------------------
# Platform-detail cell-type composition
# ------------------------------------------------------------

platform_celltype_counts_detail <- md %>%
  dplyr::filter(!is.na(platform_detail), !is.na(celltype_use)) %>%
  dplyr::count(platform_detail, celltype_use, name = "n_spatial_obs") %>%
  dplyr::group_by(platform_detail) %>%
  dplyr::mutate(
    platform_total = sum(n_spatial_obs),
    frac_of_platform = n_spatial_obs / platform_total
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(platform_detail, dplyr::desc(frac_of_platform))

readr::write_csv(
  platform_celltype_counts_detail,
  file.path(outdir, "platform_detail_celltype_counts_and_fractions.csv")
)

p_platform_comp_detail <- platform_celltype_counts_detail %>%
  dplyr::mutate(
    celltype_use = forcats::fct_reorder(celltype_use, frac_of_platform, .fun = max)
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = celltype_use, y = frac_of_platform)) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::facet_wrap(~ platform_detail, ncol = 1, scales = "free_y") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = "Fraction of spatial observations",
    title = "Predicted cell-type composition by platform and assay version"
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    text = ggplot2::element_text(family = "Arial"),
    axis.text = ggplot2::element_text(color = "black"),
    axis.title = ggplot2::element_text(color = "black"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", color = "black"),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
  )

ggplot2::ggsave(
  file.path(outdir, "Supp_platform_detail_celltype_composition.png"),
  p_platform_comp_detail,
  width = 6.8,
  height = 8.2,
  dpi = 600
)



# ------------------------------------------------------------
# 4. Cross-platform correlation of cell-type composition
# ------------------------------------------------------------
# This does not prove one-to-one cell equivalence, but supports
# broad cross-platform recovery of major compartments.

platform_comp_wide <- platform_celltype_counts %>%
  dplyr::select(platform_use, celltype_use, frac_of_platform) %>%
  pivot_wider(
    names_from = platform_use,
    values_from = frac_of_platform,
    values_fill = 0
  )

platform_cols <- setdiff(colnames(platform_comp_wide), "celltype_use")

platform_cor <- cor(
  as.matrix(platform_comp_wide[, platform_cols, drop = FALSE]),
  method = "spearman",
  use = "pairwise.complete.obs"
)

platform_cor_df <- as.data.frame(platform_cor) %>%
  rownames_to_column("platform_1") %>%
  pivot_longer(
    cols = -platform_1,
    names_to = "platform_2",
    values_to = "spearman_rho"
  )

readr::write_csv(
  platform_cor_df,
  file.path(outdir, "platform_celltype_fraction_spearman_correlations.csv")
)

p_platform_cor <- platform_cor_df %>%
  ggplot(aes(x = platform_1, y = platform_2, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 3) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Cross-platform similarity of predicted cell-type fractions"
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.ticks = element_blank()
  )

ggsave(
  file.path(outdir, "Supp_platform_celltype_fraction_correlations.png"),
  p_platform_cor,
  width = 4.8,
  height = 4.2,
  dpi = 600
)


# ------------------------------------------------------------
# Platform-detail composition similarity
# ------------------------------------------------------------

platform_comp_wide_detail <- platform_celltype_counts_detail %>%
  dplyr::select(platform_detail, celltype_use, frac_of_platform) %>%
  tidyr::pivot_wider(
    names_from = platform_detail,
    values_from = frac_of_platform,
    values_fill = 0
  )

platform_cols_detail <- setdiff(colnames(platform_comp_wide_detail), "celltype_use")

platform_cor_detail <- stats::cor(
  as.matrix(platform_comp_wide_detail[, platform_cols_detail, drop = FALSE]),
  method = "spearman",
  use = "pairwise.complete.obs"
)

platform_cor_detail_df <- as.data.frame(platform_cor_detail) %>%
  tibble::rownames_to_column("platform_1") %>%
  tidyr::pivot_longer(
    cols = -platform_1,
    names_to = "platform_2",
    values_to = "spearman_rho"
  )

readr::write_csv(
  platform_cor_detail_df,
  file.path(outdir, "platform_detail_celltype_fraction_spearman_correlations.csv")
)

p_platform_cor_detail <- platform_cor_detail_df %>%
  ggplot2::ggplot(ggplot2::aes(x = platform_1, y = platform_2, fill = spearman_rho)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", spearman_rho)), size = 3) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    title = "Cross-platform similarity of predicted cell-type fractions"
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    text = ggplot2::element_text(family = "Arial"),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, color = "black"),
    axis.text.y = ggplot2::element_text(color = "black"),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    axis.ticks = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(outdir, "Supp_platform_detail_celltype_fraction_correlations.png"),
  p_platform_cor_detail,
  width = 5.2,
  height = 4.6,
  dpi = 600
)

# ------------------------------------------------------------
# 5. Maternal-origin assignment by predicted cell type
# ------------------------------------------------------------

origin_by_celltype <- md %>%
  dplyr::filter(!is.na(celltype_use)) %>%
  dplyr::mutate(
    origin_group = dplyr::case_when(
      maternal_origin_pred == TRUE ~ "Predicted maternal origin",
      maternal_origin_pred == FALSE ~ "Predicted fetal/placental origin",
      TRUE ~ "Unassigned"
    )
  ) %>%
  dplyr::count(celltype_use, origin_group, name = "n_spatial_obs") %>%
  dplyr::group_by(celltype_use) %>%
  dplyr::mutate(
    celltype_total = sum(n_spatial_obs),
    frac = n_spatial_obs / celltype_total
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(celltype_use, dplyr::desc(frac))

readr::write_csv(
  origin_by_celltype,
  file.path(outdir, "maternal_origin_by_predicted_celltype.csv")
)

# Wide summary for text/table
origin_by_celltype_wide <- origin_by_celltype %>%
  select(celltype_use, origin_group, n_spatial_obs, frac) %>%
  pivot_wider(
    names_from = origin_group,
    values_from = c(n_spatial_obs, frac),
    values_fill = 0
  )

readr::write_csv(
  origin_by_celltype_wide,
  file.path(outdir, "maternal_origin_by_predicted_celltype_wide.csv")
)

# Plot maternal origin fraction by cell type
p_origin_celltype <- origin_by_celltype %>%
  filter(origin_group != "Unassigned") %>%
  mutate(
    celltype_use = fct_reorder(celltype_use, frac, .fun = max)
  ) %>%
  ggplot(aes(x = celltype_use, y = frac, fill = origin_group)) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Fraction of spatial observations",
    title = "Predicted maternal-origin signal by cell type",
    fill = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(
  file.path(outdir, "Supp_maternal_origin_fraction_by_celltype.png"),
  p_origin_celltype,
  width = 6.8,
  height = 5.8,
  dpi = 600
)

# ------------------------------------------------------------
# 6. Focused NK/T maternal-origin summary
# ------------------------------------------------------------

immune_origin_focus <- origin_by_celltype %>%
  filter(str_detect(celltype_use, regex("NK|T Cell|T-cell|T_", ignore_case = TRUE))) %>%
  arrange(celltype_use, origin_group)

readr::write_csv(
  immune_origin_focus,
  file.path(outdir, "maternal_origin_NK_T_focus.csv")
)

print(immune_origin_focus, n = Inf)

# Optional: subject-level NK maternal-origin fraction
nk_subject_origin <- md %>%
  dplyr::filter(stringr::str_detect(celltype_use, stringr::regex("NK", ignore_case = TRUE))) %>%
  dplyr::mutate(
    origin_group = dplyr::case_when(
      maternal_origin_pred == TRUE ~ "Predicted maternal origin",
      maternal_origin_pred == FALSE ~ "Predicted fetal/placental origin",
      TRUE ~ "Unassigned"
    )
  ) %>%
  dplyr::count(study_id, Delivery_class, MNPBurden, origin_group, name = "n_nk_obs") %>%
  dplyr::group_by(study_id, Delivery_class, MNPBurden) %>%
  dplyr::mutate(
    n_nk_total = sum(n_nk_obs),
    frac = n_nk_obs / n_nk_total
  ) %>%
  dplyr::ungroup()

readr::write_csv(
  nk_subject_origin,
  file.path(outdir, "NK_subject_level_maternal_origin_fractions.csv")
)

p_nk_origin_subject <- nk_subject_origin %>%
  filter(origin_group != "Unassigned") %>%
  ggplot(aes(x = Delivery_class, y = frac, fill = origin_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0),
    size = 1.8,
    alpha = 0.8
  ) +
  facet_wrap(~ origin_group, nrow = 1) +
  labs(
    x = NULL,
    y = "Subject-level fraction of NK observations",
    title = "Predicted maternal-origin signal among NK-cell observations"
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black"),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  file.path(outdir, "Supp_NK_subject_level_maternal_origin_fraction.png"),
  p_nk_origin_subject,
  width = 6.2,
  height = 3.8,
  dpi = 600
)

nk_scores_subject_mnp <- nk_scores_subject %>%
  dplyr::filter(!is.na(MNPBurden), MNPBurden %in% c("Low", "High")) %>%
  dplyr::mutate(
    MNPBurden = factor(MNPBurden, levels = c("Low", "High")),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm"))
  )

p_nk_scores_mnp <- ggplot2::ggplot(
  nk_scores_subject_mnp,
  ggplot2::aes(x = MNPBurden, y = score)
) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
  ggplot2::geom_point(
    ggplot2::aes(group = study_id),
    position = ggplot2::position_jitter(width = 0.08, height = 0),
    size = 1.8,
    alpha = 0.85
  ) +
  ggplot2::facet_grid(program ~ Delivery_class, scales = "free_y") +
  ggplot2::labs(
    x = "MNP burden",
    y = "Mean NK marker score",
    title = "NK-cell marker programs by MNP burden and delivery class"
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    text = ggplot2::element_text(family = "Arial"),
    axis.text = ggplot2::element_text(color = "black"),
    axis.title = ggplot2::element_text(color = "black"),
    strip.background = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_text(angle = 0, hjust = 0, color = "black", size = 7.5),
    strip.text.x = ggplot2::element_text(face = "bold", color = "black"),
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
  )

ggplot2::ggsave(
  file.path(outdir, "Supp_NK_marker_program_scores_by_MNP_delivery_MNPsamples_only.png"),
  p_nk_scores_mnp,
  width = 7.2,
  height = 5.8,
  dpi = 600
)

# ------------------------------------------------------------
# 7. Optional NK decidual-like marker scoring
# ------------------------------------------------------------
# These are broad supportive markers. Interpret cautiously.
# Adjust genes based on your platform/reference.

nk_marker_sets <- list(
  Decidual_NK_like = c("NCAM1", "KLRB1", "KLRD1", "KLRC1", "XCL1", "XCL2", "AREG", "ITGA1", "CXCR4"),
  Cytotoxic_NK = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "FGFBP2", "KLRF1"),
  Tissue_resident_NK = c("ITGA1", "CXCR6", "ZNF683", "CD69", "KLRB1"),
  Interferon_response = c("STAT1", "IRF1", "ISG15", "IFIT1", "IFIT2", "IFIT3", "CXCL10")
)

# ------------------------------------------------------------
# Robust NK marker-score block
# ------------------------------------------------------------

# Make sure celltype_use and cell are plain character vectors
md <- md %>%
  dplyr::mutate(
    cell = as.character(cell),
    celltype_use = as.character(celltype_use)
  )

nk_cells <- md %>%
  dplyr::filter(stringr::str_detect(celltype_use, stringr::regex("NK", ignore_case = TRUE))) %>%
  dplyr::pull(cell) %>%
  unlist(use.names = FALSE) %>%
  as.character() %>%
  unique()

# Keep only cells actually present in the Seurat object
nk_cells <- intersect(nk_cells, colnames(obj))

length(nk_cells)

if (length(nk_cells) > 50) {

  nk_obj <- subset(obj, cells = nk_cells)

  assay_use <- if ("SCT" %in% names(nk_obj@assays)) "SCT" else DefaultAssay(nk_obj)
  DefaultAssay(nk_obj) <- assay_use

  message("Using assay: ", assay_use)
  message("NK observations: ", ncol(nk_obj))

  for (nm in names(nk_marker_sets)) {

    genes_present <- intersect(
      as.character(nk_marker_sets[[nm]]),
      as.character(rownames(nk_obj[[assay_use]]))
    )

    message(nm, ": ", length(genes_present), " genes present: ", paste(genes_present, collapse = ", "))

    if (length(genes_present) >= 2) {

      old_cols <- colnames(nk_obj@meta.data)

      nk_obj <- AddModuleScore(
        object = nk_obj,
        features = list(genes_present),
        assay = assay_use,
        name = paste0(nm, "_score"),
        search = FALSE
      )

      new_cols <- setdiff(colnames(nk_obj@meta.data), old_cols)

      # AddModuleScore usually creates e.g. Decidual_NK_like_score1
      if (length(new_cols) >= 1) {
        nk_obj@meta.data[[paste0(nm, "_score")]] <- nk_obj@meta.data[[new_cols[1]]]
      }

    } else {
      warning("Too few genes present for NK marker set: ", nm)
    }
  }

  nk_score_cols <- paste0(names(nk_marker_sets), "_score")
  nk_score_cols <- nk_score_cols[nk_score_cols %in% colnames(nk_obj@meta.data)]

  message("NK score columns found: ", paste(nk_score_cols, collapse = ", "))

  if (length(nk_score_cols) == 0) {
    stop("No NK marker-score columns were created.")
  }

  nk_meta <- nk_obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      cell = as.character(cell),
      study_id = as.character(study_id),
      Delivery_class = as.character(Delivery_class),
      MNPBurden = as.character(MNPBurden)
    )

  # Add platform/celltype columns safely
  if ("platform_clean" %in% colnames(nk_meta)) {
    nk_meta$platform_use <- as.character(nk_meta$platform_clean)
  } else if ("platform" %in% colnames(nk_meta)) {
    nk_meta$platform_use <- as.character(nk_meta$platform)
  } else {
    nk_meta$platform_use <- NA_character_
  }

  if ("nn.predicted.celltype" %in% colnames(nk_meta)) {
    nk_meta$celltype_use <- as.character(nk_meta[["nn.predicted.celltype"]])
  } else if ("predicted.celltype" %in% colnames(nk_meta)) {
    nk_meta$celltype_use <- as.character(nk_meta$predicted.celltype)
  } else {
    nk_meta$celltype_use <- NA_character_
  }

  nk_scores_subject <- nk_meta %>%
    dplyr::group_by(study_id, Delivery_class, MNPBurden) %>%
    dplyr::summarise(
      n_nk_obs = dplyr::n(),
      dplyr::across(dplyr::all_of(nk_score_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(nk_score_cols),
      names_to = "program",
      values_to = "score"
    ) %>%
    dplyr::mutate(program = stringr::str_remove(program, "_score$"))

  readr::write_csv(
    nk_scores_subject,
    file.path(outdir, "NK_subject_level_decidual_like_marker_scores.csv")
  )

  p_nk_scores <- ggplot(
    nk_scores_subject,
    aes(x = MNPBurden, y = score)
  ) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
    geom_point(
      aes(group = study_id),
      position = position_jitter(width = 0.08, height = 0),
      size = 1.8,
      alpha = 0.85
    ) +
    facet_grid(program ~ Delivery_class, scales = "free_y") +
    labs(
      x = "MNP burden",
      y = "Mean NK marker score",
      title = "NK-cell marker programs by MNP burden and delivery class"
    ) +
    theme_classic(base_size = 9) +
    theme(
      text = element_text(family = "Arial"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text.y = element_text(angle = 0, hjust = 0, color = "black", size = 7.5),
      strip.text.x = element_text(face = "bold", color = "black"),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  ggsave(
    file.path(outdir, "Supp_NK_marker_program_scores_by_MNP_delivery.png"),
    p_nk_scores,
    width = 7.2,
    height = 5.8,
    dpi = 600
  )

} else {
  warning("Too few NK observations for NK marker scoring: ", length(nk_cells))
}
# ------------------------------------------------------------
# 8. Write a compact README / analysis log
# ------------------------------------------------------------

readr::write_lines(
  c(
    "Supplementary annotation QC outputs",
    "===================================",
    paste0("Date: ", Sys.Date()),
    "",
    "Key outputs:",
    "1. platform_celltype_counts_and_fractions.csv",
    "2. subject_platform_celltype_fractions.csv",
    "3. platform_celltype_fraction_spearman_correlations.csv",
    "4. maternal_origin_by_predicted_celltype.csv",
    "5. maternal_origin_by_predicted_celltype_wide.csv",
    "6. maternal_origin_NK_T_focus.csv",
    "7. NK_subject_level_maternal_origin_fractions.csv",
    "8. Optional NK_subject_level_decidual_like_marker_scores.csv",
    "",
    "Interpretation notes:",
    "- Cross-platform analyses evaluate recovery of broad predicted compartments, not one-to-one cell equivalence.",
    "- Maternal-origin prediction should be interpreted as supportive annotation unless independently validated.",
    "- NK decidual-like marker scores are supportive and should not be used alone to claim definitive dNK identity."
  ),
  file.path(outdir, "README_supplement_annotation_qc.txt")
)

message("Done. Outputs written to: ", normalizePath(outdir))

# ------------------------------------------------------------
# NK marker program stats: within-delivery High vs Low MNP
# ------------------------------------------------------------

nk_scores_subject_mnp <- nk_scores_subject %>%
  dplyr::filter(!is.na(MNPBurden), MNPBurden %in% c("Low", "High")) %>%
  dplyr::mutate(
    MNPBurden = factor(MNPBurden, levels = c("Low", "High")),
    Delivery_class = factor(Delivery_class, levels = c("Term", "Preterm"))
  )

nk_program_stats <- nk_scores_subject_mnp %>%
  dplyr::group_by(program, Delivery_class) %>%
  dplyr::summarise(
    n_subjects = dplyr::n_distinct(study_id),
    n_low = dplyr::n_distinct(study_id[MNPBurden == "Low"]),
    n_high = dplyr::n_distinct(study_id[MNPBurden == "High"]),
    mean_low = mean(score[MNPBurden == "Low"], na.rm = TRUE),
    mean_high = mean(score[MNPBurden == "High"], na.rm = TRUE),
    median_low = median(score[MNPBurden == "Low"], na.rm = TRUE),
    median_high = median(score[MNPBurden == "High"], na.rm = TRUE),
    estimate_high_minus_low = mean_high - mean_low,
    p_value_ttest = ifelse(
      n_low >= 2 & n_high >= 2,
      tryCatch(t.test(score ~ MNPBurden)$p.value, error = function(e) NA_real_),
      NA_real_
    ),
    p_value_wilcox = ifelse(
      n_low >= 2 & n_high >= 2,
      tryCatch(wilcox.test(score ~ MNPBurden, exact = FALSE)$p.value, error = function(e) NA_real_),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::group_by(Delivery_class) %>%
  dplyr::mutate(
    q_value_ttest = p.adjust(p_value_ttest, method = "BH"),
    q_value_wilcox = p.adjust(p_value_wilcox, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Delivery_class, p_value_ttest)

nk_program_stats

readr::write_csv(
  nk_program_stats,
  file.path(outdir, "NK_marker_program_within_delivery_MNP_stats.csv")
)

p_nk_scores_mnp <- ggplot2::ggplot(
  nk_scores_subject_mnp,
  ggplot2::aes(x = MNPBurden, y = score)
) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) +
  ggplot2::geom_point(
    ggplot2::aes(group = study_id),
    position = ggplot2::position_jitter(width = 0.08, height = 0),
    size = 1.8,
    alpha = 0.85
  ) +
  ggplot2::facet_grid(program ~ Delivery_class, scales = "free_y") +
  ggplot2::labs(
    x = "MNP burden",
    y = "Mean NK marker score",
    title = "NK-cell marker programs by MNP burden and delivery class"
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    text = ggplot2::element_text(family = "Arial"),
    axis.text = ggplot2::element_text(color = "black"),
    axis.title = ggplot2::element_text(color = "black"),
    strip.background = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_text(angle = 0, hjust = 0, color = "black", size = 7.5),
    strip.text.x = ggplot2::element_text(face = "bold", color = "black"),
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
  )

ggplot2::ggsave(
  file.path(outdir, "Supp_NK_marker_program_scores_by_MNP_delivery_MNPsamples_only.png"),
  p_nk_scores_mnp,
  width = 7.2,
  height = 5.8,
  dpi = 600
)




# ============================================================
# Publication-ready CellChat alluvial A/B panels
# For Supp Fig. 5a/b
# Uses existing CellChat objects and functions:
#   build_alluvial_tables()
#   ct_colors
#   pathways_keep
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggalluvial)
  library(CellChat)
  library(stringr)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

base_cc <- "/home/ebarrozo/visiumhd_mnp/results/audra/MNP_interactomes_CellChat_v12_full_min500_png/CellChat/pooled"

cc_preterm_high <- readRDS(file.path(base_cc, "Preterm_High", "cellchat.rds"))
cc_term_high    <- readRDS(file.path(base_cc, "Term_High", "cellchat.rds"))

# Output directory
out_ab_pub <- "/home/ebarrozo/visiumhd_mnp/results/audra/MNP_interactomes_CellChat_v12_full_min500_png/Figure_addons/alluvial_AB_core12_ordered_publication"
dir.create(out_ab_pub, showWarnings = FALSE, recursive = TRUE)

# Core pathway order
pathways_keep <- c(
  "COLLAGEN", "LAMININ", "FN1", "THBS", "TENASCIN", "SPP1",
  "ICAM", "VCAM",
  "CCL", "CXCL",
  "IL1", "TNF"
)

# ------------------------------------------------------------
# Cleaner labels
# ------------------------------------------------------------

clean_ct_label <- function(x) {
  x %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_replace("Villous Cytotrophoblast", "Villous CTB") %>%
    stringr::str_replace("Extravillous Trophoblast", "EVT") %>%
    stringr::str_replace("Syncytiotrophoblast", "STB") %>%
    stringr::str_replace("Placental Endothelial", "Endothelial") %>%
    stringr::str_replace("Perivascular Cell", "Perivascular") %>%
    stringr::str_replace("Hofbauer Cell", "Hofbauer") %>%
    stringr::str_replace("Decidual NK Cell", "Decidual NK")
}

# Optional: pathway family labels
pathway_family <- tibble::tibble(
  pathway = pathways_keep,
  family = dplyr::case_when(
    pathway %in% c("COLLAGEN", "LAMININ", "FN1", "THBS", "TENASCIN", "SPP1") ~ "ECM / adhesion",
    pathway %in% c("ICAM", "VCAM") ~ "cell adhesion",
    pathway %in% c("CCL", "CXCL") ~ "chemokine",
    pathway %in% c("IL1", "TNF") ~ "cytokine",
    TRUE ~ "other"
  )
)

# ------------------------------------------------------------
# Publication-ready plotting functions
# ------------------------------------------------------------

plot_alluvial_A_pub <- function(dfA, title = NULL, min_strength_frac = 0.00) {

  df <- dfA %>%
    dplyr::filter(pathway %in% pathways_keep) %>%
    dplyr::mutate(
      sender_raw = as.character(sender),
      sender = clean_ct_label(as.character(sender)),
      pathway = factor(pathway, levels = pathways_keep)
    )

  if (min_strength_frac > 0) {
    cutoff <- sum(df$strength, na.rm = TRUE) * min_strength_frac
    df <- df %>% dplyr::filter(strength >= cutoff)
  }

  # order sender by total strength
  sender_order <- df %>%
    dplyr::group_by(sender) %>%
    dplyr::summarise(total = sum(strength, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total)) %>%
    dplyr::pull(sender)

  df <- df %>%
    dplyr::mutate(sender = factor(sender, levels = rev(sender_order)))

  # color map with cleaned names
  ct_colors_clean <- ct_colors
  names(ct_colors_clean) <- clean_ct_label(names(ct_colors_clean))

  ggplot(df, aes(axis1 = sender, axis2 = pathway, y = strength)) +
    ggalluvial::geom_alluvium(aes(fill = sender), width = 0.18, alpha = 0.86) +
    ggalluvial::geom_stratum(width = 0.18, color = "grey25", fill = "grey95", linewidth = 0.25) +
    ggplot2::geom_text(
  stat = "stratum",
  aes(label = ifelse(after_stat(count) > 0.004, after_stat(stratum), "")),
  size = 4.0,
  lineheight = 0.9
) +
    ggplot2::scale_fill_manual(values = ct_colors_clean, drop = FALSE) +
    ggplot2::scale_x_discrete(
      limits = c("Sender cell", "Signaling pathway"),
      expand = c(0.08, 0.08)
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "Interaction strength"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold", color = "black", size = 12),
      axis.text.y = ggplot2::element_text(color = "black", size = 10),
      axis.title.y = ggplot2::element_text(color = "black", size = 11),
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
      plot.margin = ggplot2::margin(10, 20, 10, 20)
    )
}

plot_alluvial_B_pub <- function(dfB, title = NULL, min_strength_frac = 0.00) {

  df <- dfB %>%
    dplyr::filter(pathway %in% pathways_keep) %>%
    dplyr::mutate(
      receiver_raw = as.character(receiver),
      receiver = clean_ct_label(as.character(receiver)),
      pathway = factor(pathway, levels = pathways_keep)
    )

  if (min_strength_frac > 0) {
    cutoff <- sum(df$strength, na.rm = TRUE) * min_strength_frac
    df <- df %>% dplyr::filter(strength >= cutoff)
  }

  receiver_order <- df %>%
    dplyr::group_by(receiver) %>%
    dplyr::summarise(total = sum(strength, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total)) %>%
    dplyr::pull(receiver)

  df <- df %>%
    dplyr::mutate(receiver = factor(receiver, levels = rev(receiver_order)))

  ct_colors_clean <- ct_colors
  names(ct_colors_clean) <- clean_ct_label(names(ct_colors_clean))

  ggplot(df, aes(axis1 = pathway, axis2 = receiver, y = strength)) +
    ggalluvial::geom_alluvium(aes(fill = receiver), width = 0.18, alpha = 0.86) +
    ggalluvial::geom_stratum(width = 0.18, color = "grey25", fill = "grey95", linewidth = 0.25) +
    ggplot2::geom_text(
  stat = "stratum",
  aes(label = ifelse(after_stat(count) > 0.004, after_stat(stratum), "")),
  size = 4.0,
  lineheight = 0.9
) +
    ggplot2::scale_fill_manual(values = ct_colors_clean, drop = FALSE) +
    ggplot2::scale_x_discrete(
      limits = c("Signaling pathway", "Receiver cell"),
      expand = c(0.08, 0.08)
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "Interaction strength"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold", color = "black", size = 12),
      axis.text.y = ggplot2::element_text(color = "black", size = 10),
      axis.title.y = ggplot2::element_text(color = "black", size = 11),
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
      plot.margin = ggplot2::margin(10, 20, 10, 20)
    )
}

make_cellrep_alluvial_AB_pub <- function(
  cc,
  outdir,
  condition_label,
  title_prefix,
  pathways = pathways_keep,
  width_in = 13,
  height_in = 7.5,
  min_strength_frac = 0.00
) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  tabs <- build_alluvial_tables(cc, pathways)

  readr::write_csv(
    tabs$A,
    file.path(outdir, paste0("alluvial_A_sender_by_pathway_", condition_label, ".csv"))
  )
  readr::write_csv(
    tabs$B,
    file.path(outdir, paste0("alluvial_B_pathway_by_receiver_", condition_label, ".csv"))
  )
  writeLines(
    pathways,
    file.path(outdir, paste0("pathways_used_", condition_label, ".txt"))
  )

  pA <- plot_alluvial_A_pub(
    tabs$A,
    title = paste0(title_prefix, ": sender contributions to signaling pathways"),
    min_strength_frac = min_strength_frac
  )

  pB <- plot_alluvial_B_pub(
    tabs$B,
    title = paste0(title_prefix, ": signaling pathways to receiver cells"),
    min_strength_frac = min_strength_frac
  )

  ggplot2::ggsave(
    file.path(outdir, paste0("SuppFig5_", condition_label, "_A_sender_by_pathway.png")),
    pA,
    width = width_in,
    height = height_in,
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )

  ggplot2::ggsave(
    file.path(outdir, paste0("SuppFig5_", condition_label, "_A_sender_by_pathway.pdf")),
    pA,
    width = width_in,
    height = height_in,
    bg = "white",
    limitsize = FALSE
  )

  ggplot2::ggsave(
    file.path(outdir, paste0("SuppFig5_", condition_label, "_B_pathway_by_receiver.png")),
    pB,
    width = width_in,
    height = height_in,
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )

  ggplot2::ggsave(
    file.path(outdir, paste0("SuppFig5_", condition_label, "_B_pathway_by_receiver.pdf")),
    pB,
    width = width_in,
    height = height_in,
    bg = "white",
    limitsize = FALSE
  )

  invisible(list(A = pA, B = pB, tabs = tabs))
}

# ------------------------------------------------------------
# RUN: panels matching current Supp Fig. 5 
# ------------------------------------------------------------

preterm_high_panels <- make_cellrep_alluvial_AB_pub(
  cc_preterm_high,
  outdir = out_ab_pub,
  condition_label = "Preterm_High",
  title_prefix = "Preterm high-MNP",
  pathways = pathways_keep,
  width_in = 13,
  height_in = 7.5
)

term_high_panels <- make_cellrep_alluvial_AB_pub(
  cc_term_high,
  outdir = out_ab_pub,
  condition_label = "Term_High",
  title_prefix = "Term high-MNP",
  pathways = pathways_keep,
  width_in = 13,
  height_in = 7.5
)



suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(ggalluvial)
    library(readr)
})

# ============================================================
# CONFIG
# ============================================================

base_cc <- "/home/ebarrozo/visiumhd_mnp/results/audra/MNP_interactomes_CellChat_v12_full_min500_png/CellChat/pooled"

cc_preterm_high <- readRDS(file.path(base_cc, "Preterm_High", "cellchat.rds"))
cc_term_high    <- readRDS(file.path(base_cc, "Term_High", "cellchat.rds"))

out_ab_pub <- "/home/ebarrozo/visiumhd_mnp/results/audra/MNP_interactomes_CellChat_v12_full_min500_png/Figure_addons/alluvial_AB_core12_ordered_publication_final"
dir.create(out_ab_pub, showWarnings = FALSE, recursive = TRUE)

pathways_keep <- c(
    "COLLAGEN", "LAMININ", "FN1", "THBS",
    "TENASCIN", "SPP1", "ICAM", "VCAM",
    "CCL", "CXCL", "IL1", "TNF"
)

major_pathways <- c(
    "COLLAGEN", "LAMININ", "FN1", "THBS",
    "TENASCIN", "SPP1", "ICAM", "VCAM"
)

other_label <- "Other core pathways"

# ============================================================
# CLEAN LABELS
# ============================================================

clean_cell_label <- function(x) {
    x <- as.character(x)
    x <- str_replace_all(x, "_", " ")
    
    dplyr::case_when(
        x %in% c("Syncytiotrophoblast", "STB") ~ "STB",
        x %in% c("Extravillous Trophoblast", "EVT") ~ "EVT",
        x %in% c("Villous Cytotrophoblast", "Villous CTB") ~ "Villous CTB",
        x %in% c("Placental Endothelial") ~ "Endothelial",
        x %in% c("Perivascular Cell") ~ "Perivascular",
        x %in% c("Hofbauer Cell") ~ "Hofbauer",
        x %in% c("NK Cell") ~ "NK Cell",
        x %in% c("Decidual NK Cell") ~ "Decidual NK",
        TRUE ~ x
    )
}

clean_pathway_label <- function(x) {
    x <- as.character(x)
    x <- str_replace_all(x, "_", " ")
    x <- toupper(x)
    
    dplyr::case_when(
        str_detect(x, "COLLAGEN") ~ "COLLAGEN",
        str_detect(x, "LAMININ") ~ "LAMININ",
        str_detect(x, "^FN1$|FIBRONECTIN") ~ "FN1",
        str_detect(x, "THBS") ~ "THBS",
        str_detect(x, "TENASCIN") ~ "TENASCIN",
        str_detect(x, "SPP1") ~ "SPP1",
        str_detect(x, "ICAM") ~ "ICAM",
        str_detect(x, "VCAM") ~ "VCAM",
        str_detect(x, "CCL") ~ "CCL",
        str_detect(x, "CXCL") ~ "CXCL",
        str_detect(x, "IL1") ~ "IL1",
        str_detect(x, "TNF") ~ "TNF",
        TRUE ~ x
    )
}

# ============================================================
# COLLAPSE MINOR PATHWAYS
# ============================================================

collapse_minor_pathways_A <- function(dfA,
                                      min_pathway_strength = 0.004,
                                      major_pathways = major_pathways,
                                      other_label = other_label) {
    dfA <- dfA %>%
        mutate(
            sender = clean_cell_label(sender),
            pathway = clean_pathway_label(pathway),
            strength = as.numeric(strength)
        ) %>%
        filter(!is.na(sender), !is.na(pathway), !is.na(strength), strength > 0)
    
    pathway_totals <- dfA %>%
        group_by(pathway) %>%
        summarise(total_strength = sum(strength, na.rm = TRUE), .groups = "drop")
    
    keep_pathways <- pathway_totals %>%
        filter(pathway %in% major_pathways | total_strength >= min_pathway_strength) %>%
        pull(pathway)
    
    dfA %>%
        mutate(pathway = if_else(pathway %in% keep_pathways, pathway, other_label)) %>%
        group_by(sender, pathway) %>%
        summarise(strength = sum(strength, na.rm = TRUE), .groups = "drop")
}

collapse_minor_pathways_B <- function(dfB,
                                      min_pathway_strength = 0.004,
                                      major_pathways = major_pathways,
                                      other_label = other_label) {
    dfB <- dfB %>%
        mutate(
            pathway = clean_pathway_label(pathway),
            receiver = clean_cell_label(receiver),
            strength = as.numeric(strength)
        ) %>%
        filter(!is.na(pathway), !is.na(receiver), !is.na(strength), strength > 0)
    
    pathway_totals <- dfB %>%
        group_by(pathway) %>%
        summarise(total_strength = sum(strength, na.rm = TRUE), .groups = "drop")
    
    keep_pathways <- pathway_totals %>%
        filter(pathway %in% major_pathways | total_strength >= min_pathway_strength) %>%
        pull(pathway)
    
    dfB %>%
        mutate(pathway = if_else(pathway %in% keep_pathways, pathway, other_label)) %>%
        group_by(pathway, receiver) %>%
        summarise(strength = sum(strength, na.rm = TRUE), .groups = "drop")
}

# ============================================================
# PREP TABLES FROM CELLCHAT OBJECT
# ============================================================

prep_tabs_for_plot <- function(cc, min_pathway_strength = 0.004) {
    tabs <- build_alluvial_tables(cc, pathways_keep)
    
    dfA <- collapse_minor_pathways_A(
        tabs$A,
        min_pathway_strength = min_pathway_strength
    )
    
    dfB <- collapse_minor_pathways_B(
        tabs$B,
        min_pathway_strength = min_pathway_strength
    )
    
    pathway_levels <- c(
        major_pathways[major_pathways %in% union(dfA$pathway, dfB$pathway)],
        if (other_label %in% union(dfA$pathway, dfB$pathway)) other_label
    )
    
    sender_order <- dfA %>%
        group_by(sender) %>%
        summarise(total = sum(strength, na.rm = TRUE), .groups = "drop") %>%
        arrange(total) %>%
        pull(sender)
    
    receiver_order <- dfB %>%
        group_by(receiver) %>%
        summarise(total = sum(strength, na.rm = TRUE), .groups = "drop") %>%
        arrange(total) %>%
        pull(receiver)
    
    dfA <- dfA %>%
        mutate(
            sender = factor(sender, levels = sender_order),
            pathway = factor(pathway, levels = pathway_levels)
        )
    
    dfB <- dfB %>%
        mutate(
            pathway = factor(pathway, levels = pathway_levels),
            receiver = factor(receiver, levels = receiver_order)
        )
    
    list(A = dfA, B = dfB)
}

# ============================================================
# COLORS
# ============================================================

sender_cols <- c(
    "STB" = "#5050FF",
    "Villous CTB" = "#CE3D32",
    "EVT" = "#749B58",
    "Fibroblast" = "#F0E685",
    "Perivascular" = "#466983",
    "Endothelial" = "#BA6338",
    "Hofbauer" = "#5DB1DD",
    "Macrophage" = "#802268",
    "NK Cell" = "#6BD76B",
    "Decidual NK" = "#D595A7",
    "ILC3" = "#924822",
    "Plasma Cell" = "#837B8D",
    "B Cell" = "#C75127",
    "T Cell" = "#D58F5C",
    "Treg" = "#7A65A5"
)

pathway_cols <- c(
    "COLLAGEN" = "#802268",
    "LAMININ" = "#5050FF",
    "FN1" = "#6BD76B",
    "THBS" = "#CE3D32",
    "TENASCIN" = "#749B58",
    "SPP1" = "#D58F5C",
    "ICAM" = "#F0E685",
    "VCAM" = "#BA6338",
    "Other core pathways" = "#BFBFBF"
)

# ============================================================
# PLOT FUNCTIONS
# ============================================================

plot_sender_to_pathway <- function(dfA, title_text) {
    ggplot(dfA, aes(axis1 = sender, axis2 = pathway, y = strength)) +
        geom_alluvium(
            aes(fill = sender),
            width = 0.22,
            alpha = 0.88,
            knot.pos = 0.35
        ) +
        geom_stratum(
            width = 0.28,
            fill = "grey95",
            color = "grey45",
            linewidth = 0.35
        ) +
        geom_text(
            stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 4.1,
            lineheight = 0.9
        ) +
        scale_x_discrete(
            limits = c("Sender cell", "Signaling pathway"),
            expand = c(0.06, 0.06)
        ) +
        scale_fill_manual(values = sender_cols, drop = FALSE) +
        labs(
            title = title_text,
            x = NULL,
            y = "Interaction strength"
        ) +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 15) +
        theme(
            text = element_text(family = "Arial"),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
            axis.text.x = element_text(face = "bold", color = "black", size = 14),
            axis.text.y = element_text(color = "black", size = 11),
            axis.title.y = element_text(color = "black", size = 12),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none",
            plot.margin = margin(12, 36, 12, 12)
        )
}

plot_pathway_to_receiver <- function(dfB, title_text) {
    ggplot(dfB, aes(axis1 = pathway, axis2 = receiver, y = strength)) +
        geom_alluvium(
            aes(fill = pathway),
            width = 0.22,
            alpha = 0.88,
            knot.pos = 0.35
        ) +
        geom_stratum(
            width = 0.28,
            fill = "grey95",
            color = "grey45",
            linewidth = 0.35
        ) +
        geom_text(
            stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 4.1,
            lineheight = 0.9
        ) +
        scale_x_discrete(
            limits = c("Signaling pathway", "Receiver cell"),
            expand = c(0.06, 0.06)
        ) +
        scale_fill_manual(values = pathway_cols, drop = FALSE) +
        labs(
            title = title_text,
            x = NULL,
            y = "Interaction strength"
        ) +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 15) +
        theme(
            text = element_text(family = "Arial"),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
            axis.text.x = element_text(face = "bold", color = "black", size = 14),
            axis.text.y = element_text(color = "black", size = 11),
            axis.title.y = element_text(color = "black", size = 12),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none",
            plot.margin = margin(12, 36, 12, 12)
        )
}

# ============================================================
# MAKE AND SAVE PANELS
# ============================================================

make_condition_alluvials <- function(cc,
                                     condition_label,
                                     title_prefix,
                                     min_pathway_strength = 0.004,
                                     width = 13,
                                     height = 8.5) {
    tabs_plot <- prep_tabs_for_plot(
        cc,
        min_pathway_strength = min_pathway_strength
    )
    
    pA <- plot_sender_to_pathway(
        tabs_plot$A,
        paste0(title_prefix, ": sender contributions to signaling pathways")
    )
    
    pB <- plot_pathway_to_receiver(
        tabs_plot$B,
        paste0(title_prefix, ": signaling pathways to receiver cells")
    )
    
    write_csv(
        tabs_plot$A,
        file.path(out_ab_pub, paste0("final_alluvial_A_sender_by_pathway_", condition_label, ".csv"))
    )
    
    write_csv(
        tabs_plot$B,
        file.path(out_ab_pub, paste0("final_alluvial_B_pathway_by_receiver_", condition_label, ".csv"))
    )
    
    ggsave(
        file.path(out_ab_pub, paste0("SuppFig5_", condition_label, "_A_sender_to_pathway.png")),
        pA,
        width = width,
        height = height,
        dpi = 600,
        bg = "white",
        limitsize = FALSE
    )
    
    ggsave(
        file.path(out_ab_pub, paste0("SuppFig5_", condition_label, "_A_sender_to_pathway.pdf")),
        pA,
        width = width,
        height = height,
        bg = "white",
        limitsize = FALSE
    )
    
    ggsave(
        file.path(out_ab_pub, paste0("SuppFig5_", condition_label, "_B_pathway_to_receiver.png")),
        pB,
        width = width,
        height = height,
        dpi = 600,
        bg = "white",
        limitsize = FALSE
    )
    
    ggsave(
        file.path(out_ab_pub, paste0("SuppFig5_", condition_label, "_B_pathway_to_receiver.pdf")),
        pB,
        width = width,
        height = height,
        bg = "white",
        limitsize = FALSE
    )
    
    list(A = pA, B = pB, tabs = tabs_plot)
}

# ============================================================
# PATCH: avoid recursive default argument references
# ============================================================

MAJOR_PATHWAYS_SUPP5 <- c(
  "COLLAGEN", "LAMININ", "FN1", "THBS",
  "TENASCIN", "SPP1", "ICAM", "VCAM"
)

OTHER_PATHWAY_LABEL_SUPP5 <- "Other core pathways"

collapse_minor_pathways_A <- function(dfA,
                                      min_pathway_strength = 0.004,
                                      major_pathways_use = MAJOR_PATHWAYS_SUPP5,
                                      other_label_use = OTHER_PATHWAY_LABEL_SUPP5) {
  dfA <- dfA %>%
    dplyr::mutate(
      sender = clean_cell_label(sender),
      pathway = clean_pathway_label(pathway),
      strength = as.numeric(strength)
    ) %>%
    dplyr::filter(
      !is.na(sender),
      !is.na(pathway),
      !is.na(strength),
      strength > 0
    )

  pathway_totals <- dfA %>%
    dplyr::group_by(pathway) %>%
    dplyr::summarise(
      total_strength = sum(strength, na.rm = TRUE),
      .groups = "drop"
    )

  keep_pathways <- pathway_totals %>%
    dplyr::filter(
      pathway %in% major_pathways_use |
        total_strength >= min_pathway_strength
    ) %>%
    dplyr::pull(pathway)

  dfA %>%
    dplyr::mutate(
      pathway = dplyr::if_else(
        pathway %in% keep_pathways,
        pathway,
        other_label_use
      )
    ) %>%
    dplyr::group_by(sender, pathway) %>%
    dplyr::summarise(
      strength = sum(strength, na.rm = TRUE),
      .groups = "drop"
    )
}

collapse_minor_pathways_B <- function(dfB,
                                      min_pathway_strength = 0.004,
                                      major_pathways_use = MAJOR_PATHWAYS_SUPP5,
                                      other_label_use = OTHER_PATHWAY_LABEL_SUPP5) {
  dfB <- dfB %>%
    dplyr::mutate(
      pathway = clean_pathway_label(pathway),
      receiver = clean_cell_label(receiver),
      strength = as.numeric(strength)
    ) %>%
    dplyr::filter(
      !is.na(pathway),
      !is.na(receiver),
      !is.na(strength),
      strength > 0
    )

  pathway_totals <- dfB %>%
    dplyr::group_by(pathway) %>%
    dplyr::summarise(
      total_strength = sum(strength, na.rm = TRUE),
      .groups = "drop"
    )

  keep_pathways <- pathway_totals %>%
    dplyr::filter(
      pathway %in% major_pathways_use |
        total_strength >= min_pathway_strength
    ) %>%
    dplyr::pull(pathway)

  dfB %>%
    dplyr::mutate(
      pathway = dplyr::if_else(
        pathway %in% keep_pathways,
        pathway,
        other_label_use
      )
    ) %>%
    dplyr::group_by(pathway, receiver) %>%
    dplyr::summarise(
      strength = sum(strength, na.rm = TRUE),
      .groups = "drop"
    )
}

prep_tabs_for_plot <- function(cc, min_pathway_strength = 0.004) {
  tabs <- build_alluvial_tables(cc, pathways_keep)

  dfA <- collapse_minor_pathways_A(
    tabs$A,
    min_pathway_strength = min_pathway_strength,
    major_pathways_use = MAJOR_PATHWAYS_SUPP5,
    other_label_use = OTHER_PATHWAY_LABEL_SUPP5
  )

  dfB <- collapse_minor_pathways_B(
    tabs$B,
    min_pathway_strength = min_pathway_strength,
    major_pathways_use = MAJOR_PATHWAYS_SUPP5,
    other_label_use = OTHER_PATHWAY_LABEL_SUPP5
  )

  pathway_levels <- c(
    MAJOR_PATHWAYS_SUPP5[MAJOR_PATHWAYS_SUPP5 %in% union(dfA$pathway, dfB$pathway)],
    if (OTHER_PATHWAY_LABEL_SUPP5 %in% union(dfA$pathway, dfB$pathway)) {
      OTHER_PATHWAY_LABEL_SUPP5
    }
  )

  sender_order <- dfA %>%
    dplyr::group_by(sender) %>%
    dplyr::summarise(total = sum(strength, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(total) %>%
    dplyr::pull(sender)

  receiver_order <- dfB %>%
    dplyr::group_by(receiver) %>%
    dplyr::summarise(total = sum(strength, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(total) %>%
    dplyr::pull(receiver)

  dfA <- dfA %>%
    dplyr::mutate(
      sender = factor(sender, levels = sender_order),
      pathway = factor(pathway, levels = pathway_levels)
    )

  dfB <- dfB %>%
    dplyr::mutate(
      pathway = factor(pathway, levels = pathway_levels),
      receiver = factor(receiver, levels = receiver_order)
    )

  list(A = dfA, B = dfB)
}

plot_sender_to_pathway <- function(dfA, title_text) {
  ggplot(dfA, aes(axis1 = sender, axis2 = pathway, y = strength)) +
    geom_alluvium(
      aes(fill = sender),
      width = 0.22,
      alpha = 0.88,
      knot.pos = 0.35
    ) +
    geom_stratum(
      width = 0.28,
      fill = "grey95",
      color = "grey45",
      linewidth = 0.35
    ) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 4.1,
      lineheight = 0.9
    ) +
    scale_x_discrete(
      limits = c("Sender cell", "Signaling pathway"),
      expand = c(0.06, 0.06)
    ) +
    scale_fill_manual(values = sender_cols, drop = FALSE) +
    labs(
      title = title_text,
      x = NULL,
      y = "Interaction strength"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 15, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      axis.text.x = element_text(face = "bold", color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.y = element_text(color = "black", size = 12),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(12, 36, 12, 12)
    )
}

plot_pathway_to_receiver <- function(dfB, title_text) {
  ggplot(dfB, aes(axis1 = pathway, axis2 = receiver, y = strength)) +
    geom_alluvium(
      aes(fill = pathway),
      width = 0.22,
      alpha = 0.88,
      knot.pos = 0.35
    ) +
    geom_stratum(
      width = 0.28,
      fill = "grey95",
      color = "grey45",
      linewidth = 0.35
    ) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 4.1,
      lineheight = 0.9
    ) +
    scale_x_discrete(
      limits = c("Signaling pathway", "Receiver cell"),
      expand = c(0.06, 0.06)
    ) +
    scale_fill_manual(values = pathway_cols, drop = FALSE) +
    labs(
      title = title_text,
      x = NULL,
      y = "Interaction strength"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 15, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      axis.text.x = element_text(face = "bold", color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.y = element_text(color = "black", size = 12),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(12, 36, 12, 12)
    )
}

# Run final panels
preterm_high_final <- make_condition_alluvials(
    cc_preterm_high,
    condition_label = "Preterm_High",
    title_prefix = "Preterm high-MNP",
    min_pathway_strength = 0.009
)

term_high_final <- make_condition_alluvials(
    cc_term_high,
    condition_label = "Term_High",
    title_prefix = "Term high-MNP",
    min_pathway_strength = 0.009
)

# Preview
preterm_high_final$A
preterm_high_final$B
term_high_final$A
term_high_final$B







