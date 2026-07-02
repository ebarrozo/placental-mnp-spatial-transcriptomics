# Master downstream reproducibility script.
# Assumes public processed Seurat object is available at:
#   data/external/spatial_mnp_publication_processed_seurat_v1.rds

message("Loading public object...")
source("scripts/02_load_public_object.R")

message("Preparing model-ready metadata index...")
source("scripts/03_prepare_analysis_metadata.R")

message("Running cell-type composition workflow...")
source("scripts/04_celltype_composition.R")

message("Running module-score workflow...")
source("scripts/05_module_scores.R")

message("Running exploratory polymer-module association workflow...")
source("scripts/07_polymer_associations.R")

message("Running macrophage-state workflow...")
source("scripts/06_macrophage_states.R")

message("Running communication source-table workflow...")
source("scripts/07_communication_figures.R")
source("scripts/08_xenium_explorer_exports.R")

message("Running annotation QC/support workflow...")
source("scripts/10_annotation_qc_support.R")

if (tolower(Sys.getenv("RUN_PSEUDOBULK", unset = "false")) %in% c("true", "1", "yes")) {
  message("Running pseudobulk DE/pathway-enrichment workflow...")
  source("scripts/09_pseudobulk_de_enrichment.R")
} else {
  message("Skipping pseudobulk DE/pathway-enrichment by default. Set RUN_PSEUDOBULK=true to run scripts/09_pseudobulk_de_enrichment.R")
}

message("Building Supplementary Table 1 workbook/source CSV index...")
source("scripts/11_build_supplementary_table1.R")

message("Main-figure assembly and final accepted-manuscript layout remain to be finalized.")
