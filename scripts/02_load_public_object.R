public_object_path <- config$public_object_path
if (!file.exists(public_object_path)) {
  stop("Public object not found: ", public_object_path,
       "\nPlace spatial_mnp_publication_processed_seurat_v1.rds in data/external/ or update config/analysis_config.yml.")
}

obj <- readRDS(public_object_path)
validation <- validate_public_object(obj, config, stop_on_missing = FALSE)
readr::write_csv(validation, file.path(config$log_dir, "public_object_metadata_validation.csv"))

if (any(!validation$present)) {
  warning("Some configured metadata columns are missing. See outputs/logs/public_object_metadata_validation.csv")
}
