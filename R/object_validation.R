required_metadata_columns <- function(config) {
  unique(c(
    config$sample_id_col,
    config$donor_id_col,
    config$primary_group_col,
    config$primary_delivery_col,
    config$celltype_col,
    config$platform_col,
    "cumulative_mnp_ug_per_g",
    "cumulative_mnp_log10ugpg",
    config$module_score_cols,
    config$polymer_ug_g_cols
  ))
}

validate_public_object <- function(obj, config, stop_on_missing = TRUE) {
  stopifnot(inherits(obj, "Seurat"))
  present <- colnames(obj@meta.data)
  req <- required_metadata_columns(config)
  out <- data.frame(
    column = req,
    present = req %in% present,
    stringsAsFactors = FALSE
  )
  missing <- out$column[!out$present]
  if (length(missing) && isTRUE(stop_on_missing)) {
    stop("Public object is missing required metadata columns: ", paste(missing, collapse = ", "))
  }
  out
}
