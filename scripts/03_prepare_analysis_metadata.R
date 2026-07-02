meta <- obj@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::mutate(
    Delivery_class = factor(.data[[config$primary_delivery_col]], levels = c("Term", "Preterm")),
    MNPBurden = factor(.data[[config$primary_group_col]], levels = c("Low", "High")),
    celltype_display = .data[[config$celltype_col]],
    platform_clean = .data[[config$platform_col]],
    donor_id = .data[[config$donor_id_col]],
    sample_id_standardized = .data[[config$sample_id_col]]
  )

readr::write_csv(
  meta |>
    dplyr::distinct(donor_id, sample_id_standardized, Delivery_class, MNPBurden, platform_clean, .keep_all = FALSE),
  file.path(config$table_dir, "model_ready_sample_run_index.csv")
)
