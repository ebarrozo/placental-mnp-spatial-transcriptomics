subject_level_mean <- function(meta, value_col, group_cols) {
  meta |>
    dplyr::filter(!is.na(.data[[value_col]])) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
}

within_delivery_contrast <- function(df, value_col = "value", delivery_col = "Delivery_class",
                                     group_col = "MNPBurden", high_label = "High", low_label = "Low") {
  df |>
    dplyr::group_by(.data[[delivery_col]]) |>
    dplyr::summarise(
      n_high = sum(.data[[group_col]] == high_label, na.rm = TRUE),
      n_low = sum(.data[[group_col]] == low_label, na.rm = TRUE),
      mean_high = mean(.data[[value_col]][.data[[group_col]] == high_label], na.rm = TRUE),
      mean_low = mean(.data[[value_col]][.data[[group_col]] == low_label], na.rm = TRUE),
      estimate_high_minus_low = mean_high - mean_low,
      p_value = tryCatch(stats::wilcox.test(.data[[value_col]] ~ .data[[group_col]])$p.value, error = function(e) NA_real_),
      .groups = "drop"
    ) |>
    dplyr::mutate(q_value = p.adjust(p_value, method = "BH"))
}
