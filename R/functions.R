preprocess <- function(data) {
  data |>
    dplyr::select(
      -c(grade, other_liver_disease, weight, age, hu, fat_perc, fibrosis_perc, uric_acid)
    ) |>
    tidyr::pivot_longer(
      cols = -c(id, lipidosis)) |>
    dplyr::rename("metabolite" = name)
}

create_table_descriptive_stats <- function(data) {
  data |>
    preprocess() |>
    dplyr::summarise(
      dplyr::across(
        -c(id, lipidosis),
        list(mean = mean, sd = sd, median = median,
             iqr = \(x) IQR(x, na.rm = TRUE))
      ),
      .by = c(metabolite)
    ) |>
    dplyr::filter(!is.na(c(value_sd))) |>
    dplyr::mutate(across(where(is.numeric), \(x) round(x, digits = 1))) |>
    dplyr::mutate(MeanSD = glue::glue("{value_mean} ({value_sd})")) |>
    dplyr::mutate(MedianIQR = glue::glue("{value_median} ({value_iqr})")) |>
    dplyr::select(Metabolite = metabolite, "Mean (SD)" = MeanSD, "Median (IQR)" = MedianIQR)
}
