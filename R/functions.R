#' Preprocess data to make it suitable for creating a table of descriptive statistics
#'
#' @param data Lipoproteins data.
#'
#' @returns Preprocessed lipoproteins data (pivoted longer).
preprocess <- function(data) {
  data |>
    dplyr::select(
      -c(grade, other_liver_disease, weight, age, hu, fat_perc, fibrosis_perc, uric_acid)
    ) |>
    tidyr::pivot_longer(
      cols = -c(id, lipidosis)
    ) |>
    dplyr::rename("metabolite" = name)
}

#' Create a table of descriptive statistics for each metabolite
#'
#' @param data Lipoproteins data
#'
#' @returns A table of descriptive statistics for each metabolite.
create_table_descriptive_stats <- function(data) {
  data |>
    preprocess() |>
    dplyr::summarise(
      dplyr::across(
        -c(id, lipidosis),
        list(
          mean = mean, sd = sd, median = median,
          iqr = \(x) IQR(x, na.rm = TRUE)
        )
      ),
      .by = c(metabolite)
    ) |>
    dplyr::filter(!is.na(c(value_sd))) |>
    dplyr::mutate(across(where(is.numeric), \(x) round(x, digits = 1))) |>
    dplyr::mutate(MeanSD = glue::glue("{value_mean} ({value_sd})")) |>
    dplyr::mutate(MedianIQR = glue::glue("{value_median} ({value_iqr})")) |>
    dplyr::select(Metabolite = metabolite, "Mean (SD)" = MeanSD, "Median (IQR)" = MedianIQR)
}

#' Create a plot showing metabolite distributions
#'
#' @param data Lipoproteins data set.
#'
#' @returns Plot showing metabolite distributions.
create_plot_distributions <- function(data) {
  data |>
    preprocess() |>
    dplyr::filter(metabolite %in% c("insulin", "hdl", "ldl", "chol", "na", "k", "cl", "ca", "phosphate")) |>
    ggplot2::ggplot(aes(x = value)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(vars(metabolite), scales = "free") +
    ggplot2::theme(
      plot.title = element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      x = "Metabolite value",
      y = "Count",
      title = "Metabolite value distributions"
    )
}
