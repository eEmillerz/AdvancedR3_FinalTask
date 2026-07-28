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
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(ggplot2::vars(metabolite), scales = "free") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      x = "Metabolite value",
      y = "Count",
      title = "Metabolite value distributions"
    )
}

# Turning it into generalised function

#' Making data ready for running model
#'
#' @param data Data set with lipoproteins.
#'
#' @returns Ready data set for model analysis.

ready_data <- function(data) {
  data |>
    preprocess() |>
    dplyr::mutate(
      lipidosis = as.factor(lipidosis),
      value = as.vector(scale(value))
    )
}

#' Function to test model on a single metabolite
#'
#' @param data Lipoproteins data set.
#'
#' @returns Data frame with coefficients, sd and p val.
fit_model <- function(data) {
  glm(
    formula = lipidosis ~ value,
    data = data,
    family = binomial
  ) |>
    broom::tidy(exponentiate = TRUE) |>
    dplyr::mutate(
      metabolite = unique(data$metabolite),
      model = format(lipidosis ~ value),
      .before = tidyselect::everything()
    )
}

#' Create model results
#'
#' @param data Lipidomics data
#'
#' @returns Data frame with model results for insulin metabolite
create_model_results <- function(data) {
  data |>
    ready_data() |>
    dplyr::group_split(metabolite) |>
    purrr::map(fit_model) |>
    purrr::list_rbind()
}

#' Plot the estimates and standard errors of the model results
#'
#' @param results The model results.
#'
#' @returns A ggplot2 figure.
create_plot_model_results <- function(results) {
  results |>
    dplyr::filter(term == "value", std.error <= 2, estimate <= 5) |>
    dplyr::select(metabolite, model, estimate, std.error) |>
    ggplot2::ggplot(ggplot2::aes(
      x = estimate,
      y = metabolite,
      xmin = estimate - std.error,
      xmax = estimate + std.error
    )) +
    ggplot2::geom_pointrange() +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
    ggplot2::facet_grid(cols = ggplot2::vars(model))
}
