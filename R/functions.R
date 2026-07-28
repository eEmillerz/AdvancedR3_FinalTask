preprocess <- function(data) {
  data |>
    dplyr::select(
      -c(grade, other_liver_disease, weight, age, hu, fat_perc, fibrosis_perc, uric_acid)
    ) |>
    tidyr::pivot_longer(
      cols = -c(id, lipidosis)) |>
    dplyr::rename("metabolite" = name)
}
preprocess(lipoproteins)
