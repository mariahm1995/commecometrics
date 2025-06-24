#' Run an ecometric model for qualitative environmental variables
#'
#' Builds an ecometric trait space for qualitative environmental variables,
#' estimating the most probable category and the probability of each category
#' at each trait bin combination. Also calculates prediction accuracy
#' and anomalies for each point.
#'
#' @param points_df Output first element of the list from \code{summarize_traits_by_point()}. A data frame with columns: `summ_trait_1`, `summ_trait_2`, `count_trait`, and the environmental variable.
#' @param category_col Name of the column containing the categorical trait.
#' @param grid_bins_1 Number of bins for the first trait axis. If `NULL` (default),
#'   the number is calculated automatically using Scott's rule via `optimal_bins()`.
#' @param grid_bins_2 Number of bins for the second trait axis. If `NULL` (default),
#'   the number is calculated automatically using Scott's rule via `optimal_bins()`.
#' @param min_species Minimum number of species with trait data per point (default = 3).
#' @return A list containing:
#' \item{points_df}{Filtered input data frame with the following added columns:
#'   \describe{
#'     \item{bin_1}{Bin assignment code for first trait axis.}
#'     \item{bin_2}{Bin assignment code for second trait axis.}
#'     \item{prob_<category>}{Estimated probability of each environmental category per trait bin (e.g., \code{prob_1}, \code{prob_2}, etc.).}
#'     \item{observed_probability}{Probability assigned to the observed category for each point.}
#'     \item{predicted_probability}{Probability assigned to the predicted (most likely) category for each point.}
#'     \item{predicted_category}{Predicted environmental category for each point.}
#'     \item{correct_prediction}{Indicator for whether the predicted category matches the observed category (\code{"Yes"} or \code{"No"}).}
#'     \item{env_anom}{Difference between predicted and observed category probabilities.}
#'   }
#' }
#' \item{eco_space}{Raster-format data frame representing trait space bins with estimated environmental categories.}
#' \item{diagnostics}{Summary stats about bin usage and data coverage.}
#' \item{settings}{Metadata including the modeled trait.}
#' \item{prediction_accuracy}{Overall percentage of correct predictions.}
#' @importFrom stats density setNames
#' @examples
#' \dontrun{
#' # Load internal data
#' data("geoPoints", package = "commecometrics")
#' data("traits", package = "commecometrics")
#' data("spRanges", package = "commecometrics")
#'
#' # Step 1: Summarize trait values at sampling points
#' traitsByPoint <- summarize_traits_by_point(
#'   points_df = geoPoints,
#'   trait_df = traits,
#'   species_polygons = spRanges,
#'   trait_column = "RBL",
#'   species_name_col = "TaxonName",
#'   continent = FALSE,
#'   parallel = FALSE
#' )
#'
#' # Step 2: Run ecometric model using land cover class as qualitative variable
#' modelResult <- ecometric_model_qual(
#'   points_df = traitsByPoint$points,
#'   category_col = "vegetation",
#'   min_species = 3
#' )
#'
#' # View the percentage of correctly predicted categories
#' print(modelResult$prediction_accuracy)
#' }
#' @export
ecometric_model_qual <- function(points_df,
                                 category_col,
                                 grid_bins_1 = NULL,
                                 grid_bins_2 = NULL,
                                 min_species = 3) {
  # Remove NAs in trait category
  points_df <- points_df %>% dplyr::filter(!is.na(.data[[category_col]]))

  # Coerce category column to character if needed
  if (!is.character(points_df[[category_col]])) {
    message("Converting '", category_col, "' to character.")
    points_df[[category_col]] <- as.character(points_df[[category_col]])
  }

  if (!all(c("summ_trait_1", "summ_trait_2", "count_trait") %in% names(points_df))) {
    stop("points_df must contain 'summ_trait_1', 'summ_trait_2', and 'count_trait'.")
  }

  # Filter low-coverage points
  message("Filtering points with at least ", min_species, " species...")
  filtered_df <- dplyr::filter(points_df, count_trait >= min_species)
  prop_retained <- round(nrow(filtered_df) / nrow(points_df) * 100, 1)
  message("Retained ", nrow(filtered_df), " of ", nrow(points_df), " points (", prop_retained, "%)")
  if (nrow(filtered_df) == 0) stop("No points retained. Try lowering 'min_species'.")

  # Determine bin numbers if not provided
  if (is.null(grid_bins_1)) {
    grid_bins_1 <- optimal_bins(filtered_df$summ_trait_1)
    message("Optimal number of bins for the first trait summary metric (Scott's rule): ", grid_bins_1)
  }
  if (is.null(grid_bins_2)) {
    grid_bins_2 <- optimal_bins(filtered_df$summ_trait_2)
    message("Optimal number of bins for the second trait summary metric (Scott's rule): ", grid_bins_2)
  }

  # Define bin breaks
  mrange <- range(filtered_df$summ_trait_1, na.rm = TRUE)
  sdrange <- range(filtered_df$summ_trait_2, na.rm = TRUE)
  mbrks <- seq(mrange[1] - 0.001, mrange[2] + 0.001, length.out = grid_bins_1 + 1)
  sdbrks <- seq(sdrange[1] - 0.001, sdrange[2] + 0.001, length.out = grid_bins_2 + 1)

  filtered_df <- filtered_df %>%
    dplyr::mutate(
      bin_1 = .bincode(summ_trait_1, mbrks),
      bin_2 = .bincode(summ_trait_2, sdbrks)
    )

  categories <- sort(unique(filtered_df[[category_col]]))

  # Build bin grid
  eco_list <- list()
  bin_counts <- matrix(0, nrow = grid_bins_2, ncol = grid_bins_1)

  for (i in 1:grid_bins_1) {
    for (j in 1:grid_bins_2) {
      idx <- which(filtered_df$bin_1 == i & filtered_df$bin_2 == j)
      dat <- filtered_df[[category_col]][idx]
      bin_counts[j, i] <- length(dat)

      if (length(dat) > 0) {
        tab <- table(factor(dat, levels = categories))
        probs <- as.numeric(tab) / sum(tab)
        mode_cat <- categories[which.max(probs)]

        eco_list[[length(eco_list) + 1]] <- tibble::tibble(
          bin_1 = i,
          bin_2 = j,
          env_est = mode_cat,
          !!!setNames(as.list(probs), paste0("prob_", categories))
        )
      }
    }
  }

  eco_space <- dplyr::bind_rows(eco_list)

  # Bin diagnostics
  bin_counts_flipped <- bin_counts[grid_bins_2:1, , drop = FALSE]
  used_bins <- sum(bin_counts_flipped > 0)
  total_bins <- grid_bins_1 * grid_bins_2
  message("Used ", used_bins, " of ", total_bins, " bins (", round(used_bins / total_bins * 100, 1), "%)")

  # Map predictions back to points
  filtered_df <- dplyr::left_join(filtered_df, eco_space, by = c("bin_1", "bin_2"))

  # Compute probabilities for observed/predicted
  observed_prob_col <- paste0("prob_", filtered_df[[category_col]])
  predicted_prob_col <- paste0("prob_", filtered_df$env_est)

  observed_prob_val <- purrr::map2_dbl(observed_prob_col, seq_len(nrow(filtered_df)), function(colname, idx) {
    if (colname %in% names(filtered_df)) filtered_df[[colname]][idx] else NA_real_
  })

  predicted_prob_val <- purrr::map2_dbl(predicted_prob_col, seq_len(nrow(filtered_df)), function(colname, idx) {
    if (colname %in% names(filtered_df)) filtered_df[[colname]][idx] else NA_real_
  })

  filtered_df$observed_probability <- observed_prob_val
  filtered_df$predicted_probability <- predicted_prob_val

  filtered_df <- filtered_df %>%
    dplyr::mutate(
      predicted_category = env_est,
      correct_prediction = ifelse(.data[[category_col]] == env_est, "Yes", "No"),
      env_anom = dplyr::case_when(
        is.na(predicted_category) ~ NA_real_,
        TRUE ~ predicted_probability - observed_probability
      )
    )

  correct_predictions <- sum(filtered_df$correct_prediction == "Yes", na.rm = TRUE)
  total_points <- nrow(filtered_df)
  percent_correct <- round((correct_predictions / total_points) * 100, 2)
  message("Prediction accuracy: ", percent_correct, "% correct predictions.")

  return(list(
    points_df = filtered_df,
    eco_space = eco_space,
    diagnostics = list(
      total_points = nrow(points_df),
      retained_points = nrow(filtered_df),
      proportion_retained = prop_retained,
      used_bins = used_bins,
      total_bins = total_bins,
      bin_counts = bin_counts,
      brks_1 = mbrks,
      brks_2 = sdbrks,
      categories = categories
    ),
    settings = list(
      trait_var = category_col
    ),
    prediction_accuracy = percent_correct
  ))
}
