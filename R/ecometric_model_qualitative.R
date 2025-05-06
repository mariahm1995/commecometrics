#' Run an ecometric model for qualitative environmental variables
#'
#' Builds an ecometric trait space for qualitative environmental variables,
#' estimating the most probable category and the probability of each category
#' at each trait bin combination. Also calculates prediction accuracy
#' and anomalies based on observed category values.
#'
#' @param points_df A data frame with columns: `mean_trait`, `sd_trait`, `count_trait`, and the categorical trait.
#' @param category_col Name of the column containing the categorical trait.
#' @param grid_bins_mean Number of bins for the trait mean axis (default = Scott's rule).
#' @param grid_bins_sd Number of bins for the trait SD axis (default = Scott's rule).
#' @param min_species Minimum number of species with trait data per point (default = 3).
#'
#' @return A list containing:
#' \item{points_df}{Filtered input data frame with predicted and observed environmental values, bin assignments, and anomalies.}
#' \item{eco_space}{Raster-format data frame representing trait space bins with estimated environmental categories.}
#' \item{diagnostics}{Summary stats about bin usage and data coverage.}
#' \item{settings}{Metadata including the modeled trait.}
#' \item{prediction_accuracy}{Overall percentage of correct predictions.}
#' @importFrom stats density setNames
#' @export
ecometric_model_qual <- function(points_df,
                                        category_col,
                                        grid_bins_mean = NULL,
                                        grid_bins_sd = NULL,
                                        min_species = 3) {

  # Remove NAs in trait category
  points_df <- points_df %>% dplyr::filter(!is.na(.data[[category_col]]))

  if (!all(c("mean_trait", "sd_trait", "count_trait") %in% names(points_df))) {
    stop("points_df must contain 'mean_trait', 'sd_trait', and 'count_trait'.")
  }

  # Filter low-coverage points
  message("Filtering points with at least ", min_species, " species...")
  filtered_df <- dplyr::filter(points_df, count_trait >= min_species)
  prop_retained <- round(nrow(filtered_df) / nrow(points_df) * 100, 1)
  message("Retained ", nrow(filtered_df), " of ", nrow(points_df), " points (", prop_retained, "%)")
  if (nrow(filtered_df) == 0) stop("No points retained. Try lowering 'min_species'.")

  # Determine bin numbers if not provided
  if (is.null(grid_bins_mean)) {
    grid_bins_mean <- optimal_bins_scott(filtered_df$mean_trait)
    message("Optimal bins for mean trait (Scott): ", grid_bins_mean)
  }
  if (is.null(grid_bins_sd)) {
    grid_bins_sd <- optimal_bins_scott(filtered_df$sd_trait)
    message("Optimal bins for SD trait (Scott): ", grid_bins_sd)
  }

  # Define bin breaks
  mrange <- range(filtered_df$mean_trait, na.rm = TRUE)
  sdrange <- range(filtered_df$sd_trait, na.rm = TRUE)
  mbrks <- seq(mrange[1] - 0.001, mrange[2] + 0.001, length.out = grid_bins_mean + 1)
  sdbrks <- seq(sdrange[1] - 0.001, sdrange[2] + 0.001, length.out = grid_bins_sd + 1)

  filtered_df <- filtered_df %>%
    dplyr::mutate(
      mbc = .bincode(mean_trait, mbrks),
      sdc = .bincode(sd_trait, sdbrks)
    )

  categories <- sort(unique(filtered_df[[category_col]]))

  # Build bin grid
  eco_list <- list()
  bin_counts <- matrix(0, nrow = grid_bins_sd, ncol = grid_bins_mean)

  for (i in 1:grid_bins_mean) {
    for (j in 1:grid_bins_sd) {
      idx <- which(filtered_df$mbc == i & filtered_df$sdc == j)
      dat <- filtered_df[[category_col]][idx]
      bin_counts[j, i] <- length(dat)

      if (length(dat) > 0) {
        tab <- table(factor(dat, levels = categories))
        probs <- as.numeric(tab) / sum(tab)
        mode_cat <- categories[which.max(probs)]

        eco_list[[length(eco_list) + 1]] <- tibble::tibble(
          mbc = i,
          sdc = j,
          env_est = mode_cat,
          !!!setNames(as.list(probs), paste0("prob_", categories))
        )
      }
    }
  }

  eco_space <- dplyr::bind_rows(eco_list)

  # Bin diagnostics
  bin_counts_flipped <- bin_counts[grid_bins_sd:1, , drop = FALSE]
  used_bins <- sum(bin_counts_flipped > 0)
  total_bins <- grid_bins_mean * grid_bins_sd
  message("Used ", used_bins, " of ", total_bins, " bins (", round(used_bins / total_bins * 100, 1), "%)")

  # Map predictions back to points
  filtered_df <- dplyr::left_join(filtered_df, eco_space, by = c("mbc", "sdc"))

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
      anomaly = dplyr::case_when(
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
      mbrks = mbrks,
      sdbrks = sdbrks,
      categories = categories
    ),
    settings = list(
      trait_var = category_col
    ),
    prediction_accuracy = percent_correct
  ))
}
