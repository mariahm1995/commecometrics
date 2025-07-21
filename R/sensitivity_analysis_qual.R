#' Perform sensitivity analysis on ecometric models (qualitative environmental variables)
#'
#' Evaluates how varying sample sizes affect the performance of ecometric models,
#' focusing on two aspects:
#' \itemize{
#'   \item \strong{Sensitivity (internal consistency)}: How accurately the model predicts environmental conditions
#'         on the same data it was trained on.
#'   \item \strong{Transferability (external applicability)}: How well the model performs on unseen data.
#' }
#' It tests different sample sizes by resampling the data multiple times (bootstrap iterations),
#' training an ecometric model on each subset, and evaluating prediction error and correlation.
#'
#' Two plots are generated:
#' \enumerate{
#'   \item \strong{Training Accuracy vs. Sample size:} Reflects internal model consistency.
#'   \item \strong{Testing Accuracy vs. Sample size:} Reflects external model performance.
#' }
#'
#' Parallel processing is supported to speed up the analysis.
#'
#' @param points_df Output first element of the list from \code{summarize_traits_by_point()}. A data frame with columns: `summ_trait_1`, `summ_trait_2`, `count_trait`, and the environmental variable.
#' @param category_col Name of the column containing the categorical trait.
#' @param sample_sizes Numeric vector specifying the number of communities (sampling points)
#'   to evaluate in the sensitivity analysis. For each value, a random subset of the data of that
#'   size is drawn without replacement and then split into training and testing sets using the
#'   proportion defined by `test_split` (default is 80% training, 20% testing).
#'   All values in `sample_sizes` must be less than or equal to the number of rows in `points_df`,
#'   and large enough to allow splitting based on `test_split` (i.e., both the training and testing
#'   sets must contain at 30 communities).
#' @param iterations Number of bootstrap iterations per sample size (default = 20).
#' @param test_split Proportion of data to use for testing (default = 0.2).
#' @param grid_bins_1  Number of bins for the first trait axis. If `NULL` (default),
#'   the number is calculated automatically using Scott's rule via `optimal_bins()`.
#' @param grid_bins_2  Number of bins for the second trait axis. If `NULL` (default),
#'   the number is calculated automatically using Scott's rule via `optimal_bins()`.
#' @param parallel Logical; whether to run iterations in parallel (default = TRUE).
#' @param n_cores Number of cores for parallelization (default = detectCores() - 1).
#'
#' @return A list containing:
#'   \item{combined_results}{All raw iteration results.}
#'   \item{summary_results}{Mean accuracy per sample size.}
#'
#' @examples
#' \donttest{
#' # Load internal data
#' data("geoPoints", package = "commecometrics")
#' data("traits", package = "commecometrics")
#' data("spRanges", package = "commecometrics")
#'
#' # Summarize trait values at sampling points
#' traitsByPoint <- summarize_traits_by_point(
#'   points_df = geoPoints,
#'   trait_df = traits,
#'   species_polygons = spRanges,
#'   trait_column = "RBL",
#'   species_name_col = "sci_name",
#'   continent = FALSE,
#'   parallel = FALSE
#' )
#'
#' # Run sensitivity analysis for dominant land cover class
#' sensitivityQual <- sensitivity_analysis_qual(
#'   points_df = traitsByPoint$points,
#'   category_col = "vegetation",
#'   sample_sizes = seq(40, 90, 10),
#'   iterations = 5,
#'   parallel = FALSE
#' )
#'
#' # View results
#' head(sensitivityQual$summary_results)
#' }
#' @export
sensitivity_analysis_qual <- function(points_df,
                                      category_col,
                                      sample_sizes,
                                      iterations = 20,
                                      test_split = 0.2,
                                      grid_bins_1 = NULL,
                                      grid_bins_2 = NULL,
                                      parallel = TRUE,
                                      n_cores = parallel::detectCores() - 1) {
  # Helper: mode function
  get_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  # Load required package
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required but not installed.")
  }

  # Check for required columns
  required_cols <- c("summ_trait_1", "summ_trait_2", "count_trait", category_col)
  missing_cols <- setdiff(required_cols, names(points_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in 'points_df': ", paste(missing_cols, collapse = ", "))
  }

  # Remove points with missing traits or environment
  points_df <- points_df %>%
    dplyr::filter(!is.na(summ_trait_1), !is.na(summ_trait_2), !is.na(.data[[category_col]]))

  # Check that all sample_sizes are smaller than available data
  if (any(sample_sizes > nrow(points_df))) {
    stop("One or more sample sizes exceed the number of available data points (", nrow(points_df), ").")
  }

  # Enforce minimum valid sample size to ensure meaningful binning and correlations
  min_valid_sample_size <- 30

  if (any(sample_sizes < min_valid_sample_size)) {
    stop("All sample sizes must be >= ", min_valid_sample_size,
         " to ensure valid model training and prediction. ",
         "Invalid sample sizes: ",
         paste(sample_sizes[sample_sizes < min_valid_sample_size], collapse = ", "))
  }

  # Define single iteration function
  single_iteration <- function(sample_size, iteration) {

    sampled_indices <- sample(1:nrow(points_df), size = sample_size, replace = FALSE)
    sampled_data <- points_df[sampled_indices, ]

    # Split into training/testing
    train_indices <- sample(1:nrow(sampled_data), size = floor((1 - test_split) * nrow(sampled_data)))
    test_indices <- setdiff(1:nrow(sampled_data), train_indices)

    training_data <- sampled_data[train_indices, ]
    testing_data <- sampled_data[test_indices, ]

    # Determine bin numbers if not provided
    if (is.null(grid_bins_1)) {
      grid_bins_1_train <- optimal_bins(training_data$summ_trait_1)
    }
    if (is.null(grid_bins_2)) {
      grid_bins_2_train <- optimal_bins(training_data$summ_trait_2)
    }

    # Binning function
    bin_codes <- function(traits, grid_bins) {
      trait_range <- range(traits, na.rm = TRUE)
      trait_range[1] <- trait_range[1] - 1e-6
      trait_range[2] <- trait_range[2] + 1e-6
      breaks <- seq(trait_range[1], trait_range[2], length.out = grid_bins + 1)
      codes <- .bincode(traits, breaks = breaks, include.lowest = TRUE)
      return(list(codes = codes, breaks = breaks))
    }

    # Binning for training
    mean_bin_train <- bin_codes(training_data$summ_trait_1, grid_bins_1_train)
    sd_bin_train <- bin_codes(training_data$summ_trait_2, grid_bins_2_train)
    training_data$mbc <- mean_bin_train$codes
    training_data$sdc <- sd_bin_train$codes


    # Predict category per bin
    bin_predictions <- training_data %>%
      dplyr::group_by(mbc, sdc) %>%
      dplyr::summarise(
        env_est = as.character(na.omit(get_mode(.data[[category_col]]))),
        .groups = "drop"
      )

    # Predict training
    training_data <- dplyr::left_join(training_data, bin_predictions, by = c("mbc", "sdc"))
    training_accuracy <- mean(training_data[[category_col]] == training_data$env_est, na.rm = TRUE)

    # Predict testing
    testing_data$mbc <- .bincode(testing_data$summ_trait_1, breaks = mean_bin_train$breaks, include.lowest = TRUE)
    testing_data$sdc <- .bincode(testing_data$summ_trait_2, breaks = sd_bin_train$breaks, include.lowest = TRUE)
    testing_data <- dplyr::left_join(testing_data, bin_predictions, by = c("mbc", "sdc"))
    testing_accuracy <- mean(testing_data[[category_col]] == testing_data$env_est, na.rm = TRUE)

    return(data.frame(
      SampleSize = sample_size,
      Iteration = iteration,
      Training_Accuracy = training_accuracy,
      Testing_Accuracy = testing_accuracy
    ))
  }

  # Run all iterations
  if (parallel) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, varlist = c(
      "points_df", "category_col", "test_split",
      "iterations", "single_iteration", "get_mode",
      "optimal_bins", "grid_bins_1", "grid_bins_2"
    ), envir = environment())
    parallel::clusterEvalQ(cl, library(dplyr))

    results <- parallel::parLapply(cl, sample_sizes, function(samp_size) {
      do.call(rbind, lapply(1:iterations, function(iter) {
        single_iteration(samp_size, iter)
      }))
    })

    parallel::stopCluster(cl)
  } else {
    results <- lapply(sample_sizes, function(samp_size) {
      do.call(rbind, lapply(1:iterations, function(iter) {
        single_iteration(samp_size, iter)
      }))
    })
  }

  combined_results <- do.call(rbind, results)

  # Plot
  combined_results_clean <- na.omit(combined_results)
  transp_black <- rgb(0, 0, 0, alpha = 0.3)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1, 2))

  with(combined_results_clean, {
    # Training accuracy
    plot(SampleSize, Training_Accuracy,
      pch = 16, col = transp_black,
      xlab = "Sample size", ylab = "Training Accuracy",
      ylim = c(0, 1)
    )
    loess_fit <- loess(Training_Accuracy ~ SampleSize)
    lines(sort(SampleSize), predict(loess_fit)[order(SampleSize)], lwd = 2)

    # Testing accuracy
    plot(SampleSize, Testing_Accuracy,
      pch = 16, col = transp_black,
      xlab = "Sample size", ylab = "Testing Accuracy",
      ylim = c(0, 1)
    )
    loess_fit <- loess(Testing_Accuracy ~ SampleSize)
    lines(sort(SampleSize), predict(loess_fit)[order(SampleSize)], lwd = 2)
  })

  # Summarize
  summary_results <- combined_results %>%
    dplyr::group_by(SampleSize) %>%
    dplyr::summarise(
      Training_Accuracy = mean(Training_Accuracy, na.rm = TRUE),
      Testing_Accuracy = mean(Testing_Accuracy, na.rm = TRUE)
    )

  return(list(
    combined_results = combined_results,
    summary_results = summary_results
  ))
}
