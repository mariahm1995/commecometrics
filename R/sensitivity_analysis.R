#' Perform sensitivity analysis on ecometric models (quantitative environmental variables)
#'
#' This function evaluates how varying sample sizes affect the performance of ecometric models,
#' focusing on two aspects:
#' \itemize{
#'   \item \strong{Sensitivity (internal consistency)}: How accurately the model predicts environmental conditions
#'         on the same data it was trained on.
#'   \item \strong{Transferability (external applicability)}: How well the model performs on unseen data.
#' }
#' It tests different sample sizes by resampling the data multiple times (bootstrap iterations),
#' training an ecometric model on each subset, and evaluating prediction error and correlation.
#'
#' Four base R plots are generated to visualize model performance as a function of sample size:
#' \enumerate{
#'   \item \strong{Training correlation vs. Sample size:} Shows how well the model fits training data.
#'   \item \strong{Testing correlation vs. Sample size:} Shows generalizability to new data.
#'   \item \strong{Training mean anomaly vs. Sample size:} Shows average prediction error on training data.
#'   \item \strong{Testing mean anomaly vs. Sample size:} Shows average prediction error on test data.
#' }
#'
#' Parallel processing is supported to speed up the analysis.
#'
#' @param points_df Output first element of the list from \code{summarize_traits_by_point()}. A data frame with columns: `mean_trait`, `sd_trait`, `count_trait`, and the environmental variable.
#' @param env_var Name of the environmental variable column in points_df (e.g., "BIO12").
#' @param sample_sizes Vector of sample sizes to evaluate (default: seq(100, 10000, 1000)).
#' @param iterations Number of bootstrap iterations per sample size (default: 20).
#' @param test_split Proportion of data to use for testing (default: 0.2).
#' @param grid_bins_mean Number of bins for the mean trait axis. If `NULL` (default),
#'   the number is calculated automatically using Scott's rule via `optimal_bins()`.
#' @param grid_bins_sd Number of bins for the SD trait axis. If `NULL` (default),
#'   the number is calculated automatically using Scott's rule via `optimal_bins()`.
#' @param transform_fun Function to transform the environmental variable (default: NULL = no transformation).
#' @param parallel Logical; whether to use parallel processing (default: TRUE).
#' @param n_cores Number of cores to use for parallel processing (default: parallel::detectCores() - 1).
#'
#' @return A list containing:
#'   \item{combined_results}{A data frame with mean absolute anomalies and correlations for each sample size and iteration.}
#'   \item{summary_results}{A data frame summarizing the mean anomalies and correlations across sample sizes.}
#'
#' @examples
#' \dontrun{
#' # Load internal data
#' data("points", package = "commecometrics")
#' data("traits", package = "commecometrics")
#' data("polygons", package = "commecometrics")
#'
#' # Summarize trait values at sampling points
#' traitsByPoint <- summarize_traits_by_point(
#'   points_df = points,
#'   trait_df = traits,
#'   species_polygons = polygons,
#'   trait_column = "RBL",
#'   species_name_col = "sci_name",
#'   continent = FALSE,
#'   parallel = FALSE
#' )
#'
#' # Run sensitivity analysis using annual precipitation (BIO12)
#' sensitivityResults <- sensitivity_analysis(
#'   points_df = traitsByPoint$points,
#'   env_var = "BIO12",
#'   sample_sizes = seq(100, 1000, 300),
#'   iterations = 5,
#'   transform_fun = function(x) log(x + 1),
#'   parallel = FALSE  # Set to TRUE for faster performance on multicore machines
#' )
#'
#' # View results
#' head(sensitivityResults$summary_results)
#' }
#' @export
sensitivity_analysis <- function(points_df,
                                 env_var,
                                 sample_sizes = seq(100, 10000, 1000),
                                 iterations = 20,
                                 test_split = 0.2,
                                 grid_bins_mean = NULL,
                                 grid_bins_sd = NULL,
                                 transform_fun = NULL,
                                 parallel = TRUE,
                                 n_cores = parallel::detectCores() - 1) {
  # Load required package
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required but not installed.")
  }

  # Check for required columns
  required_cols <- c("mean_trait", "sd_trait", "count_trait", env_var)
  missing_cols <- setdiff(required_cols, names(points_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in 'points_df': ", paste(missing_cols, collapse = ", "))
  }

  # Remove points with missing traits or environment
  points_df <- points_df %>%
    dplyr::filter(!is.na(mean_trait) & !is.na(sd_trait)) %>%
    dplyr::filter(!is.na(.data[[env_var]]))

  # Check that all sample_sizes are smaller than available data
  if (any(sample_sizes > nrow(points_df))) {
    stop("One or more sample sizes exceed the number of available data points (", nrow(points_df), ").")
  }

  # Transform environmental variable if needed
  if (!is.null(transform_fun)) {
    points_df$env_trans <- transform_fun(points_df[[env_var]])
  } else {
    points_df$env_trans <- points_df[[env_var]]
  }

  # Define single iteration function
  single_iteration <- function(sample_size, iteration) {
    set.seed(iteration)

    sampled_indices <- sample(1:nrow(points_df), size = sample_size, replace = FALSE)
    sampled_data <- points_df[sampled_indices, ]

    # Split into training and testing
    train_indices <- sample(1:nrow(sampled_data), size = floor((1 - test_split) * nrow(sampled_data)))
    test_indices <- setdiff(1:nrow(sampled_data), train_indices)

    training_data <- sampled_data[train_indices, ]
    testing_data <- sampled_data[test_indices, ]

    # Determine bin numbers if not provided
    if (is.null(grid_bins_mean)) {
      grid_bins_mean_train <- optimal_bins(training_data$mean_trait)
    }
    if (is.null(grid_bins_sd)) {
      grid_bins_sd_train <- optimal_bins(training_data$sd_trait)
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
    mean_bin_train <- bin_codes(training_data$mean_trait, grid_bins_mean_train)
    sd_bin_train <- bin_codes(training_data$sd_trait, grid_bins_sd_train)
    training_data$mbc <- mean_bin_train$codes
    training_data$sdc <- sd_bin_train$codes


    # Estimate environment per bin
    bin_env_estimates <- matrix(NA, nrow = grid_bins_mean_train, ncol = grid_bins_sd_train)
    for (i in 1:grid_bins_mean_train) {
      for (j in 1:grid_bins_sd_train) {
        bin_data <- training_data$env_trans[training_data$mbc == i & training_data$sdc == j]
        if (length(bin_data) > 0) {
          dens <- density(bin_data, bw = 1, na.rm = TRUE)
          bin_env_estimates[i, j] <- dens$x[which.max(dens$y)]
        }
      }
    }

    # Predict for training data
    training_preds <- mapply(function(m, s) {
      if (!is.na(m) && !is.na(s) && !is.na(bin_env_estimates[m, s])) {
        return(bin_env_estimates[m, s])
      } else {
        return(NA)
      }
    }, training_data$mbc, training_data$sdc)

    training_anom <- abs(training_data$env_trans - training_preds)
    training_cor <- cor(training_data$env_trans, training_preds, use = "complete.obs")

    # Predict for testing data
    testing_data$mbc <- .bincode(testing_data$mean_trait, breaks = mean_bin_train$breaks, include.lowest = TRUE)
    testing_data$sdc <- .bincode(testing_data$sd_trait, breaks = sd_bin_train$breaks, include.lowest = TRUE)

    testing_preds <- mapply(function(m, s) {
      if (!is.na(m) && !is.na(s) && !is.na(bin_env_estimates[m, s])) {
        return(bin_env_estimates[m, s])
      } else {
        return(NA)
      }
    }, testing_data$mbc, testing_data$sdc)

    testing_anom <- abs(testing_data$env_trans - testing_preds)
    testing_cor <- cor(testing_data$env_trans, testing_preds, use = "complete.obs")

    return(data.frame(
      SampleSize = sample_size,
      Iteration = iteration,
      Training_Mean_Anomaly = mean(training_anom, na.rm = TRUE),
      Training_Correlation = training_cor,
      Testing_Mean_Anomaly = mean(testing_anom, na.rm = TRUE),
      Testing_Correlation = testing_cor
    ))
  }

  # Perform iterations
  if (parallel) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, varlist = c(
      "points_df", "env_var", "transform_fun", "test_split",
      "iterations", "single_iteration", "optimal_bins"
    ), envir = environment())
    parallel::clusterEvalQ(cl, library(stats))

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

  # Combine results
  combined_results <- do.call(rbind, results)

  # Plotting
  combined_results_clean <- na.omit(combined_results)
  transp_black <- rgb(0, 0, 0, alpha = 0.3)
  par(mfrow = c(2, 2))

  with(combined_results_clean, {
    plot(SampleSize, Training_Correlation,
      pch = 16, col = transp_black,
      xlab = "Sample size", ylab = "Training correlation",
      main = "Training correlation vs Sample size"
    )
    loess_fit <- loess(Training_Correlation ~ SampleSize)
    lines(sort(SampleSize), predict(loess_fit)[order(SampleSize)], lwd = 2)

    plot(SampleSize, Testing_Correlation,
      pch = 16, col = transp_black,
      xlab = "Sample size", ylab = "Testing correlation",
      main = "Testing correlation vs Sample size"
    )
    loess_fit <- loess(Testing_Correlation ~ SampleSize)
    lines(sort(SampleSize), predict(loess_fit)[order(SampleSize)], lwd = 2)

    plot(SampleSize, Training_Mean_Anomaly,
      pch = 16, col = transp_black,
      xlab = "Sample size", ylab = "Training mean anomaly",
      main = "Training anomaly vs Sample size"
    )
    loess_fit <- loess(Training_Mean_Anomaly ~ SampleSize)
    lines(sort(SampleSize), predict(loess_fit)[order(SampleSize)], lwd = 2)

    plot(SampleSize, Testing_Mean_Anomaly,
      pch = 16, col = transp_black,
      xlab = "Sample size", ylab = "Testing mean anomaly",
      main = "Testing anomaly vs Sample size"
    )
    loess_fit <- loess(Testing_Mean_Anomaly ~ SampleSize)
    lines(sort(SampleSize), predict(loess_fit)[order(SampleSize)], lwd = 2)
  })

  # Summarizing
  summary_results <- combined_results %>%
    dplyr::group_by(SampleSize) %>%
    dplyr::summarise(
      Training_Mean_Anomaly = mean(Training_Mean_Anomaly, na.rm = TRUE),
      Training_Correlation = mean(Training_Correlation, na.rm = TRUE),
      Testing_Mean_Anomaly = mean(Testing_Mean_Anomaly, na.rm = TRUE),
      Testing_Correlation = mean(Testing_Correlation, na.rm = TRUE)
    )

  return(list(
    combined_results = combined_results,
    summary_results = summary_results
  ))
}
