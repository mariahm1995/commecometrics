#' Run an ecometric model for quantitative environmental variables
#'
#' Builds an ecometric trait space for quantitative environmental variables,
#' estimating environmental values in binned ecometric trait space. Also
#' calculates anomalies based on observed values.
#'
#' @param points_df A data frame with columns: `mean_trait`, `sd_trait`, `count_trait`, and the environmental variable.
#' @param env_var Name of the environmental variable (e.g., "BIO12").
#' @param transform_fun Optional transformation function for environmental variable (e.g., \code{log(x + 1)}).
#' @param inv_transform_fun Optional inverse transformation for environmental variable (e.g., \code{exp(x) - 1}).
#' @param grid_bins_mean Number of bins for the mean trait axis (default: computed via Scott's rule).
#' @param grid_bins_sd Number of bins for the SD trait axis (default: computed via Scott's rule).
#' @param min_species Minimum number of species per point (default = 3).
#'
#' @return A list containing:
#' \item{points_df}{Filtered input data frame with predicted and observed environmental values, bin assignments, and anomalies.}
#' \item{eco_space}{Raster-format data frame representing trait space bins with estimated environmental values.}
#' \item{model}{Linear model object (`lm`) relating predicted environmental values to observed environmental values (transformed scale when used).}
#' \item{correlation}{Output from `cor.test`, reporting the Pearson correlation between predicted and observed environmental values (transformed scale when used).}
#' \item{diagnostics}{Summary stats about bin usage and data coverage.}
#' \item{settings}{Metadata including the modeled trait and transformation functions.}
#' \item{prediction_accuracy}{Overall percentage of correct predictions.}
#' @importFrom stats density lm cor.test
#' @export
ecometric_model <- function(points_df,
                            env_var = "BIO12",
                            transform_fun = function(x) log(x + 1),
                            inv_transform_fun = function(x) exp(x) - 1,
                            grid_bins_mean = NULL,
                            grid_bins_sd = NULL,
                            min_species = 3) {
  # Check required columns
  required_cols <- c("mean_trait", "sd_trait", "count_trait", env_var)
  missing_cols <- setdiff(required_cols, names(points_df))
  if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

  # Transform environmental variable
  if (!is.null(transform_fun)) {
    message("Applying transformation to environmental variable...")
    points_df$env_trans <- transform_fun(points_df[[env_var]])
  } else {
    message("No transformation applied to environmental variable.")
    points_df$env_trans <- points_df[[env_var]]
  }

  # Filter by species count
  message("Filtering points with at least ", min_species, " species...")
  filtered_df <- points_df %>% dplyr::filter(count_trait >= min_species)
  prop_retained <- round(nrow(filtered_df) / nrow(points_df) * 100, 1)
  message("Retained ", nrow(filtered_df), " of ", nrow(points_df), " points (", prop_retained, "%)")
  if (nrow(filtered_df) == 0) stop("No points retained. Try lowering 'min_species'.")

  # Determine bin numbers if not provided
  if (is.null(grid_bins_mean)) {
    grid_bins_mean <- optimal_bins_scott(filtered_df$mean_trait)
    message("Optimal number of mean trait bins (Scott's rule): ", grid_bins_mean)
  }
  if (is.null(grid_bins_sd)) {
    grid_bins_sd <- optimal_bins_scott(filtered_df$sd_trait)
    message("Optimal number of SD trait bins (Scott's rule): ", grid_bins_sd)
  }

  # Define bin breaks
  mrange <- range(filtered_df$mean_trait, na.rm = TRUE)
  sdrange <- range(filtered_df$sd_trait, na.rm = TRUE)
  mbrks <- seq(mrange[1] - 0.001, mrange[2] + 0.001, length.out = grid_bins_mean + 1)
  sdbrks <- seq(sdrange[1] - 0.001, sdrange[2] + 0.001, length.out = grid_bins_sd + 1)

  # Assign bin codes
  filtered_df <- filtered_df %>%
    dplyr::mutate(
      mbc = .bincode(mean_trait, breaks = mbrks),
      sdc = .bincode(sd_trait, breaks = sdbrks)
    )

  # Estimate likelihood in trait bins
  message("Estimating maximum likelihood environmental value per trait bin...")
  grid_vals <- matrix(NA, nrow = grid_bins_sd, ncol = grid_bins_mean)
  bin_counts <- matrix(0, nrow = grid_bins_sd, ncol = grid_bins_mean)

  for (i in 1:grid_bins_mean) {
    for (j in 1:grid_bins_sd) {
      idx <- which(filtered_df$mbc == i & filtered_df$sdc == j)
      bin_counts[j, i] <- length(idx)
      if (length(idx) > 0) {
        dens <- density(filtered_df$env_trans[idx], bw = 1, na.rm = TRUE)
        grid_vals[grid_bins_sd + 1 - j, i] <- dens$x[which.max(dens$y)]
      }
    }
  }

  # Raster output
  r <- raster::raster(raster::extent(0, grid_bins_mean, 0, grid_bins_sd), resolution = 1)
  r <- raster::setValues(r, as.vector(t(grid_vals)))
  raster_df <- as.data.frame(r, xy = TRUE)

  # Predict environmental values for each point
  filtered_df$env_est <- purrr::map_dbl(seq_len(nrow(filtered_df)), function(i) {
    mi <- filtered_df$mbc[i]
    si <- filtered_df$sdc[i]
    idx <- which(filtered_df$mbc == mi & filtered_df$sdc == si)
    if (length(idx) > 0) {
      dens <- density(filtered_df$env_trans[idx], bw = 1, na.rm = TRUE)
      return(dens$x[which.max(dens$y)])
    } else {
      return(NA_real_)
    }
  })

  # Anomalies
  filtered_df$env_anom <- filtered_df$env_trans - filtered_df$env_est
  if (!is.null(transform_fun) && !is.null(inv_transform_fun)) {
    filtered_df$env_anom_UN <- inv_transform_fun(filtered_df$env_anom)
    filtered_df$env_est_UN <- inv_transform_fun(filtered_df$env_est)
  } else {
    filtered_df$env_anom_UN <- NULL
    filtered_df$env_est_UN <- NULL
  }

  model <- lm(env_est ~ env_trans, data = filtered_df)
  corr <- cor.test(filtered_df$env_est, filtered_df$env_trans)

  # Bin diagnostics
  bin_counts_flipped <- bin_counts[grid_bins_sd:1, , drop = FALSE]
  used_bins <- sum(bin_counts_flipped > 0)
  total_bins <- grid_bins_mean * grid_bins_sd
  message("Used ", used_bins, " of ", total_bins, " bins (", round(used_bins / total_bins * 100, 1), "%)")

  if (is.null(transform_fun)) {
    filtered_df <- filtered_df %>% dplyr::select(-env_trans)
  }

  return(list(
    points_df = filtered_df,
    eco_space = raster_df,
    model = model,
    correlation = corr,
    diagnostics = list(
      total_points = nrow(points_df),
      retained_points = nrow(filtered_df),
      proportion_retained = prop_retained,
      used_bins = used_bins,
      total_bins = total_bins,
      bin_counts = bin_counts,
      mbrks = mbrks,
      sdbrks = sdbrks
    ),
    settings = list(
      env_var = env_var,
      transformed = !is.null(transform_fun),
      transform_fun = transform_fun,
      inv_transform_fun = inv_transform_fun
    )
  ))
}
