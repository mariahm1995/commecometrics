#' Run an ecometric model for quantitative environmental variables
#'
#' Builds an ecometric trait space for quantitative environmental variables,
#' estimating environmental values of each category at each trait bin combination.
#' Also calculates anomalies based on observed values for each point.
#'
#' @param points_df Output first element of the list from \code{summarize_traits_by_point()}. A data frame with columns: `summ_trait_1`, `summ_trait_2`, `count_trait`, and the environmental variable.
#' @param env_var Name of the environmental variable (e.g., "precip”).
#' @param transform_fun Optional transformation function for environmental variable (e.g., \code{log(x + 1)}).
#' @param inv_transform_fun Optional inverse transformation for environmental variable (e.g., \code{exp(x) - 1}).
#' @param grid_bins_1 Number of bins for the first trait axis. If `NULL` (default),
#'   the number is calculated automatically using Scott's rule via `optimal_bins()`.
#' @param grid_bins_2 Number of bins for the second trait axis. If `NULL` (default),
#'   the number is calculated automatically using Scott's rule via `optimal_bins()`.
#' @param min_species Minimum number of species per point (default = 3).
#'
#' @return A list containing:
#' \item{points_df}{Filtered input data frame with the following added columns:
#'   \describe{
#'     \item{env_trans}{Transformed environmental variable (if a transformation function is used).}
#'     \item{bin_1}{Bin assignment code for mean trait value.}
#'     \item{bin_2}{Bin assignment code for standard deviation of trait.}
#'     \item{env_est}{Predicted (maximum likelihood) environmental value on transformed scale.}
#'     \item{env_anom}{Difference between observed and predicted environmental values (transformed scale).}
#'     \item{env_est_UN}{Inverse-transformed predicted value (if `inv_transform_fun` is provided).}
#'     \item{env_anom_UN}{Inverse-transformed anomaly value (if `inv_transform_fun` is provided).}
#'   }
#' }
#' \item{eco_space}{Raster-format data frame representing trait space bins with estimated environmental values.}
#' \item{model}{Linear model object (`lm`) relating predicted environmental values to observed environmental values (transformed scale when used).}
#' \item{correlation}{Output from `cor.test`, reporting the Pearson correlation between predicted and observed environmental values (transformed scale when used).}
#' \item{diagnostics}{Summary stats about bin usage and data coverage.}
#' \item{settings}{Metadata including the modeled trait and transformation functions.}
#' \item{prediction_accuracy}{Overall percentage of correct predictions.}
#' @importFrom stats density lm cor.test
#' @examples
#' \dontrun{
#' # Load internal dataset
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
#'   species_name_col = "TaxonName",
#'   continent = FALSE,
#'   parallel = FALSE
#' )
#'
#' # Fit an ecometric model using annual precipitation (BIO12)
#' modelResult <- ecometric_model(
#'   points_df = traitsByPoint$points,
#'   env_var = "precip",
#'   transform_fun = function(x) log(x + 1),
#'   inv_transform_fun = function(x) exp(x) - 1,
#'   min_species = 3
#' )
#'
#' # View correlation between predicted and observed values
#' print(modelResult$correlation)
#'
#' # View summary of the linear model fit
#' summary(modelResult$model)
#' }
#' @export
ecometric_model <- function(points_df,
                            env_var = "env_var",
                            transform_fun = function(x) x,
                            inv_transform_fun = function(x) x,
                            grid_bins_1 = NULL,
                            grid_bins_2 = NULL,
                            min_species = 3) {
  # Check required columns
  required_cols <- c("summ_trait_1", "summ_trait_2", "count_trait", env_var)
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

  # Assign bin codes
  filtered_df <- filtered_df %>%
    dplyr::mutate(
      bin_1 = .bincode(summ_trait_1, breaks = mbrks),
      bin_2 = .bincode(summ_trait_2, breaks = sdbrks)
    )

  # Estimate likelihood in trait bins
  message("Estimating maximum likelihood environmental value per trait bin...")
  grid_vals <- matrix(NA, nrow = grid_bins_2, ncol = grid_bins_1)
  bin_counts <- matrix(0, nrow = grid_bins_2, ncol = grid_bins_1)

  for (i in 1:grid_bins_1) {
    for (j in 1:grid_bins_2) {
      idx <- which(filtered_df$bin_1 == i & filtered_df$bin_2 == j)
      bin_counts[j, i] <- length(idx)
      if (length(idx) > 0) {
        dens <- density(filtered_df$env_trans[idx], bw = 1, na.rm = TRUE)
        grid_vals[grid_bins_2 + 1 - j, i] <- dens$x[which.max(dens$y)]
      }
    }
  }

  # Raster output
  r <- raster::raster(raster::extent(0, grid_bins_1, 0, grid_bins_2), resolution = 1)
  r <- raster::setValues(r, as.vector(t(grid_vals)))
  raster_df <- raster::as.data.frame(r, xy = TRUE)

  # Predict environmental values for each point
  filtered_df$env_est <- purrr::map_dbl(seq_len(nrow(filtered_df)), function(i) {
    mi <- filtered_df$bin_1[i]
    si <- filtered_df$bin_2[i]
    idx <- which(filtered_df$bin_1 == mi & filtered_df$bin_2 == si)
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
  bin_counts_flipped <- bin_counts[grid_bins_2:1, , drop = FALSE]
  used_bins <- sum(bin_counts_flipped > 0)
  total_bins <- grid_bins_1 * grid_bins_2
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
      brks_1 = mbrks,
      brks_2 = sdbrks
    ),
    settings = list(
      env_var = env_var,
      transformed = !is.null(transform_fun),
      transform_fun = transform_fun,
      inv_transform_fun = inv_transform_fun
    )
  ))
}
