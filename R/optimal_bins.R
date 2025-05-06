#' Determine optimal number of bins using Scott's Rule
#'
#' Calculates the optimal number of bins for histogramming a numeric vector
#' based on Scott's rule, which minimizes the integrated mean squared error.
#'
#' @param x Numeric vector of trait values (e.g., mean_trait or sd_trait).
#' @return Integer representing the optimal number of bins.

#' @export
optimal_bins <- function(x) {
  x <- na.omit(x)
  n <- length(x)
  if (n < 2) {
    warning("Not enough data points to compute bins.")
    return(1)
  }
  sigma <- stats::sd(x)
  h <- 3.49 * sigma / (n^(1 / 3))
  if (h <= 0) {
    warning("Non-positive bin width calculated; defaulting to 1 bin.")
    return(1)
  }
  range_x <- max(x) - min(x)
  bins <- ceiling(range_x / h)
  return(max(1, bins))
}
