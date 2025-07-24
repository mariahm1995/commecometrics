#' Determine optimal number of bins using Scott's Rule
#'
#' Calculates the optimal number of bins for a numeric vector
#' based on Scott's rule. For more details see Scott (1979) <doi:10.1093/biomet/66.3.605>
#'
#' @param x Numeric vector.
#' @return Integer representing the optimal number of bins.
#'
#' @examples
#' # Simple example
#' # Example with normally distributed data
#' optimal_bins(rnorm(100))
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
    warning("Non-positive bin width calculated. Defaulting to 1 bin.")
    return(1)
  }
  range_x <- max(x) - min(x)
  bins <- ceiling(range_x / h)
  return(max(1, bins))
}
