#' @title Internal utilities and global variables
#' @description Internal utilities and variable declarations to support NSE and ggplot2 piping.
#' @name commecometrics-utils
#' @importFrom grDevices rgb
#' @importFrom graphics par
#' @importFrom stats cor na.omit sd
#' @importFrom ggplot2 ggplot aes geom_raster geom_rect scale_fill_gradientn
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous coord_fixed theme_bw
NULL

# Suppress CMD check notes for non-standard evaluation variables used in pipes and ggplot2
utils::globalVariables(c(
  ".data", "env_est", "envest", "envestUN", "layer", "x", "y", "fossil_bin_1", "fossil_bin_2",
  "bin_1", "bin_2", "Probability", "SampleSize", "Training_Mean_Anomaly", "Training_Correlation",
  "Testing_Mean_Anomaly", "Testing_Correlation", "Training_Accuracy", "Testing_Accuracy",
  "fossilmbc", "fossil_summ_trait_2", "fossil_summ_trait_1", "near.point", "maxlimit", "minlimit", "brks_1", "brks_2",
  "maxlimitUN", "minlimitUN", "TaxonName", "count_trait", "summ_trait_1", "summ_trait_2", "env_trans",
  "fossilsdbc", "fossilsdc", "mbc", "points_df", "sdc"
))
