#' @importFrom grDevices rgb
#' @importFrom graphics par
#' @importFrom stats cor na.omit sd
#' @importFrom ggplot2 ggplot aes geom_raster geom_rect scale_fill_gradientn
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous coord_fixed theme_bw
NULL
utils::globalVariables(c(
".data", "env_est", "envest", "envestUN", "layer", "x", "y", "fossil_mbc", "fossil_sdc",
"mbc", "sdc", "Probability", "SampleSize", "Training_Mean_Anomaly", "Training_Correlation",
"Testing_Mean_Anomaly", "Testing_Correlation", "Training_Accuracy", "Testing_Accuracy",
"fossilmbc", "fossilsdc", "fossilsdbc", "near.point", "maxlimit", "minlimit",
"maxlimitUN", "minlimitUN", "TaxonName", "count_trait", "mean_trait", "sd_trait", "env_trans"
))
