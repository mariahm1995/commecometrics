#' Plot ecometric space for quantitative environmental variables
#'
#' Visualizes the ecometric space for quantitative environmental variables
#' based on the output from \code{ecometric_model()}.
#'
#' @param model_out Output from \code{ecometric_model()}, containing environmental estimates in trait space.
#' @param env_name Name to display for the environmental variable (used in the legend title).
#' @param fossil_data Optional. Output from \code{reconstruct_env()}.
#' @param fossil_color Outline color for fossil data bins (default: \code{"#000000"}).
#' @param modern_color Outline color for modern data bins (default: \code{"#bc4749"}).
#' @param palette Vector of colors to use for the gradient scale representing environmental values.
#' @param x_label Label for the x-axis in the output plots (default: "Summary metric 1").
#' @param y_label Label for the y-axis in the output plots (default: "Summary metric 2").
#'
#' @return A ggplot2 object visualizing the ecometric trait-environment surface.
#' @examples
#' \dontrun{
#' # Load internal data
#' data("geoPoints", package = "commecometrics")
#' data("traits", package = "commecometrics")
#' data("spRanges", package = "commecometrics")
#' data("fossils", package = "commecometrics")
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
#' # Run ecometric model
#' ecoModel <- ecometric_model(
#'   points_df = traitsByPoint$points,
#'   env_var = "precip",
#'   transform_fun = function(x) log(x + 1),
#'   inv_transform_fun = function(x) exp(x) - 1,
#'   min_species = 3
#' )
#'
#' # Reconstruct environments for fossil sites
#' recon <- reconstruct_env(
#'   fossildata = fossils,
#'   model_out = ecoModel,
#'   match_nearest = TRUE,
#'   fossil_lon = "Long",
#'   fossil_lat = "Lat",
#'   modern_id = "ID",
#'   modern_lon = "Longitude",
#'   modern_lat = "Latitude"
#' )
#'
#' # Plot the ecometric trait–environment space
#' ecometricPlot <- ecometric_space(
#'   model_out = ecoModel,
#'   env_name = "Precipitation (log mm)",
#'   fossil_data = recon
#' )
#'
#' # Display plot
#' print(ecometricPlot)
#' }
#' @export

ecometric_space <- function(model_out,
                            env_name = "env_var",
                            fossil_data = NULL,
                            fossil_color = "#000000",
                            modern_color = "#bc4749",
                            palette = c("#bc6c25", "#fefae0", "#606c38"),
                            x_label = "Summary metric 1",
                            y_label = "Summary metric 2") {
  # Extract model outputs
  raster_df <- model_out$eco_space
  mbreaks <- model_out$diagnostics$brks_1
  sd_breaks <- model_out$diagnostics$brks_2
  grid_bins_x <- length(mbreaks) - 1
  grid_bins_y <- length(sd_breaks) - 1

  middle_idx_x <- if (grid_bins_x %% 2 == 0) {
    (grid_bins_x / 2) + 1
  } else {
    ceiling(grid_bins_x / 2)
  }

  middle_idx_y <- if (grid_bins_y %% 2 == 0) {
    (grid_bins_y / 2) + 1
  } else {
    ceiling(grid_bins_y / 2)
  }

  # X axis labels and breaks
  x_breaks <- c(mbreaks[1], mbreaks[middle_idx_x], mbreaks[grid_bins_x])
  x_labels <- round(x_breaks, 2)
  x_pos <- c(0.5, middle_idx_x - 0.5, grid_bins_x - 0.5)

  # Y axis labels and breaks
  y_breaks <- c(sd_breaks[1], sd_breaks[middle_idx_y], sd_breaks[grid_bins_y])
  y_labels <- round(y_breaks, 2)
  y_pos <- c(0.5, middle_idx_y - 0.5, grid_bins_y - 0.5)

  # Base plot
  ecospace <- ggplot(raster_df, aes(x = x, y = y, fill = layer)) +
    geom_raster() +
    scale_fill_gradientn(colors = palette, name = env_name, na.value = "transparent") +
    scale_x_continuous(name = x_label, breaks = x_pos, labels = x_labels, expand = c(0, 0), limits = c(0, grid_bins_x)) +
    scale_y_continuous(name = y_label, breaks = y_pos, labels = y_labels, expand = c(0, 0), limits = c(0, grid_bins_y)) +
    coord_fixed() +
    theme_bw()

  # Optional fossil overlay
  if (!is.null(fossil_data)) {
    ecospace <- ecospace +
      geom_rect(
        data = fossil_data,
        aes(
          xmin = as.numeric(fossil_bin_1) - 1,
          xmax = as.numeric(fossil_bin_1),
          ymin = as.numeric(fossil_bin_2) - 1,
          ymax = as.numeric(fossil_bin_2)
        ),
        inherit.aes = FALSE,
        colour = fossil_color, alpha = 0, linewidth = 1
      ) +
      geom_rect(
        data = fossil_data,
        aes(
          xmin = as.numeric(bin_1) - 1,
          xmax = as.numeric(bin_1),
          ymin = as.numeric(bin_2) - 1,
          ymax = as.numeric(bin_2)
        ),
        inherit.aes = FALSE,
        colour = modern_color, alpha = 0, linewidth = 1
      )
  }

  return(ecospace)
}
