#' Plot ecometric space for qualitative environmental variables
#'
#' Visualizes the predicted ecometric space (predicted category) and probability maps
#' for each category based on the output from \code{ecometric_model_qualitative()}.
#'
#' @param model_out Output from \code{ecometric_model_qualitative()}, containing environmental estimates in trait space.
#' @param category_labels Optional named vector for category labels (used in the legend title). If \code{NULL}, the unique strings in the predicted category column (\code{env_est}) will be used as-is.
#' @param palette Optional color vector for categories (must match number of categories).
#' @param fossil_data Optional. Output from \code{reconstruct_env_qual()}.
#' @param fossil_color Outline color for fossil data bins (default = "#000000").
#' @param modern_color Outline color for modern data bins (default: \code{"#bc4749"}).
#'
#' @return A list containing:
#'   \item{ecometric_space_plot}{ggplot showing the predicted category across trait space.}
#'   \item{probability_maps}{List of ggplots showing probability surfaces for each category.}
#'
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
#' # Run ecometric model for qualitative variable
#' ecoModelQual <- ecometric_model_qual(
#'   points_df = traitsByPoint$points,
#'   category_col = "vegetation",
#'   min_species = 3
#' )
#'
#' # Reconstruct fossil environmental categories
#' reconQual <- reconstruct_env_qual(
#'   fossildata = fossils,
#'   model_out = ecoModelQual,
#'   match_nearest = TRUE,
#'   fossil_lon = "Long",
#'   fossil_lat = "Lat",
#'   modern_id = "ID",
#'   modern_lon = "Longitude",
#'   modern_lat = "Latitude"
#' )
#'
#' # Plot qualitative ecometric space
#' ecoPlotQual <- ecometric_space_qual(
#'   model_out = ecoModelQual,
#'   fossil_data = reconQual
#' )
#'
#' # Display predicted category map
#' print(ecoPlotQual$ecometric_space_plot)
#'
#' # Display one of the probability maps
#' print(ecoPlotQual$probability_maps[["1"]])
#' }
#' @export
#'
ecometric_space_qual <- function(model_out,
                                 category_labels = NULL,
                                 palette = NULL,
                                 fossil_data = NULL,
                                 fossil_color = "#000000",
                                 modern_color = "#bc4749") {
  requireNamespace("ggplot2")
  requireNamespace("dplyr")
  requireNamespace("viridis")

  eco_space <- model_out$eco_space

  all_prob_cols <- grep("^prob_", names(eco_space), value = TRUE)
  categories <- gsub("^prob_", "", all_prob_cols)

  # Check if palette matches
  if (!is.null(palette)) {
    if (length(palette) != length(categories)) {
      stop("The number of colors in 'palette' must match the number of categories (", length(categories), ").")
    }
  } else {
    palette <- viridis::viridis(length(categories), direction = -1)
  }

  # Axis breakpoints
  mbreaks <- model_out$diagnostics$brks_1
  sd_breaks <- model_out$diagnostics$brks_2
  grid_bins_x <- length(mbreaks) - 1
  grid_bins_y <- length(sd_breaks) - 1

  # Axis midpoints
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

  # Default category labels if not provided
  if (is.null(category_labels)) {
    category_labels <- categories
    names(category_labels) <- categories
  }

  ## 1. Predicted Category Map
  p1 <- ggplot2::ggplot(eco_space, ggplot2::aes(x = bin_1 - 0.5, y = bin_2 - 0.5, fill = as.factor(env_est))) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::scale_fill_manual(values = palette, labels = category_labels, name = "Predicted") +
    ggplot2::scale_x_continuous(name = "Summary metric 1", breaks = x_pos, labels = x_labels, expand = c(0, 0), limits = c(0, grid_bins_x)) +
    ggplot2::scale_y_continuous(name = "Summary metric 2", breaks = y_pos, labels = y_labels, expand = c(0, 0), limits = c(0, grid_bins_y)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()

  # Optional fossil overlay
  if (!is.null(fossil_data)) {
    p1 <- p1 +
      ggplot2::geom_rect(
        data = fossil_data,
        ggplot2::aes(
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

  ## 2. Probability Maps for Each Category
  probability_maps <- list()

  for (cat in categories) {
    prob_col <- paste0("prob_", cat)

    if (prob_col %in% names(eco_space)) {
      rdf <- eco_space %>%
        dplyr::select(bin_1, bin_2, !!prob_col) %>%
        dplyr::rename(Probability = !!prob_col) %>%
        dplyr::mutate(Probability = ifelse(Probability == 0, NA_real_, Probability))

      p <- ggplot2::ggplot(rdf, ggplot2::aes(x = bin_1 - 0.5, y = bin_2 - 0.5, fill = Probability)) +
        ggplot2::geom_tile(color = NA) +
        ggplot2::scale_fill_viridis_c(name = "Probability", limits = c(0, 1), na.value = "transparent") +
        ggplot2::scale_x_continuous(name = "Summary metric 1", breaks = x_pos, labels = x_labels, expand = c(0, 0), limits = c(0, grid_bins_x)) +
        ggplot2::scale_y_continuous(name = "Summary metric 2", breaks = y_pos, labels = y_labels, expand = c(0, 0), limits = c(0, grid_bins_y)) +
        ggplot2::coord_fixed() +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(category_labels[as.character(cat)])

      # Optional fossil overlay on probability maps
      if (!is.null(fossil_data)) {
        p <- p +
          ggplot2::geom_rect(
            data = fossil_data,
            ggplot2::aes(
              xmin = as.numeric(fossil_bin_1) - 1,
              xmax = as.numeric(fossil_bin_1),
              ymin = as.numeric(fossil_bin_2) - 1,
              ymax = as.numeric(fossil_bin_2)
            ),
            inherit.aes = FALSE,
            colour = fossil_color, alpha = 0, size = 0.8
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
            colour = modern_color, alpha = 0, size = 1
          )
      }

      probability_maps[[cat]] <- p
    }
  }

  return(list(
    ecometric_space_plot = p1,
    probability_maps = probability_maps
  ))
}
