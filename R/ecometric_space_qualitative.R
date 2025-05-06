#' Plot Ecometric Space for Qualitative Traits with Optional Fossil Overlay
#'
#' Visualizes the predicted ecometric space (dominant category) and probability maps
#' for each category based on the output from \code{ecometric_model_qualitative()}.
#'
#' @param model_out Output from \code{ecometric_model_qualitative()}, containing environmental estimates in trait space.
#' @param category_labels Optional named vector for category labels (for prettier plots).
#' @param palette Optional color vector for categories (must match number of categories).
#' @param fossil_data Optional. Fossil data from \code{reconstruct_env_qualitative()}, including \code{fossil_mbc} and \code{fossil_sdc}.
#' @param fossil_color Outline color for fossil data bins (default = "#c44536").
#' @param modern_color Outline color for modern data bins (default: \code{"#bc4749"}).
#'
#' @return A list containing:
#'   \item{ecometric_space_plot}{ggplot showing the predicted category across trait space.}
#'   \item{probability_maps}{List of ggplots showing probability surfaces for each category.}
#'
#' @export
#'
ecometric_space_qual <- function(model_out,
                                        category_labels = NULL,
                                        palette = NULL,
                                        fossil_data = NULL,
                                        fossil_color = "black",
                                        modern_color = "#bc4749") {

  requireNamespace("ggplot2")
  requireNamespace("dplyr")
  requireNamespace("viridis")

  eco_space <- model_out$eco_space
  categories <- model_out$diagnostics$categories

  # Check if palette matches
  if (!is.null(palette)) {
    if (length(palette) != length(categories)) {
      stop("The number of colors in 'palette' must match the number of categories (", length(categories), ").")
    }
  } else {
    palette <- viridis::viridis(length(categories), direction = -1)
  }

  # Axis breakpoints
  mbreaks <- model_out$diagnostics$mbrks
  sd_breaks <- model_out$diagnostics$sdbrks
  grid_bins <- length(mbreaks) - 1

  middle_idx <- if (grid_bins %% 2 == 0) { (grid_bins / 2) + 1 } else { ceiling(grid_bins / 2) }

  x_breaks <- c(mbreaks[1], mbreaks[middle_idx], mbreaks[grid_bins - 1])
  x_labels <- round(x_breaks, 2)
  x_pos <- c(0.5, middle_idx - 0.5, grid_bins - 0.5)

  y_breaks <- c(sd_breaks[1], sd_breaks[middle_idx], sd_breaks[grid_bins - 1])
  y_labels <- round(y_breaks, 2)
  y_pos <- c(0.5, middle_idx - 0.5, grid_bins - 0.5)

  # Default category labels if not provided
  if (is.null(category_labels)) {
    category_labels <- categories
    names(category_labels) <- categories
  }

  ## 1. Predicted Category Map
  p1 <- ggplot2::ggplot(eco_space, ggplot2::aes(x = mbc - 0.5, y = sdc - 0.5, fill = as.factor(env_est))) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::scale_fill_manual(values = palette, labels = category_labels, name = "Predicted") +
    ggplot2::scale_x_continuous(name = "Mean", breaks = x_pos, labels = x_labels, expand = c(0, 0), limits = c(0, grid_bins)) +
    ggplot2::scale_y_continuous(name = "SD", breaks = y_pos, labels = y_labels, expand = c(0, 0), limits = c(0, grid_bins)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()

  # Optional fossil overlay
  if (!is.null(fossil_data)) {
    p1 <- p1 +
      ggplot2::geom_rect(data = fossil_data,
                         ggplot2::aes(xmin = as.numeric(fossil_mbc) - 1,
                                      xmax = as.numeric(fossil_mbc),
                                      ymin = as.numeric(fossil_sdc) - 1,
                                      ymax = as.numeric(fossil_sdc)),
                         inherit.aes = FALSE,
                         colour = fossil_color, alpha = 0, size = 1) +
      geom_rect(data = fossil_data,
                aes(xmin = as.numeric(mbc) - 1,
                    xmax = as.numeric(mbc),
                    ymin = as.numeric(sdc) - 1,
                    ymax = as.numeric(sdc)),
                inherit.aes = FALSE,
                colour = modern_color, alpha = 0, size = 1)
  }

  ## 2. Probability Maps for Each Category
  probability_maps <- list()

  for (cat in categories) {

    prob_col <- paste0("prob_", cat)

    if (prob_col %in% names(eco_space)) {

      rdf <- eco_space %>%
        dplyr::select(mbc, sdc, !!prob_col) %>%
        dplyr::rename(Probability = !!prob_col) %>%
        dplyr::mutate(Probability = ifelse(Probability == 0, NA_real_, Probability))

      p <- ggplot2::ggplot(rdf, ggplot2::aes(x = mbc - 0.5, y = sdc - 0.5, fill = Probability)) +
        ggplot2::geom_tile(color = NA) +
        ggplot2::scale_fill_viridis_c(name = "Probability", limits = c(0, 1), na.value = "transparent") +
        ggplot2::scale_x_continuous(name = "Mean", breaks = x_pos, labels = x_labels, expand = c(0, 0), limits = c(0, grid_bins)) +
        ggplot2::scale_y_continuous(name = "SD", breaks = y_pos, labels = y_labels, expand = c(0, 0), limits = c(0, grid_bins)) +
        ggplot2::coord_fixed() +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(category_labels[as.character(cat)])

      # Optional fossil overlay on probability maps
      if (!is.null(fossil_data)) {
        p <- p +
          ggplot2::geom_rect(data = fossil_data,
                             ggplot2::aes(xmin = as.numeric(fossil_mbc) - 1,
                                          xmax = as.numeric(fossil_mbc),
                                          ymin = as.numeric(fossil_sdc) - 1,
                                          ymax = as.numeric(fossil_sdc)),
                             inherit.aes = FALSE,
                             colour = fossil_color, alpha = 0, size = 0.8) +
          geom_rect(data = fossil_data,
                    aes(xmin = as.numeric(mbc) - 1,
                        xmax = as.numeric(mbc),
                        ymin = as.numeric(sdc) - 1,
                        ymax = as.numeric(sdc)),
                    inherit.aes = FALSE,
                    colour = modern_color, alpha = 0, size = 1)
      }

      probability_maps[[cat]] <- p
    }
  }

  return(list(
    ecometric_space_plot = p1,
    probability_maps = probability_maps
  ))
}
