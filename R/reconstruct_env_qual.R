#' Reconstruct past qualitative environmental categories using ecometric models
#'
#' Uses fossil community trait summaries to reconstruct the most likely
#' environmental category by projecting them onto a qualitative ecometric space
#' built from modern data. Optionally, it assigns each fossil point to the nearest modern sampling point.
#'
#' @param fossildata A data frame containing fossil trait summaries per fossil site.
#'                   Must include columns corresponding to the same two summary metrics used for modern communities,
#'                   using the column names specified by `fossil_summ_trait_1` and `fossil_summ_trait_2`.
#' @param model_out Output list from \code{ecometric_model_qual()}, containing modern data, diagnostics, and model settings.
#' @param match_nearest Logical; if TRUE, matches each fossil to the nearest modern point (default = TRUE).
#' @param fossil_lon Name of the longitude column in `fossildata`. Required if \code{match_nearest = TRUE}.
#' @param fossil_lat Name of the latitude column in `fossildata`. Required if \code{match_nearest = TRUE}.
#' @param modern_id Name of the unique ID column in modern points (e.g., "GlobalID").
#' @param modern_lon Name of the longitude column in modern points. Required if \code{match_nearest = TRUE}.
#' @param modern_lat Name of the latitude column in modern points. Required if \code{match_nearest = TRUE}.
#' @param crs_proj Coordinate reference system to use when converting fossil and modern data to sf format (default = EPSG:4326)
#'
#' @return A data frame (`fossildata`) with reconstructed environmental values and optional nearest modern point data. Includes the following additional columns:
#' \describe{
#'   \item{fossil_bin_1}{Numeric bin index for the first trait axis (based on first summary metric of trait distribution of fossil communities).}
#'   \item{fossil_bin_2}{Numeric bin index for the second trait axis (based on second summary metric of trait distribution of fossil communities).}
#'   \item{fossil_env_est}{Predicted environmental category based on trait bin.}
#'   \item{fossil_prob_*}{Probability of each environmental category for the assigned bin.}
#'   \item{nearest_modern_point}{(Optional) ID of the nearest modern sampling point (if \code{match_nearest = TRUE}).}
#'   \item{...}{Additional columns from the matched modern site if \code{match_nearest = TRUE}.}
#' }
#' @examples
#' \donttest{
#' # Load internal data
#' data("geoPoints", package = "commecometrics")
#' data("traits", package = "commecometrics")
#' data("spRanges", package = "commecometrics")
#' data("fossils", package = "commecometrics")
#'
#' # Step 1: Summarize trait values at sampling points
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
#' # Step 2: Run a qualitative ecometric model (e.g., land cover class)
#' ecoModelQual <- ecometric_model_qual(
#'   points_df = traitsByPoint$points,
#'   category_col = "vegetation",
#'   min_species = 3
#' )
#'
#' # Step 3: Reconstruct qualitative environments for fossil data
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
#' }
#' @export
reconstruct_env_qual <- function(fossildata,
                                 model_out,
                                 match_nearest = TRUE,
                                 fossil_lon = NULL,
                                 fossil_lat = NULL,
                                 modern_id = NULL,
                                 modern_lon = NULL,
                                 modern_lat = NULL,
                                 crs_proj = 4326) {
  # Check for required columns
  required_cols <- c("fossil_summ_trait_1", "fossil_summ_trait_2")
  missing_cols <- setdiff(required_cols, names(fossildata))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in 'fossildata': ", paste(missing_cols, collapse = ", "))
  }

  message("Binning fossil points into trait space...")

  mbrks <- model_out$diagnostics$brks_1
  sdbrks <- model_out$diagnostics$brks_2
  eco_space <- model_out$eco_space
  modern_points <- model_out$points_df

  # Assign bins to fossils
  fossildata$fossilmbc <- .bincode(fossildata$fossil_summ_trait_1, breaks = mbrks)
  fossildata$fossilsdc <- .bincode(fossildata$fossil_summ_trait_2, breaks = sdbrks)

  # Predict environmental category and probabilities for each fossil
  predictions <- vector("list", nrow(fossildata))

  for (i in seq_len(nrow(fossildata))) {
    mb <- fossildata$fossilmbc[i]
    sd <- fossildata$fossilsdc[i]

    if (!is.na(mb) && !is.na(sd)) {
      pred_row <- eco_space %>% dplyr::filter(x == mb, y == sd)

      if (nrow(pred_row) >= 1) {
        predictions[[i]] <- pred_row[1, , drop = FALSE]
      } else {
        empty <- tibble::tibble(x = mb, y = sd, env_est = NA)
        prob_cols <- grep("^prob_", names(eco_space), value = TRUE)
        for (col in prob_cols) empty[[col]] <- NA_real_
        predictions[[i]] <- empty
      }
    } else {
      empty <- tibble::tibble(x = NA_integer_, y = NA_integer_, env_est = NA)
      prob_cols <- grep("^prob_", names(eco_space), value = TRUE)
      for (col in prob_cols) empty[[col]] <- NA_real_
      predictions[[i]] <- empty
    }
  }

  fossil_predictions <- dplyr::bind_rows(predictions)

  # Combine fossil data with predictions
  fossildata <- dplyr::bind_cols(fossildata, fossil_predictions[, -(1:2)])

  # Rename for clarity
  fossildata <- fossildata %>%
    dplyr::rename(
      fossil_bin_1 = fossilmbc,
      fossil_bin_2 = fossilsdc,
      fossil_env_est = env_est
    )

  # Rename probability columns
  prob_cols <- grep("^prob_", names(fossildata), value = TRUE)
  if (length(prob_cols) > 0) {
    names(fossildata)[names(fossildata) %in% prob_cols] <- paste0("fossil_", prob_cols)
  }

  # Match to nearest modern point (optional)
  if (match_nearest) {
    if (is.null(fossil_lon) || is.null(fossil_lat) || is.null(modern_id) ||
        is.null(modern_lon) || is.null(modern_lat)) {
      stop("When matching nearest, please provide fossil_lon, fossil_lat, modern_id, modern_lon, and modern_lat.")
    }

    fossil.sf <- sf::st_as_sf(fossildata, coords = c(fossil_lon, fossil_lat), crs = crs_proj)
    modern.sf <- sf::st_as_sf(modern_points, coords = c(modern_lon, modern_lat), crs = crs_proj)

    nearest_idx <- purrr::map_int(seq_len(nrow(fossil.sf)), function(i) {
      which.min(sf::st_distance(fossil.sf[i, ], modern.sf))
    })

    fossildata$nearest_modern_point <- modern_points[[modern_id]][nearest_idx]

    fossildata <- merge(fossildata, modern_points, by.x = "nearest_modern_point", by.y = modern_id)
  }

  return(fossildata)
}
