#' Reconstruct past qualitative environmental categories using ecometric models
#'
#' Uses fossil community trait summaries (mean and SD) to reconstruct
#' the most probable environmental category by projecting them onto a qualitative ecometric space
#' built from modern data. Optionally, it assigns each fossil point to the nearest modern sampling point.
#'
#' @param fossildata A data frame containing fossil trait summaries per fossil site.
#'                   Must include columns for `Mean` and `SD` of the trait.
#' @param model_out Output list from \code{ecometric_model_qual()}, containing modern data, diagnostics, and model settings.
#' @param match_nearest Logical; if TRUE, matches each fossil to the nearest modern point (default = TRUE).
#' @param fossil_lon Name of the longitude column in `fossildata`. Required if \code{match_nearest = TRUE}.
#' @param fossil_lat Name of the latitude column in `fossildata`. Required if \code{match_nearest = TRUE}.
#' @param modern_id Name of the unique ID column in modern points (optional for metadata merging).
#' @param modern_lon Name of the longitude column in modern points. Required if \code{match_nearest = TRUE}.
#' @param modern_lat Name of the latitude column in modern points. Required if \code{match_nearest = TRUE}.
#' @param crs_proj Coordinate reference system for sf operations (default = EPSG:4326).
#'
#' @return A data frame (`fossildata`) updated with:
#' \describe{
#'   \item{fossil_mbc}{Assigned bin number for mean trait (based on fossil trait mean).}
#'   \item{fossil_sdc}{Assigned bin number for SD trait (based on fossil trait SD).}
#'   \item{fossil_env_est}{Predicted environmental category based on trait bin.}
#'   \item{fossil_prob_*}{Probability of each environmental category for the assigned bin.}
#'   \item{nearest_modern_point}{(Optional) ID of the nearest modern sampling point (if \code{match_nearest = TRUE}).}
#'   \item{...}{Additional columns from the matched modern site if \code{match_nearest = TRUE}.}
#' }
#' @examples
#' \dontrun{
#' # Load internal data
#' data("points", package = "commecometrics")
#' data("traits", package = "commecometrics")
#' data("polygons", package = "commecometrics")
#' data("fossils", package = "commecometrics")
#'
#' # Step 1: Summarize trait values at sampling points
#' traitsByPoint <- summarize_traits_by_point(
#'   points_df = points,
#'   trait_df = traits,
#'   species_polygons = polygons,
#'   trait_column = "RBL",
#'   species_name_col = "sci_name",
#'   continent = FALSE,
#'   parallel = FALSE
#' )
#'
#' # Step 2: Run a qualitative ecometric model (e.g., land cover class)
#' ecoModelQual <- ecometric_model_qual(
#'   points_df = traitsByPoint$points,
#'   category_col = "DOM_NUM",
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
#'   modern_id = "GlobalID",
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
  message("Binning fossil points into trait space...")

  mbrks <- model_out$diagnostics$mbrks
  sdbrks <- model_out$diagnostics$sdbrks
  eco_space <- model_out$eco_space
  modern_points <- model_out$points_df

  # Assign bins to fossils
  fossildata$fossilmbc <- .bincode(fossildata$Mean, breaks = mbrks)
  fossildata$fossilsdc <- .bincode(fossildata$SD, breaks = sdbrks)

  # Predict environmental category and probabilities for each fossil
  predictions <- vector("list", nrow(fossildata))

  for (i in seq_len(nrow(fossildata))) {
    mb <- fossildata$fossilmbc[i]
    sd <- fossildata$fossilsdc[i]

    if (!is.na(mb) && !is.na(sd)) {
      pred_row <- eco_space %>% dplyr::filter(mbc == mb, sdc == sd)

      if (nrow(pred_row) >= 1) {
        predictions[[i]] <- pred_row[1, , drop = FALSE]
      } else {
        empty <- tibble::tibble(mbc = mb, sdc = sd, env_est = NA)
        prob_cols <- grep("^prob_", names(eco_space), value = TRUE)
        for (col in prob_cols) empty[[col]] <- NA_real_
        predictions[[i]] <- empty
      }
    } else {
      empty <- tibble::tibble(mbc = NA_integer_, sdc = NA_integer_, env_est = NA)
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
      fossil_mbc = fossilmbc,
      fossil_sdc = fossilsdc,
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
