#' Reconstruct Past Qualitative Environmental Categories Using Ecometric Trait Bins
#'
#' This function uses fossil community trait summaries (mean and SD) to reconstruct
#' the most probable environmental category by projecting them onto a qualitative ecometric space
#' built from modern data. Optionally, it assigns each fossil point to the nearest modern sampling point.
#'
#' @param fossildata A data frame containing fossil trait summaries (must include columns `Mean` and `SD`).
#' @param model_out Output from \code{ecometric_model_qualitative()}.
#' @param match_nearest Logical; if TRUE, matches each fossil to the nearest modern point (default = TRUE).
#' @param fossil_lon Name of the longitude column in `fossildata`. Required if \code{match_nearest = TRUE}.
#' @param fossil_lat Name of the latitude column in `fossildata`. Required if \code{match_nearest = TRUE}.
#' @param modern_id Name of the unique ID column in modern points (optional for metadata merging).
#' @param modern_lon Name of the longitude column in modern points. Required if \code{match_nearest = TRUE}.
#' @param modern_lat Name of the latitude column in modern points. Required if \code{match_nearest = TRUE}.
#' @param crs_proj Coordinate reference system for sf operations (default = EPSG:4326).
#'
#' @return A data frame (`fossildata`) updated with:
#' \itemize{
#'   \item Fossil bin assignments (fossil_mbc, fossil_sdc).
#'   \item Predicted environmental category (fossil_env_est).
#'   \item Probabilities for each category (fossil_prob_*).
#'   \item Nearest modern site metadata (optional).
#' }
#'
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

  predictions <- list()

  for (i in seq_len(nrow(fossildata))) {
    mb <- fossildata$fossilmbc[i]
    sd <- fossildata$fossilsdc[i]

    if (!is.na(mb) && !is.na(sd)) {
      pred_row <- eco_space %>% dplyr::filter(mbc == mb, sdc == sd)

      if (nrow(pred_row) == 1) {
        predictions[[i]] <- pred_row
      } else {
        predictions[[i]] <- tibble::tibble(env_est = NA)
      }
    } else {
      predictions[[i]] <- tibble::tibble(env_est = NA)
    }
  }

  fossil_predictions <- dplyr::bind_rows(predictions)

  # Combine fossil with prediction results, excluding redundant bin columns
  fossildata <- dplyr::bind_cols(fossildata, fossil_predictions[, -(1:2)])

  # Rename fossil columns for clean naming
  fossildata <- fossildata %>%
    dplyr::rename(
      fossil_mbc = fossilmbc,
      fossil_sdc = fossilsdc,
      fossil_env_est = env_est
    )

  # Rename probability columns too
  prob_cols <- grep("^prob_", names(fossildata), value = TRUE)

  if (length(prob_cols) > 0) {
    names(fossildata)[names(fossildata) %in% prob_cols] <- paste0("fossil_", prob_cols)
  }

  # Optional: match fossils to nearest modern points
  if (match_nearest) {
    if (is.null(fossil_lon) || is.null(fossil_lat) || is.null(modern_id) || is.null(modern_lon) || is.null(modern_lat)) {
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
