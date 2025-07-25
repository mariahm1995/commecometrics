#' Inspect overlapping species at sampling points
#'
#' Creates an interactive map to verify species overlap at selected points.
#'
#' @param traits_summary A list output from `summarize_traits_by_point()`, containing summarized trait values (`$points`) and species overlaps (`$overlap`).
#' @param point_ids Optional. A vector of specific point IDs to inspect. If NULL, selects \code{n_random} points at random.
#' @param n_random Number of random points to inspect if \code{point_ids} not provided (default = 10).
#' @param lon_col Name of the longitude column in \code{points} (default = "Longitude").
#' @param lat_col Name of the latitude column in \code{points} (default = "Latitude").
#' @param ID_col Name of the ID column in \code{points} (default = "ID").
#' @param min_species_valid Minimum number of species with trait data to consider a point valid (default = 3).
#' @param env_var Optional. Name of the environmental variable column in \code{points} to include in popup.
#'
#' @return An interactive leaflet map showing selected points with species list popups.
#'
#' @import leaflet
#'
#' @examples
#' \donttest{
#' # Load sample data from the package
#' data("geoPoints", package = "commecometrics")
#' data("traits", package = "commecometrics")
#' data("spRanges", package = "commecometrics")
#'
#' # Summarize traits at points
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
#' # Visualize a random sample of 10 points
#' inspect_point_species(
#'   traits_summary = traitsByPoint,
#'   n_random = 10,
#'   min_species_valid = 3
#' )
#' }
#' @export
#'
inspect_point_species <- function(traits_summary,
                                  point_ids = NULL,
                                  n_random = 10,
                                  lon_col = "Longitude",
                                  lat_col = "Latitude",
                                  ID_col = "ID",
                                  min_species_valid = 3,
                                  env_var = NULL) {
  points_df <- traits_summary$points
  species_overlap <- traits_summary$overlap

  required_cols <- c("summ_trait_1", "summ_trait_2", "richness", "count_trait")
  missing_cols <- setdiff(required_cols, names(points_df))
  if (length(missing_cols) > 0) {
    stop("Missing expected columns in traits_summary$points: ", paste(missing_cols, collapse = ", "))
  }

  # Determine which points to inspect
  if (is.null(point_ids)) {
    selected_idx <- sample(seq_len(nrow(points_df)), n_random)
  } else {
    selected_idx <- which(points_df[[ID_col]] %in% point_ids)
    if (length(selected_idx) == 0) stop("Provided point IDs not found.")
  }

  selected_points <- points_df[selected_idx, ]
  selected_species <- species_overlap[selected_idx]

  # Build popup text
  popup_texts <- purrr::map2_chr(
    seq_len(nrow(selected_points)),
    selected_species,
    function(i, sp) {
      pt <- selected_points[i, ]
      species_list <- if (length(sp) > 0 && !all(is.na(sp))) {
        paste(sp, collapse = "<br>")
      } else {
        "No species recorded"
      }

      paste0(
        "<b>Point ID:</b> ", pt[[ID_col]], "<br><br>",
        "<b>Summary metric 1:</b> ", round(pt$summ_trait_1, 2), "<br>",
        "<b>Summary metric 2:</b> ", round(pt$summ_trait_2, 2), "<br>",
        "<b>Richness:</b> ", pt$richness, "<br>",
        "<b>Richness (species with trait):</b> ", pt$count_trait, "<br>",
        if (!is.null(env_var)) paste0("<b>Environmental variable value:</b> ", round(pt[[env_var]], 2), "<br>") else "",
        "<br><b>Species list:</b><br>", species_list
      )
    }
  )

  # Assign color based on user-defined threshold
  point_colors <- ifelse(selected_points$count_trait >= min_species_valid, "#3182bd", "#de2d26")

  # Precompute coordinates
  lon_vals <- selected_points[[lon_col]]
  lat_vals <- selected_points[[lat_col]]

  # Create leaflet map
  m <- leaflet::leaflet(data = selected_points) %>%
    leaflet::addTiles() %>%
    leaflet::addCircleMarkers(
      lng = lon_vals,
      lat = lat_vals,
      popup = popup_texts,
      radius = 5,
      color = point_colors,
      stroke = TRUE,
      fillOpacity = 0.8
    ) %>%
    leaflet::addLegend(
      "bottomright",
      colors = c("#3182bd", "#de2d26"),
      labels = c(
        paste0("Valid (equal or more than ", min_species_valid, " species)"),
        paste0("Not valid (less than ", min_species_valid, " species)")
      ),
      title = "Point Status",
      opacity = 1
    )

  return(m)
}
