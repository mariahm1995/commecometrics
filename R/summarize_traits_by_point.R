#' Summarize trait distributions at sampling points with optional continent assignment
#'
#' For each spatial sampling point, this function calculates the mean, standard deviation,
#' of a specified trait across and richness all overlapping species polygons.
#' Optionally, it assigns each point to a continent using Natural Earth data.
#'
#' @param points_df A data frame containing sampling points with columns for longitude and latitude.
#' @param trait_df A data frame of trait data. Must include a column for species names ('TaxonName')
#'                and the trait of interest (default = "trait_name").
#' @param species_polygons An `sf` object containing species distribution polygons. Must include a species name column.
#' @param trait_column The name of the trait column in `trait_df` to summarize.
#' @param species_name_col The name of the column in `species_polygons` that contains species names (default = "sci_name").
#' @param continent Logical. If \code{TRUE}, assigns each sampling point to a continent using the Natural Earth shapefile via \code{rnaturalearth::ne_countries()}. If \code{FALSE} (default), no continent assignment is performed.
#' @param lon_col Name of the longitude column in `points_df`. Default is 'Longitude'.
#' @param lat_col Name of the latitude column in `points_df`. Default is 'Latitude'.
#' @param parallel Logical; whether to parallelize the summarization step (default TRUE).
#' @param n_cores Number of cores to use if parallelizing (default: detectCores() - 1).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{points}{A data frame identical to `points_df` but with additional columns:
#'     \describe{
#'       \item{mean_trait}{Mean value of the specified trait across overlapping species.}
#'       \item{sd_trait}{Standard deviation of the trait across overlapping species.}
#'       \item{richness}{Number of species overlapping the point (regardless of trait availability).}
#'       \item{count_trait}{Number of species with non-missing trait values at the point.}
#'       \item{continent}{(Optional) Continent name assigned from Natural Earth data, if `continent = TRUE`.}
#'     }
#'   }
#'   \item{overlap}{A list of character vectors, each containing the names of species whose distribution polygons overlap a given sampling point.}
#' }
#'
#' @examples
#' \dontrun{
#' # Load sample data from the package
#' data("points", package = "commecometrics")
#' data("traits", package = "commecometrics")
#' data("polygons", package = "commecometrics")
#'
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
#' head(traitsByPoint$points)
#' }
#'
#' @export
summarize_traits_by_point <- function(points_df,
                                      trait_df,
                                      species_polygons,
                                      trait_column = "trait_name",
                                      species_name_col = "sci_name",
                                      continent = FALSE,
                                      lon_col = "Longitude",
                                      lat_col = "Latitude",
                                      parallel = TRUE,
                                      n_cores = parallel::detectCores() - 1) {
  if (isTRUE(continent)) {
    continent_shp <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  }

  # Input checks
  if (!inherits(species_polygons, "sf")) {
    stop("'species_polygons' must be an sf object (already converted).")
  }
  if (!all(c(lon_col, lat_col) %in% names(points_df))) {
    stop(paste0("Missing columns in 'points_df': expected '", lon_col, "' and '", lat_col, "'"))
  }
  if (!trait_column %in% names(trait_df)) {
    stop(paste0("Trait column '", trait_column, "' not found in 'trait_df'."))
  }
  if (!species_name_col %in% names(species_polygons)) {
    stop(paste0("Species name column '", species_name_col, "' not found in 'species_polygons'."))
  }

  message("Identifying overlapping species...")

  species_polygons <- sf::st_make_valid(species_polygons)

  points_sf <- sf::st_as_sf(points_df, coords = c(lon_col, lat_col), crs = sf::st_crs(species_polygons))

  if (!sf::st_crs(points_sf) == sf::st_crs(species_polygons)) {
    message("CRS mismatch detected, transforming points to match polygons CRS...")
    points_sf <- sf::st_transform(points_sf, sf::st_crs(species_polygons))
  }

  points_df_coords <- sf::st_coordinates(points_sf) %>% as.data.frame()
  colnames(points_df_coords) <- c(lon_col, lat_col)

  sf::sf_use_s2(FALSE)
  species_overlap_idx <- sf::st_intersects(points_sf, species_polygons)
  sf::sf_use_s2(TRUE)

  species_overlap <- lapply(species_overlap_idx, function(idx) {
    if (length(idx) > 0) {
      species_polygons[[species_name_col]][idx]
    } else {
      NA_character_
    }
  })

  points_df <- sf::st_drop_geometry(points_sf)
  points_df <- dplyr::bind_cols(points_df, points_df_coords)

  trait_lookup <- trait_df %>%
    dplyr::select(TaxonName, {{ trait_column }}) %>%
    tibble::deframe()

  message("Summarizing trait values for each point...")

  summarize_point <- function(overlap_species) {
    if (is.null(overlap_species) || length(overlap_species) == 0 || all(is.na(overlap_species))) {
      return(tibble::tibble(
        mean_trait = NA_real_,
        sd_trait = NA_real_,
        richness = 0,
        count_trait = 0
      ))
    } else {
      trait_values <- trait_lookup[overlap_species]
      tibble::tibble(
        mean_trait = mean(trait_values, na.rm = TRUE),
        sd_trait = sd(trait_values, na.rm = TRUE),
        richness = length(overlap_species),
        count_trait = sum(!is.na(trait_values))
      )
    }
  }

  if (parallel) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, varlist = c("trait_lookup", "summarize_point"), envir = environment())
    trait_summary_list <- parallel::parLapply(cl, species_overlap, summarize_point)
    parallel::stopCluster(cl)
  } else {
    trait_summary_list <- lapply(species_overlap, summarize_point)
  }

  trait_summaries <- dplyr::bind_rows(trait_summary_list)

  points_df <- dplyr::bind_cols(points_df, trait_summaries)

  if (isTRUE(continent)) {
    message("Assigning continent to each point...")

    continent_shp <- sf::st_make_valid(continent_shp)

    points_sf <- sf::st_as_sf(points_df, coords = c(lon_col, lat_col), crs = sf::st_crs(continent_shp))
    joined <- sf::st_join(points_sf, continent_shp, left = TRUE)
    points_df$continent <- joined$continent[match(sf::st_coordinates(points_sf)[, 1], sf::st_coordinates(joined)[, 1])]
    # Correct French Guiana: if labeled as Europe but located in South America by coordinates
    # Typical coordinates for French Guiana are ~longitude -54, latitude ~4 to 6
    points_df$continent[
      points_df$continent == "Europe" &
        points_df[[lon_col]] < -50 & points_df[[lat_col]] < 10
    ] <- "South America"
  }
  return(list(
    points = points_df,
    overlap = species_overlap
  ))
}
