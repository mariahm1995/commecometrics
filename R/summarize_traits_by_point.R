#' Summarize Trait Distributions at Sampling Points with Optional Continent Assignment
#'
#' For each spatial sampling point, this function calculates the mean, standard deviation,
#' and richness of a specified trait across all overlapping species polygons.
#' Optionally, it assigns each point to a continent based on a shapefile.
#'
#' @param points_df A data frame containing sampling points with columns for longitude and latitude.
#' @param trait_df A data frame of trait data. Must include a column for species names (typically 'TaxonName')
#'                and the trait of interest (e.g., 'RBL').
#' @param species_polygons An `sf` object containing species distribution polygons. Must include a species name column.
#' @param trait_column The name of the trait column in `trait_df` to summarize.
#' @param species_name_col The name of the column in `species_polygons` that contains species names (default = "binomial").
#' @param continent_shp Optional. A shapefile (`sf` object) with continental boundaries and a column named 'CONTINENT'.
#' @param lon_col Name of the longitude column in `points_df`. Default is 'Longitude'.
#' @param lat_col Name of the latitude column in `points_df`. Default is 'Latitude'.
#' @param parallel Logical; whether to parallelize the summarization step (default TRUE).
#' @param n_cores Number of cores to use if parallelizing (default: detectCores() - 1).
#'
#' @return A list with:
#'   \item{points}{The input `points_df`, updated with trait summary columns and optionally continent column.}
#'   \item{overlap}{A list of species names overlapping each point.}
#'
#' @export
summarize_traits_by_point <- function(points_df,
                                      trait_df,
                                      species_polygons,
                                      trait_column = NULL,
                                      species_name_col = "binomial",
                                      continent_shp = FALSE,
                                      lon_col = "Longitude",
                                      lat_col = "Latitude",
                                      parallel = TRUE,
                                      n_cores = parallel::detectCores() - 1) {
  if(continent_shp){
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

  if (continent_shp) {
    message("Assigning continent to each point...")

    if (inherits(continent_shp, "Spatial")) {
      continent_shp <- sf::st_as_sf(continent_shp)
    }
    continent_shp <- sf::st_make_valid(continent_shp)

    points_sf <- sf::st_as_sf(points_df, coords = c(lon_col, lat_col), crs = sf::st_crs(continent_shp))
    joined <- sf::st_join(points_sf, continent_shp, left = TRUE)

    if ("continent" %in% names(joined)) {
      points_df$continent <- joined$continent
    } else {
      warning("The 'continent' column was not found in the shapefile.")
      points_df$continent <- NA
    }
  }
  return(list(
    points = points_df,
    overlap = species_overlap
  ))
}
