#' Example climate sampling points
#'
#' A subset of 100 global sampling points with associated bioclimatic and vegetation variables.
#' All points overlap with species ranges in the `spRanges` dataset.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{ID}{Unique identifier for each point}
#'   \item{Longitude}{Longitude coordinate (decimal degrees)}
#'   \item{Latitude}{Latitude coordinate (decimal degrees)}
#'   \item{temp}{Mean annual temperature (°C × 10)}
#'   \item{precip}{Annual precipitation (mm)}
#'   \item{vegetation}{Vegetation units (integer code)}
#' }
#' @source Derived from Siciliano-Martina et al. (2024), filtered for overlap with IUCN polygons.
"geoPoints"

#' Relative blade length trait data for Carnivora
#'
#' A dataset of relative blade length (RBL) values for five species in the order Carnivora.
#' These species match those in the `spRanges` dataset.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{TaxonName}{Species name (binomial)}
#'   \item{RBL}{Relative blade length (unitless ratio)}
#' }
#' @source Siciliano-Martina et al. (2024). *Ecology and Evolution*, 14(10), e70214.
"traits"

#' Species distribution polygons for five Carnivora species
#'
#' A spatial dataset of species range polygons matching the species in the `traits` dataset.
#'
#' @format An `sf` object with the following columns:
#' \describe{
#'   \item{TaxonName}{Species name (matching the `traits` table)}
#'   \item{geometry}{Polygon geometry representing species distribution}
#' }
#' @source Download from the IUCN Red List webpage (IUCN, 2025).
"spRanges"

#' Fossil trait data for projection onto ecometric space
#'
#' A dataset of fossil sites with estimated trait distribution and geographic coordinates,
#' used to project extinct communities onto modern ecometric space.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Site}{Unique identifier for the fossil sample}
#'   \item{Mean}{Estimated mean of relative blade length for the fossil site}
#'   \item{SD}{Estimated sd of relative blade length for the fossil site}
#'   \item{Long}{Longitude coordinate (decimal degrees)}
#'   \item{Lat}{Latitude coordinate (decimal degrees)}
#' }
#' @source Siciliano-Martina et al. (2024). *Ecology and Evolution*, 14(10), e70214.
"fossils"
