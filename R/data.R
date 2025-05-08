#' Example climate sampling points
#'
#' A subset of 500 global sampling points with associated bioclimatic and land cover variables.
#' All points overlap with species ranges in the `polygons` dataset.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{GlobalID}{Unique identifier for each point}
#'   \item{Longitude}{Longitude coordinate (decimal degrees)}
#'   \item{Latitude}{Latitude coordinate (decimal degrees)}
#'   \item{BIO1}{Mean annual temperature (°C × 10)}
#'   \item{BIO12}{Annual precipitation (mm)}
#'   \item{DOM_NUM}{Dominant land cover class (integer code)}
#' }
#' @source Derived from Siciliano-Martina et al. (2021), filtered for overlap with IUCN polygons.
"points"

#' Relative blade length trait data for Carnivora
#'
#' A dataset of relative blade length (RBL) values for 10 species in the order Carnivora.
#' These species match those in the `polygons` dataset.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{TaxonName}{Species name (binomial)}
#'   \item{RBL}{Relative blade length (unitless ratio)}
#' }
#' @source Siciliano-Martina et al. (2024). *Ecology and Evolution*, 14(10), e70214.
"traits"

#' Species distribution polygons for 4 Carnivora species
#'
#' A spatial dataset of species range polygons matching the species in the `traits` dataset.
#'
#' @format An `sf` object with the following columns:
#' \describe{
#'   \item{TaxonName}{Species name (matching the `traits` table)}
#'   \item{geometry}{Polygon geometry representing species distribution}
#' }
#' @source Donwload from the IUCN Red List webpage (IUCN, 2025).
"polygons"

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
