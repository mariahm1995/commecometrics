#' Reconstruct past environmental conditions using ecometric models
#'
#' Uses fossil community trait summaries (mean and SD) to reconstruct
#' past environmental conditions by projecting them onto a binned ecometric trait space
#' built from modern data. Optionally, it also assigns each fossil point to the nearest
#' modern sampling site to retrieve observed environmental data.
#'
#' @param fossildata A data frame containing fossil trait summaries per fossil site.
#'                   Must include columns for `Mean` and `SD` of the trait.
#' @param model_out Output list from \code{run_ecometric_model()}, containing modern data, diagnostics, and model settings.
#' @param inv_transform A function to back-transform environmental estimates to the original scale.
#'                      Default is \code{exp(x) - 1}. If \code{NULL}, the inverse transform stored in \code{model_out} is used if available.
#' @param ci The width of the interval to calculate around the maximum likelihood estimate (default = 0.05).
#' @param match_nearest Logical; if TRUE, the function matches each fossil to its nearest modern point based on coordinates (default = TRUE).
#' @param fossil_lon Name of the longitude column in `fossildata`. Required if \code{match_nearest = TRUE}.
#' @param fossil_lat Name of the latitude column in `fossildata`. Required if \code{match_nearest = TRUE}.
#' @param modern_id Name of the unique ID column in modern points (e.g., "GlobalID").
#' @param modern_lon Name of the longitude column in modern points. Required if \code{match_nearest = TRUE}.
#' @param modern_lat Name of the latitude column in modern points. Required if \code{match_nearest = TRUE}.
#' @param crs_proj Coordinate reference system to use when converting fossil and modern data to sf format (default = EPSG:4326).
#'
#' @return A data frame (`fossildata`) with reconstructed environmental values and optional nearest modern point data. Includes the following additional columns:
#' \describe{
#'   \item{fossil_mbc}{Assigned bin number for mean trait (based on fossil trait mean).}
#'   \item{fossil_sdc}{Assigned bin number for SD trait (based on fossil trait SD).}
#'   \item{fossil_env_est}{Maximum likelihood estimate of the environmental variable (on transformed scale if applicable).}
#'   \item{fossil_minlimit}{Lower bound of the confidence interval around the environmental estimate (transformed scale).}
#'   \item{fossil_maxlimit}{Upper bound of the confidence interval around the environmental estimate (transformed scale).}
#'   \item{fossil_env_est_UN}{(Optional) Inverse-transformed environmental estimate, on the original scale.}
#'   \item{fossil_minlimit_UN}{(Optional) Inverse-transformed lower bound of the confidence interval.}
#'   \item{fossil_maxlimit_UN}{(Optional) Inverse-transformed upper bound of the confidence interval.}
#'   \item{nearest_modern_point}{(Optional) ID of the nearest modern sampling point (if \code{match_nearest = TRUE}).}
#'   \item{...}{Additional columns from the matched modern site if \code{match_nearest = TRUE} (e.g., observed environmental values).}
#' }
#' @examples
#' \dontrun{
#' # Load internal data
#' data("points", package = "commecometrics")
#' data("traits", package = "commecometrics")
#' data("polygons", package = "commecometrics")
#' data("fossils", package = "commecometrics")
#'
#' # Step 1: Summarize modern trait values at sampling points
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
#' # Step 2: Run an ecometric model with BIO12 (precipitation)
#' ecoModel <- ecometric_model(
#'   points_df = traitsByPoint$points,
#'   env_var = "BIO12",
#'   transform_fun = function(x) log(x + 1),
#'   inv_transform_fun = function(x) exp(x) - 1,
#'   min_species = 3
#' )
#'
#' # Step 3: Reconstruct fossil environments
#' recon <- reconstruct_env(
#'   fossildata = fossils,
#'   model_out = ecoModel,
#'   match_nearest = TRUE,
#'   fossil_lon = "Long",
#'   fossil_lat = "Lat",
#'   modern_id = "GlobalID",
#'   modern_lon = "Longitude",
#'   modern_lat = "Latitude"
#' )
#'
#' }
#' @export
reconstruct_env <- function(fossildata,
                            model_out,
                            inv_transform = NULL,
                            ci = 0.05,
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
  modern_points <- model_out$points_df

  if (model_out$settings$transformed) {
    env_col <- "env_trans"
  } else {
    env_col <- model_out$settings$env_var
  }

  mbc_col <- "mbc"
  sdc_col <- "sdc"

  # Use model's inverse transform if inv_transform not provided
  if (missing(inv_transform) || is.null(inv_transform)) {
    if (!is.null(model_out$settings$inv_transform_fun)) {
      inv_transform <- model_out$settings$inv_transform_fun
      message("Using inverse transformation function stored in model settings.")
    } else {
      message("No inverse transformation function found. Estimates will remain in original or transformed space.")
      inv_transform <- NULL
    }
  }

  # Assign bins to fossil points
  fossildata$fossilmbc <- .bincode(fossildata$Mean, breaks = mbrks)
  fossildata$fossilsdbc <- .bincode(fossildata$SD, breaks = sdbrks)

  fossilmodmax <- list()

  for (i in seq_len(nrow(fossildata))) {
    mb <- fossildata$fossilmbc[i]
    sd <- fossildata$fossilsdbc[i]

    if (!is.na(mb) && !is.na(sd)) {
      idx <- which(modern_points[[mbc_col]] == mb & modern_points[[sdc_col]] == sd)
      if (length(idx) > 0) {
        dat <- modern_points[[env_col]][idx]
        dens <- density(dat[!is.na(dat)], bw = 1)
        mode_val <- dens$x[which.max(dens$y)]

        n <- length(dens$x)
        mode_idx <- which.max(dens$y)
        bound_n <- floor(n * ci)

        lower_idx <- max(1, mode_idx - bound_n)
        upper_idx <- min(n, mode_idx + bound_n)

        fossilmodmax$envest[i] <- mode_val
        fossilmodmax$minlimit[i] <- dens$x[lower_idx]
        fossilmodmax$maxlimit[i] <- dens$x[upper_idx]
      } else {
        fossilmodmax$envest[i] <- NA_real_
        fossilmodmax$minlimit[i] <- NA_real_
        fossilmodmax$maxlimit[i] <- NA_real_
      }
    } else {
      fossilmodmax$envest[i] <- NA_real_
      fossilmodmax$minlimit[i] <- NA_real_
      fossilmodmax$maxlimit[i] <- NA_real_
    }
  }

  fossildata$envest <- fossilmodmax$envest
  fossildata$minlimit <- fossilmodmax$minlimit
  fossildata$maxlimit <- fossilmodmax$maxlimit

  if (!is.null(inv_transform)) {
    fossildata$envestUN <- inv_transform(fossildata$envest)
    fossildata$minlimitUN <- inv_transform(fossildata$minlimit)
    fossildata$maxlimitUN <- inv_transform(fossildata$maxlimit)
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

    fossildata$near.point <- modern_points[[modern_id]][nearest_idx]

    fossildata <- merge(fossildata, modern_points, by.x = "near.point", by.y = modern_id)
  }

  # Final step: Rename fossil-specific columns for clarity
  if (!is.null(inv_transform)) {
    fossildata <- fossildata %>%
      dplyr::rename(
        fossil_mbc = fossilmbc,
        fossil_sdc = fossilsdbc,
        fossil_env_est = envest,
        fossil_minlimit = minlimit,
        fossil_maxlimit = maxlimit,
        fossil_env_est_UN = envestUN,
        fossil_minlimit_UN = minlimitUN,
        fossil_maxlimit_UN = maxlimitUN,
        nearest_modern_point = near.point
      )
  } else {
    fossildata <- fossildata %>%
      dplyr::rename(
        fossil_mbc = fossilmbc,
        fossil_sdc = fossilsdbc,
        fossil_env_est = envest,
        fossil_minlimit = minlimit,
        fossil_maxlimit = maxlimit,
        nearest_modern_point = near.point
      )
  }

  return(fossildata)
}
