test_that("ecometric_space returns a ggplot object", {
  skip_on_cran()

  data("geoPoints", package = "commecometrics")
  data("traits", package = "commecometrics")
  data("spRanges", package = "commecometrics")
  data("fossils", package = "commecometrics")

  traitsByPoint <- summarize_traits_by_point(
    points_df = geoPoints,
    trait_df = traits,
    species_polygons = spRanges,
    trait_column = "RBL",
    species_name_col = "sci_name",
    continent = FALSE,
    parallel = FALSE
  )

  ecoModel <- ecometric_model(
    points_df = traitsByPoint$points,
    env_var = "precip",
    transform_fun = function(x) log(x + 1),
    inv_transform_fun = function(x) exp(x) - 1,
    min_species = 3
  )

  recon <- reconstruct_env(
    fossildata = fossils,
    model_out = ecoModel,
    match_nearest = TRUE,
    fossil_lon = "Long",
    fossil_lat = "Lat",
    modern_id = "ID",
    modern_lon = "Longitude",
    modern_lat = "Latitude"
  )

  plot <- ecometric_space(
    model_out = ecoModel,
    env_name = "Precipitation (log mm)",
    fossil_data = recon
  )

  expect_s3_class(plot, "gg")
})
