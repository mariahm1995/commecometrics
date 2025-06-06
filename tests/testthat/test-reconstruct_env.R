test_that("reconstruct_env returns expected structure with nearest match", {
  skip_on_cran()
  skip_if_not_installed("sf")
  skip_if_offline()

  data("geoPoints", package = "commecometrics")
  data("traits", package = "commecometrics")
  data("spRanges", package = "commecometrics")
  data("fossils", package = "commecometrics")

  traits_summary <- summarize_traits_by_point(
    points_df = geoPoints,
    trait_df = traits,
    species_polygons = spRanges,
    trait_column = "RBL",
    species_name_col = "sci_name",
    continent = FALSE,
    parallel = FALSE
  )

  model_out <- ecometric_model(
    points_df = traits_summary$points,
    env_var = "precip",
    transform_fun = function(x) log(x + 1),
    inv_transform_fun = function(x) exp(x) - 1,
    min_species = 3
  )

  recon <- reconstruct_env(
    fossildata = fossils,
    model_out = model_out,
    match_nearest = TRUE,
    fossil_lon = "Long",
    fossil_lat = "Lat",
    modern_id = "ID",
    modern_lon = "Longitude",
    modern_lat = "Latitude"
  )

  expect_s3_class(recon, "data.frame")
  expect_true(all(c("fossil_bin_1", "fossil_bin_2", "fossil_env_est", "fossil_env_est_UN") %in% colnames(recon)))
  expect_true("nearest_modern_point" %in% colnames(recon))
  expect_type(recon$fossil_env_est, "double")
})
