test_that("reconstruct_env_qual returns expected structure with nearest match", {
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

  eco_model_qual <- ecometric_model_qual(
    points_df = traitsByPoint$points,
    category_col = "vegetation",
    min_species = 3
  )

  recon_qual <- reconstruct_env_qual(
    fossildata = fossils,
    model_out = eco_model_qual,
    match_nearest = TRUE,
    fossil_lon = "Long",
    fossil_lat = "Lat",
    modern_id = "ID",
    modern_lon = "Longitude",
    modern_lat = "Latitude"
  )

  expect_s3_class(recon_qual, "data.frame")
  expect_true("fossil_env_est" %in% names(recon_qual))
  expect_true("fossil_bin_1" %in% names(recon_qual))
  expect_true("fossil_bin_2" %in% names(recon_qual))
  expect_true(any(grepl("^fossil_prob_", names(recon_qual))))
  expect_true("nearest_modern_point" %in% names(recon_qual))
})
