test_that("ecometric_model runs correctly and returns expected structure", {
  skip_if_not_installed("sf")
  skip_if_not_installed("raster")
  skip_on_cran()

  # Load sample data
  data("geoPoints", package = "commecometrics")
  data("traits", package = "commecometrics")
  data("spRanges", package = "commecometrics")

  # Step 1: summarize traits
  traitsByPoint <- summarize_traits_by_point(
    points_df = geoPoints,
    trait_df = traits,
    species_polygons = spRanges,
    trait_column = "RBL",
    species_name_col = "sci_name",
    continent = FALSE,
    parallel = FALSE
  )

  # Step 2: run ecometric model
  modelResult <- ecometric_model(
    points_df = traitsByPoint$points,
    env_var = "precip",
    transform_fun = function(x) log(x + 1),
    inv_transform_fun = function(x) exp(x) - 1,
    min_species = 3
  )

  # Step 3: check structure
  expect_type(modelResult, "list")
  expect_named(modelResult, c(
    "points_df", "eco_space", "model", "correlation",
    "diagnostics", "settings"
  ), ignore.order = TRUE)

  expect_s3_class(modelResult$model, "lm")
  expect_s3_class(modelResult$correlation, "htest")
  expect_s3_class(modelResult$eco_space, "data.frame")
  expect_s3_class(modelResult$points_df, "data.frame")

  # Check a few specific elements
  expect_true("env_est" %in% colnames(modelResult$points_df))
  expect_true("env_anom" %in% colnames(modelResult$points_df))
  expect_true(modelResult$diagnostics$retained_points > 0)
})
