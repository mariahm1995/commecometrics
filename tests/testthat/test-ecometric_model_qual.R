test_that("ecometric_model_qual runs correctly and returns expected structure", {
  skip_on_cran()
  skip_if_not_installed("sf")

  # Load sample data from package
  data("geoPoints", package = "commecometrics")
  data("traits", package = "commecometrics")
  data("spRanges", package = "commecometrics")

  # Summarize traits by point
  traitsByPoint <- summarize_traits_by_point(
    points_df = geoPoints,
    trait_df = traits,
    species_polygons = spRanges,
    trait_column = "RBL",
    species_name_col = "sci_name",
    continent = FALSE,
    parallel = FALSE
  )

  # Run ecometric model for qualitative variable
  model_out <- ecometric_model_qual(
    points_df = traitsByPoint$points,
    category_col = "vegetation",
    min_species = 3
  )

  # Check output structure
  expect_type(model_out, "list")
  expect_named(model_out, c("points_df", "eco_space", "diagnostics", "settings", "prediction_accuracy"))

  # Check that points_df contains expected new columns
  expect_s3_class(model_out$points_df, "data.frame")
  expect_true(all(c("bin_1", "bin_2", "observed_probability",
                    "predicted_probability", "predicted_category",
                    "correct_prediction", "env_anom") %in% colnames(model_out$points_df)))

  # Check that eco_space is a data frame with mode predictions
  expect_s3_class(model_out$eco_space, "data.frame")
  expect_true(all(c("bin_1", "bin_2", "env_est") %in% colnames(model_out$eco_space)))

  # Check prediction accuracy is numeric between 0 and 100
  expect_type(model_out$prediction_accuracy, "double")
  expect_gte(model_out$prediction_accuracy, 0)
  expect_lte(model_out$prediction_accuracy, 100)
})
