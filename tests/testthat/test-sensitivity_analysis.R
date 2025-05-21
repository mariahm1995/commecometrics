test_that("sensitivity_analysis returns valid structure", {
  skip_on_cran()

  data("geoPoints", package = "commecometrics")
  data("traits", package = "commecometrics")
  data("spRanges", package = "commecometrics")

  traitsByPoint <- summarize_traits_by_point(
    points_df = geoPoints,
    trait_df = traits,
    species_polygons = spRanges,
    trait_column = "RBL",
    species_name_col = "sci_name",
    continent = FALSE,
    parallel = FALSE
  )

  # Run a minimal version for testing
  result <- sensitivity_analysis(
    points_df = traitsByPoint$points,
    env_var = "precip",
    sample_sizes = c(20, 30),
    iterations = 2,
    transform_fun = function(x) log(x + 1),
    parallel = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("combined_results", "summary_results"))

  expect_s3_class(result$combined_results, "data.frame")
  expect_s3_class(result$summary_results, "data.frame")

  expect_true(all(c("SampleSize", "Iteration", "Training_Mean_Anomaly",
                    "Training_Correlation", "Testing_Mean_Anomaly", "Testing_Correlation") %in%
                    names(result$combined_results)))
})
