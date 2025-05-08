test_that("sensitivity_analysis_qual returns expected structure", {
  skip_on_cran()
  skip_if_not_installed("commecometrics")

  data("points", package = "commecometrics")
  data("traits", package = "commecometrics")
  data("polygons", package = "commecometrics")

  traitsByPoint <- summarize_traits_by_point(
    points_df = points,
    trait_df = traits,
    species_polygons = polygons,
    trait_column = "RBL",
    species_name_col = "sci_name",
    continent = FALSE,
    parallel = FALSE
  )

  out <- sensitivity_analysis_qual(
    points_df = traitsByPoint$points,
    category_col = "DOM_NUM",
    sample_sizes = c(20, 30),
    iterations = 2,
    parallel = FALSE
  )

  expect_type(out, "list")
  expect_named(out, c("combined_results", "summary_results"))

  expect_s3_class(out$combined_results, "data.frame")
  expect_s3_class(out$summary_results, "data.frame")

  expect_true(all(c("SampleSize", "Iteration", "Training_Accuracy", "Testing_Accuracy") %in% colnames(out$combined_results)))
  expect_true(all(c("SampleSize", "Training_Accuracy", "Testing_Accuracy") %in% colnames(out$summary_results)))

  expect_true(all(out$combined_results$Training_Accuracy >= 0 & out$combined_results$Training_Accuracy <= 1))
  expect_true(all(out$combined_results$Testing_Accuracy >= 0 & out$combined_results$Testing_Accuracy <= 1))
})
