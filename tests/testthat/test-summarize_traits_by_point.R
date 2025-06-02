test_that("summarize_traits_by_point returns expected structure", {
  skip_if_not_installed("sf")
  skip_on_cran()

  # Load internal data
  data("geoPoints", package = "commecometrics")
  data("traits", package = "commecometrics")
  data("spRanges", package = "commecometrics")

  result <- summarize_traits_by_point(
    points_df = geoPoints[1:10, ],
    trait_df = traits,
    species_polygons = spRanges,
    trait_column = "RBL",
    species_name_col = "sci_name",
    continent = FALSE,
    parallel = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("points", "overlap"))
  expect_s3_class(result$points, "data.frame")
  expect_true(all(c("summ_trait_1", "summ_trait_2", "richness", "count_trait") %in% colnames(result$points)))
  expect_type(result$overlap, "list")
  expect_length(result$overlap, 10)
})

test_that("function handles empty overlap correctly", {
  skip_if_not_installed("sf")
  skip_on_cran()

  # Load internal data
  data("geoPoints", package = "commecometrics")
  data("traits", package = "commecometrics")
  data("spRanges", package = "commecometrics")

  dummy_points <- data.frame(Longitude = 0, Latitude = 0)

  result <- summarize_traits_by_point(
    points_df = dummy_points,
    trait_df = traits,
    species_polygons = spRanges,
    trait_column = "RBL",
    species_name_col = "sci_name",
    continent = FALSE,
    parallel = FALSE
  )

  expect_equal(result$points$richness, 0)
  expect_true(is.na(result$points$summ_trait_1))
})
