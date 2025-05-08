test_that("ecometric_space_qual returns a valid ggplot list structure", {
  data("points", package = "commecometrics")
  data("traits", package = "commecometrics")
  data("polygons", package = "commecometrics")
  data("fossils", package = "commecometrics")

  traitsByPoint <- summarize_traits_by_point(
    points_df = points,
    trait_df = traits,
    species_polygons = polygons,
    trait_column = "RBL",
    species_name_col = "sci_name",
    continent = FALSE,
    parallel = FALSE
  )

  ecoModelQual <- ecometric_model_qual(
    points_df = traitsByPoint$points,
    category_col = "DOM_NUM",
    min_species = 3
  )

  reconQual <- reconstruct_env_qual(
    fossildata = fossils,
    model_out = ecoModelQual,
    match_nearest = TRUE,
    fossil_lon = "Long",
    fossil_lat = "Lat",
    modern_id = "GlobalID",
    modern_lon = "Longitude",
    modern_lat = "Latitude"
  )

  plotList <- ecometric_space_qual(
    model_out = ecoModelQual,
    fossil_data = reconQual
  )

  # Check top-level output
  expect_type(plotList, "list")
  expect_named(plotList, c("ecometric_space_plot", "probability_maps"))
  expect_s3_class(plotList$ecometric_space_plot, "ggplot")

  # Check only non-null elements of probability_maps
  nonnull_probs <- Filter(Negate(is.null), plotList$probability_maps)
  expect_true(all(vapply(nonnull_probs, inherits, logical(1), what = "ggplot")))
})
