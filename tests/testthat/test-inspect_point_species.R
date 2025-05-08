test_that("inspect_point_species returns a leaflet map", {
  skip_if_not_installed("sf")
  skip_if_not_installed("leaflet")
  skip_on_cran()

  # Load internal sample data
  data("points", package = "commecometrics")
  data("traits", package = "commecometrics")
  data("polygons", package = "commecometrics")

  # Step 1: Run the summarization function
  traitsByPoint <- summarize_traits_by_point(
    points_df = points,
    trait_df = traits,
    species_polygons = polygons,
    trait_column = "RBL",
    species_name_col = "sci_name",
    continent = FALSE,
    parallel = FALSE
  )

  # Step 2: Generate map with random points
  map_obj <- inspect_point_species(
    traits_summary = traitsByPoint,
    n_random = 5,
    min_species_valid = 2
  )

  expect_s3_class(map_obj, "leaflet")

  # Step 3: Generate map using specific point IDs
  specific_ids <- traitsByPoint$points$GlobalID[1:3]
  map_obj2 <- inspect_point_species(
    traits_summary = traitsByPoint,
    point_ids = specific_ids,
    min_species_valid = 1
  )

  expect_s3_class(map_obj2, "leaflet")
})
