test_that("DRC WorldPop data is appropriate", {
  rast(downloadWorldPopData("COD"))
  expect_equal(2 * 2, 4)
})

test_that("Population SpatRaster is as described", {
  ## TODO: integrate this into a larger fixture...
  CongoSpatRaster <- getCountryPopulation.SpatRaster("COD")

  expect_type(CongoSpatRaster, "SpatRaster")
  expect_false(any(is.na(CongoSpatRaster)))
  expect_named(CongoSpatRaster, "Susceptible")
})

test_that("Average Euclidean distance is correct", {
})

test_that("Transmission likelihood weightings are corect", {
})

test_that("The Linear Interpolation Operator, H, is correct", {
})

"Forecast error covariance matrix, Q, is correct"


