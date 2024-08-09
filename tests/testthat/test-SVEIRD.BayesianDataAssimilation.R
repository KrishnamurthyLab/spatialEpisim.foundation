test_that("DRC WorldPop data is appropriate", {
  expect_no_error(terra::rast(downloadWorldPopData("COD")))
})

test_that("Population SpatRaster is as described", {
  ## TODO: integrate the following data creation into a larger test fixture...
  CongoSpatRaster <- getCountryPopulation.SpatRaster("COD")

  expect_type(CongoSpatRaster, "SpatRaster")
  expect_false(any(is.na(CongoSpatRaster)))
  expect_named(CongoSpatRaster, "Susceptible")
})

## test_that("Average Euclidean distance is correct", {
## })

## test_that("Transmission likelihood weightings are corect", {
## })

## test_that("The Linear Interpolation Operator, H, is correct", {
## })

## "Forecast error covariance matrix, Q, is correct"

## "SVEIRD compartmental epidemic model with Bayesian data assimilation is correct"

## "Setup of Bayesian data assimilation is correct"

## "Bayesian data assimilation works correctly"

## "Casting seed data in a Moore neighbourhood works correctly"
