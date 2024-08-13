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

test_that("Casting seed data in a Queen's neighbourhood of a given order works correctly", {
  expect_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, -10))

  ## NOTE: it should also work for higher orders, though these will never be used. Primarily, zero or one will be used.
  expect_no_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 0))
  expect_no_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 1))
  expect_no_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 2))
})

## Create this stuff in a fixture that all other tests can use
subregionsSpatVector <-
  system.file("extdata", "subregionsSpatVector", package = "spatialEpisim.foundation", mustWork = TRUE) %>%
  terra::vect()
susceptibleSpatRaster <-
  system.file("extdata", "susceptibleSpatRaster.tif", package = "spatialEpisim.foundation", mustWork = TRUE) %>%
  terra::rast()
layers <- getSVEIRD.SpatRaster(subregionsSpatVector, susceptibleSpatRaster, aggregationFactor = 35)
data("InitialInfections.fourCities")
