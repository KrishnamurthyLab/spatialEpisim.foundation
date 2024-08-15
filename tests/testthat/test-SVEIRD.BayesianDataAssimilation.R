test_that("DRC WorldPop data is appropriate", {
  expect_no_error(terra::rast(downloadWorldPopData("COD")))
})

test_that("Population SpatRaster is as described", {
  ## Create this stuff in a fixture that all other tests can use
    subregionsSpatVector <-
      system.file("extdata",
                  "subregionsSpatVector",
                  package = "spatialEpisim.foundation",
                  mustWork = TRUE) %>%
      terra::vect()
    susceptibleSpatRaster <-
      system.file("extdata",
                  "susceptibleSpatRaster.tif",
                  package = "spatialEpisim.foundation",
                  mustWork = TRUE) %>%
      terra::rast()
    layers <- getSVEIRD.SpatRaster(subregionsSpatVector,
                                   susceptibleSpatRaster,
                                   aggregationFactor = 35)
    data("initialInfections.fourCities")

  CongoSpatRaster <- getCountryPopulation.SpatRaster("COD")

  expect_true(is(CongoSpatRaster, "SpatRaster"))
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
    ## Create this stuff in a fixture that all other tests can use
    subregionsSpatVector <-
      system.file("extdata",
                  "subregionsSpatVector",
                  package = "spatialEpisim.foundation",
                  mustWork = TRUE) %>%
      terra::vect()
    susceptibleSpatRaster <-
      system.file("extdata",
                  "susceptibleSpatRaster.tif",
                  package = "spatialEpisim.foundation",
                  mustWork = TRUE) %>%
      terra::rast()
    layers <- getSVEIRD.SpatRaster(subregionsSpatVector,
                                   susceptibleSpatRaster,
                                   aggregationFactor = 35)
    data("initialInfections.fourCities")

  expect_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, -10))

  ## NOTE: it should also work for higher orders, though these will never be
  ## used. Primarily, zero or one will be used.
  expect_no_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 0))
  expect_no_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 1))
  expect_no_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 2))
  expect_no_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 3))
})

test_that("Population remains the same after masking and aggregation", {
    ## Create this stuff in a fixture that all other tests can use
    subregionsSpatVector <-
      system.file("extdata",
                  "subregionsSpatVector",
                  package = "spatialEpisim.foundation",
                  mustWork = TRUE) %>%
      terra::vect()
    susceptibleSpatRaster <-
      system.file("extdata",
                  "susceptibleSpatRaster.tif",
                  package = "spatialEpisim.foundation",
                  mustWork = TRUE) %>%
      terra::rast()
    layers <- getSVEIRD.SpatRaster(subregionsSpatVector,
                                   susceptibleSpatRaster,
                                   aggregationFactor = 35)
    data("initialInfections.fourCities")

  options("geodata_default_path" = "/tmp")
  Congo <- getCountryPopulation.SpatRaster("COD")
  Congo.GADM <- getCountrySubregions.SpatVector("COD", c("Nord-Kivu", "Ituri"))

  expected <- terra::global(Congo, sum, na.rm = TRUE)

  expect_equal(terra::global(maskAndClassifySusceptibleSpatRaster(Congo.GADM, Congo), sum, na.rm = TRUE),
               expected = expected)

  expect_equal(terra::global(getSVEIRD.SpatRaster(Congo.GADM, Congo, 10), "sum", na.rm = TRUE),
               expected = expected)
})
