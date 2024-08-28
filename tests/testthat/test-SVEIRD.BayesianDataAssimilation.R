test_that("DRC WorldPop data is appropriate", {
  expect_no_error(terra::rast(downloadWorldPopData("COD")))
})

test_that("Population SpatRaster is as described", {
  CongoSpatRaster <- getCountryPopulation.SpatRaster("COD")
  expect_true(is(CongoSpatRaster, "SpatRaster"))
  expect_named(CongoSpatRaster, "Population")
})

test_that("Casting seed data in a Queen's neighbourhood of a given order works correctly", {
  layers <- getSVEIRD.SpatRaster(subregionsSpatVector,
                                 susceptibleSpatRaster,
                                 aggregationFactor = 10)
  populationBeforeSeeding <- sum(terra::global(layers, "sum", na.rm = TRUE))

  ## neighbourhoods other than zero and one are prohibited.
  expect_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, -10))
  expect_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 2))

  ## NOTE: it should also work for higher orders, though these will never be
  ## used. Primarily, zero or one will be used.
  expect_no_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 0))
  expect_no_error(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 1))

  ## Seeding should not alter the sum total of the population in the raster.
  castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 0) %>%
    terra::global("sum", na.rm = TRUE) %>%
    sum(na.rm = TRUE) %>%
    expect_equal(expected = populationBeforeSeeding)
})

test_that("Population remains the same after cropping and aggregation", {
  layers <- getSVEIRD.SpatRaster(subregionsSpatVector,
                                 susceptibleSpatRaster,
                                 aggregationFactor = 35)

  options("geodata_default_path" = "/tmp")
  Congo <- getCountryPopulation.SpatRaster("COD")
  expect_equal(sum(terra::global(terra::aggregate(Congo, 10, "sum", na.rm = TRUE), "sum", na.rm = TRUE), na.rm = TRUE),
               expected = sum(terra::global(Congo, "sum", na.rm = TRUE), na.rm = TRUE))
})
