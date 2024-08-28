## NOTE: used in both the originalCodeResults and newCodeResults
rasterAggregationFactor <- 20

originalCodeResults <- {
  generateLIO2 <- function(rasterStack, sitRepData, states_observable = 2) {
    nrows <- nrow(rasterStack)
    ncols <- ncol(rasterStack)
    p <- ncell(rasterStack)

    Locations <- read.csv(file = sitRepData, header = T)
    nHealthZones <- dim(Locations)[1]

    cellFromXY(rasterStack, cbind(27.13, 3.72)) # 1. Note: This is the top left corner cell
    cellFromXY(rasterStack, cbind(29.47306, 0.49113)) # 1929. Note: This is the Lon, Lat for Beni
    cellFromXY(rasterStack, cbind(0.49113, 29.47306)) # NA will be produced if you flip the (Lon, Lat) to (Lat, Lon)
    cellFromXY(rasterStack, cbind(31.29, -2.19)) # 3550. Note This is the bottom righ corner cell

    Hpos <- cellFromXY(rasterStack, as.matrix(Locations[, 3:2]))
    message("Hpos is")
    message(Hpos)

    ########## FOR DIAGNOSTICS##########################
    rows <- rowFromY(rasterStack, as.vector(Locations[, 2]))
    rows
    cols <- colFromX(rasterStack, as.vector(Locations[, 3]))
    cols

    # message('A test:')
    # message(Hpos[5])
    # message(rows[5])
    # message(cols[5])

    # Hpos for Beni is obtained as = (rows - 1)*ncols + cols = (39-1)*50 + 29 = 1929
    ###################################################

    if (!(anyDuplicated(Hpos) == 0)) {
      message("Warning: Duplicate Indices")
    }

    #------------------------------------------------------------------------#
    # Compute H matrix, the linear interpolation operator of dimension q x p #
    #------------------------------------------------------------------------#

    H <- H0 <- matrix(0, nHealthZones, p)

    for (i in 1:nHealthZones) {
      # message(paste("Hposition:", Hpos[i]))
      H[i, Hpos[i] - 1] <- 0.08
      H[i, Hpos[i] + 2] <- 0.04
      H[i, Hpos[i] - 2] <- 0.04
      H[i, Hpos[i] + 2 * ncols] <- 0.04
      H[i, Hpos[i] + 2 * ncols + 1] <- 0.04
      H[i, Hpos[i] + 2 * ncols - 1] <- 0.04
      H[i, Hpos[i] - 2 * ncols] <- 0.04
      H[i, Hpos[i] + 2 + ncols] <- 0.04
      H[i, Hpos[i] + 2 + 2 * ncols] <- 0.04
      H[i, Hpos[i] + 2 - 2 * ncols] <- 0.04
      H[i, Hpos[i] + 2 - ncols] <- 0.04
      H[i, Hpos[i] - 2 + 2 * ncols] <- 0.04
      H[i, Hpos[i] - 1 + 2 * ncols] <- 0.04
      H[i, Hpos[i] + 1 + 2 * ncols] <- 0.04
      H[i, Hpos[i] - 2 - ncols] <- 0.04
      H[i, Hpos[i] - 2 - 2 * ncols] <- 0.04
      H[i, Hpos[i] - 2 + ncols] <- 0.04
      H[i, Hpos[i] + 1] <- 0.08
      H[i, Hpos[i] - ncols] <- 0.08
      H[i, Hpos[i] + ncols] <- 0.08
      H[i, Hpos[i] - ncols - 1] <- 0.08
      H[i, Hpos[i] - ncols + 1] <- 0.08
      H[i, Hpos[i] + ncols - 1] <- 0.08
      H[i, Hpos[i] + ncols + 1] <- 0.08
      H[i, Hpos[i]] <- 0.12
    }
    H <- H * 5 / 7

    message(paste("Number of Health Zones:", nHealthZones))

    Hmat <- H

    # if (states_observable == 2)
    # {
    #   Htop <- cbind(H, H0)
    #   Hbottom <- cbind(H0, H)
    #
    #   Hmat <- rbind(Htop, Hbottom)
    # } else {
    #   Hmat <- H
    # }

    if (states_observable > 1) {
      for (n in seq(from = 1, to = states_observable - 1, by = 1)) {
        Hmat <- cbind(Hmat, H0)
      }

      for (n in seq(from = 1, to = states_observable - 1, by = 1)) {
        Htop <- matrix(0, nHealthZones, n * p)
        if ((n + 1 - states_observable) != 0) {
          Hbottom <- matrix(0, nHealthZones, (states_observable - n - 1) * p)
          # message(dim(Htop))
          # message(dim(H))
          # message(dim(Hbottom))
          Hmat <- rbind(Hmat, cbind(Htop, H, Hbottom))
        } else {
          # message(dim(H))
          # message(dim(Htop))
          Hmat <- rbind(Hmat, cbind(Htop, H))
        }
      }
    }

    message(paste("Dimension of the linear interpolation operator, H:", dim(Hmat)))
    message(paste("Row sums of H matrix:", Hmat))

    # message(table(colSums(Hmat)))
    # message(sum(Hmat))
    # message(table(Hmat))

    return(list(
      "Hmat" = Hmat,
      "Locations" = Locations,
      "rasterStack" = rasterStack,
      "states_observable" = states_observable
    ))
  }

  COD <- "Democratic Republic of Congo"
  ## TODO: change this to the 
  sitRepData <- test_path("observeddata", "Ebola_Health_Zones_LatLon.csv")
  rs <- createRasterStack(
    selectedCountry = COD,
    rasterAgg = rasterAggregationFactor,
    isCropped = TRUE,
    level1Names = c("Ituri", "Nord-Kivu"),
    susceptibleLayer = createSusceptibleLayer(COD, rasterAgg = rasterAggregationFactor)
  )$rasterStack

  list(
    oneStateObserved = suppressMessages(generateLIO2(
      rs,
      sitRepData,
      states_observable = 1
    )$Hmat),
    twoStateObserved = suppressMessages(generateLIO2(
      rs,
      sitRepData,
      states_observable = 2
    )$Hmat)
  )
}

test_that("Linear interpolation operator (LIO2) is created without issue, acknowledging expected warnings about over-aggregation.", {
  ## NOTE: when the rasterAggregationFactor changes the expected warnings will
  ## dissappear and so that will need adjustment if this is changed, so it is
  ## tested!
  expect_equal(rasterAggregationFactor, 20)

  ## NOTE: other tests are done within linearInterpolationOperator() which check
  ## that the row sums are appropriate and that the sum of all row sums is also
  ## appropriate, regardless of the ultimate size of the operator matrix. These
  ## are expected not to fail, so we can expect no error (stop() would raise
  ## errors).
  expect_no_error({
    sveirdLayers <<- getSVEIRD.SpatRaster(subregionsSpatVector,
                                          susceptibleSpatRaster,
                                          aggregationFactor = rasterAggregationFactor)
    newCodeResults <<- list(
      oneStateObserved =
        expect_warning(linearInterpolationOperator(
          layers = sveirdLayers,
          healthZoneCoordinates = healthZonesCongo,
          neighbourhood.order = 0,
          compartmentsReported = 1
        )),
      twoStateObserved =
        expect_warning(linearInterpolationOperator(
          layers = sveirdLayers,
          healthZoneCoordinates = healthZonesCongo,
          neighbourhood.order = 0,
          compartmentsReported = 2
        ))
    )
  })

  ## DONE: the intention is that the geographical attributes of the SpatRasters
  ## should be equal, as a minimum prerequisite for constructing the linear
  ## interpolation operators.
  expect_true(identical(crs(rs), crs(sveirdLayers)))

  ## TODO: these should be removed, because they aren't MEANINGFUL other than an
  ## expression of what is already known, that they differ (but the intention is
  ## not that they SHOULD).
  ## expect_false(identical(nrow(rs), nrow(sveirdLayers)))
  ## expect_false(identical(ncol(rs), ncol(sveirdLayers)))
  ## expect_false(identical(res(rs), res(sveirdLayers)))
  ## expect_false(identical(ext(rs), ext(sveirdLayers)))
  ## expect_false(identical(dim(originalCodeResults$oneStateObserved),
  ##                        dim(newCodeResults$oneStateObserved)))
  ## expect_false(identical(dim(originalCodeResults$twoStateObserved),
  ##                        dim(newCodeResults$twoStateObserved)))

})
