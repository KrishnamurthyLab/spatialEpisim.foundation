############################################################################
## TODO: the states_observable/compartmentsReported needs to be tested as ##
## well; with one state oberved (one compartment reported) the resulting  ##
## variance-covariance (Q) matrix is not a block matrix, while with two   ##
## compartments reported it is a block matrix.                            ##
############################################################################

## These variables are treated as constants to be used throughout the testing
## conducted in this file.
COD <- "Democratic Republic of Congo"
QVar <- 2
QCorrLength <- 0.8
max.print.old <- getOption("max.print")
statesObservable <- 2
rasterAggregationFactor <- 35
variableCovarianceFunctions <- c("DBD", "Balgovind", "Exponential", "Gaussian", "Spherical")

## This compound expression (defining originalCodeResults) contains the
## necessary definitions to run the original Q matrix function with all of its
## unique arguments. Only neighborhood (nbhd) and the variance-covariance
## function's different values are changed. QVar, QCorrLength,
## states_observable, and the input raster are all kept the same between the
## original code and the new code to allow a simpler but still useful
## comparison.
originalCodeResults <- {
  suppressPackageStartupMessages({
    library(Matrix)
    library(raster)
    library(countrycode)
    library(terra)
  })

  genQ <- function(nrows, ncols, varCovarFunc, QVar, QCorrLength, nbhd, states_observable = 2) {
    p <- nrows * ncols
    message(p)

    Q <- matrix(0, p, p)

    # uninhabitableCells <- c()
    #
    # for (a in 1:nrows) {
    #   for (b in 1:ncols){
    #     if (rs$rasterStack[["Inhabitable"]][a,b] == 0){
    #       uninhabitableCells <- append(uninhabitableCells, cellFromRowCol(rs$rasterStack, a, b))
    #     }
    #   }
    # }

    rows <- rep(1:(p / ncols), each = ncols)
    cols <- rep(1:ncols, times = (p / ncols))

    irow <- matrix(rep(rows, length(rows)), nrow = p, byrow = TRUE)
    icol <- matrix(rep(cols, length(cols)), nrow = p, byrow = TRUE)
    jrow <- t(irow) # Transpose of irow matrix
    jcol <- t(icol) # Transpose of icol matrix
    d <- sqrt((irow - jrow)^2 + (icol - jcol)^2)

    if (varCovarFunc == "DBD") {
      varMatrix <- QVar * QCorrLength^d
    } else if (varCovarFunc == "Balgovind") {
      varMatrix <- QVar * (1 + (d / QCorrLength)) * exp(-d / QCorrLength)
    } else if (varCovarFunc == "Exponential") {
      varMatrix <- QVar * exp(-d / QCorrLength)
    } else if (varCovarFunc == "Gaussian") {
      varMatrix <- QVar * exp(-(d^2) / 2 * (QCorrLength^2))
    } else if (varCovarFunc == "Spherical") {
      # Note that "QCorrLength" actually refers to the radius for the
      # spherical variance-covariance function
      varMatrix <- QVar * ((3 * d) / (2 * QCorrLength) - (d^3) / (2 * (QCorrLength^3)))
      varMatrix[d >= QCorrLength] <- 0
    } else {
      stop("Invalid variance-covariance function selected. Currently supported functions are: DBD, Balgovind, Exponential, Gaussian and Spherical")
      # Error if selected variance-covariance function is invalid
    }

    Q[d < nbhd] <- varMatrix[d < nbhd]


    # for (i in 1:p){
    #   #if (!(i %in% uninhabitableCells)) {
    #     icol <- ceiling(i/ncols)
    #     message(paste("icol =", icol))
    #     irow <- (i-1)%%ncols + 1
    #     message(paste("irow =", irow))
    #     for (j in 1:p) {
    #       #if (!(i %in% uninhabitableCells)) {
    #         jcol <- ceiling(j/ncols)
    #         message(paste("jcol =", jcol))
    #         jrow <- (j-1)%%ncols + 1
    #         message(paste("jrow =", jrow))
    #         d <- sqrt((irow-jrow)^2 + (icol - jcol)^2)
    #         message(d)
    #
    #         if (d < nbhd) {
    #           if (varCovarFunc == "DBD"){
    #             val <- QVar*QCorrLength^d
    #           } else if (varCovarFunc == "Balgovind"){
    #             val <- QVar*(1 + (d/QCorrLength))*exp(-d/QCorrLength)
    #           } else if (varCovarFunc == "Exponential"){
    #             val <- QVar*exp(-d/QCorrLength)
    #           } else if (varCovarFunc == "Gaussian"){
    #             val <- QVar*exp(-(d^2)/2*(QCorrLength^2))
    #           } else if (varCovarFunc == "Spherical") {
    #             # Note that "QCorrLength" actually refers to the radius for the
    #             # spherical variance-covariance function
    #             if (d < QCorrLength) {
    #               val <- QVar*((3*d)/(2*QCorrLength) - (d^3)/(2*(QCorrLength^3)))
    #             } else {
    #               val <- 0
    #             }
    #           } else {
    #             stop('Invalid variance-covariance function selected. Currently supported functions are: DBD, Balgovind, Exponential, Gaussian and Spherical')
    #             # Error if selected variance-covariance function is invalid
    #           }
    #          Q[i,j] = val
    #          message(val)
    #          }
    #       #}
    #     }
    #   #}
    # }

    diag(Q) <- ifelse(diag(Q) == 0, QVar, diag(Q))

    message(dim(Q))
    QFull <- Q

    if (states_observable == 2) {
      Q0 <- matrix(0, p, p)

      Qtop <- cbind(QFull, Q0)
      # message(dim(Qtop))

      Qbottom <- cbind(Q0, QFull)
      # message(dim(Qbottom))

      QFull <- rbind(Qtop, Qbottom)
    }

    message(dim(QFull))
    message(QFull[1:5, 1:5])

    return(list("Q" = Q, "QFull" = QFull))
  }

  ##################################################################
  ## Call the above function to produce the "original" Q matrix.  ##
  ##################################################################
  rs <- createRasterStack(
    selectedCountry = COD,
    rasterAgg = rasterAggregationFactor,
    susceptibleLayer = createSusceptibleLayer(COD, rasterAgg = rasterAggregationFactor)
  )

  options("max.print" = 3)

  print("Generating results of neighbourhood of zero with original code...")
  results.nbhd.zero <- lapply(
    X = variableCovarianceFunctions,
    function(x) {
      suppressMessages(genQ(
        varCovarFunc = x,
        nrows = nrow(rs$rasterStack),
        ncols = ncol(rs$rasterStack),
        nbhd = 0,
        QVar = QVar,
        QCorrLength = QCorrLength,
        states_observable = statesObservable
      ))
    }
  )

  print("Generating results of neighbourhood of one with original code...")
  results.nbhd.one <- lapply(
    X = variableCovarianceFunctions,
    function(x) {
      suppressMessages(genQ(
        varCovarFunc = x,
        nrows = nrow(rs$rasterStack),
        ncols = ncol(rs$rasterStack),
        nbhd = 1,
        QVar = 0,
        QCorrLength = 0,
        states_observable = statesObservable
      ))
    }
  )

  max.print.old <- getOption("max.print")

  list(
    No.neighborhood = results.nbhd.zero,
    Moore.neighborhood = results.nbhd.one
  )
}

## microbenchmark::microbenchmark of 100 runs of the following assignment of a
## compound expression
## milliseconds
##       min       lq    mean   median       uq      max neval
##  249.4146 274.3294 317.667 282.7732 311.4679 516.9176   100
newCodeResults <- {
  layers <- getSVEIRD.SpatRaster(subregionsSpatVector,
    susceptibleSpatRaster,
    aggregationFactor = rasterAggregationFactor
  )

  results.nbhd.zero <- lapply(
    X = variableCovarianceFunctions,
    function(x) {
      forecastError.cov(
        layers,
        variableCovarianceFunction = x,
        forecastError.cov.sdBackground = QVar,
        forecastError.cor.length = QCorrLength,
        neighbourhood = 0,
        compartmentsReported = 2
      )
    }
  )

  results.nbhd.one <- lapply(
    X = variableCovarianceFunctions,
    function(x) {
      forecastError.cov(
        layers,
        variableCovarianceFunction = x,
        forecastError.cov.sdBackground = QVar,
        forecastError.cor.length = QCorrLength,
        neighbourhood = 1,
        compartmentsReported = 2
      )
    }
  )

  list(
    No.neighborhood = results.nbhd.zero,
    Moore.neighborhood = results.nbhd.one
  )
}

test_that("Forecast error covariance (Q) matrix is correct", {
  ## NOTE: the objects must be available and have the correct names for the test
  ## to proceed. If either object doesn't have the correct names, then something
  ## is obviously wrong.
  expect_named(originalCodeResults, c("No.neighborhood", "Moore.neighborhood"))
  expect_named(newCodeResults, c("No.neighborhood", "Moore.neighborhood"))
})
