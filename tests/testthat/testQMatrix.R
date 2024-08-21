## These variables are treated as constants to be used throughout the testing
## conducted in this file.
COD <- "Democratic Republic of Congo"
QVar <- 0
QCorrLength <- 0
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
  library(countrycode)
  library(terra)

  createSusceptibleLayer <- function(selectedCountry, rasterAgg = 0) {
    #----------------------------------------------------------------#
    # Source 1: WorldPop UN-Adjusted Population Count GeoTIFF raster #
    #----------------------------------------------------------------#

    inputISO <- countrycode(selectedCountry, origin = "country.name", destination = "iso3c") # Converts country name to ISO Alpha
    inputISOLower <- tolower(inputISO)

    url <- paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/2020/", inputISO, "/", inputISOLower, "_ppp_2020_1km_Aggregated_UNadj.tif")

    tifFileName <- basename(url) # name of the .tif file
    tifFolder <- "tif/" # .tif files should be stored in local tif/ folder

    if (!file.exists(paste0(tifFolder, tifFileName))) {
      download.file(url, paste0(tifFolder, tifFileName), mode = "wb")
    }

    WorldPop <- rast(paste0(tifFolder, tifFileName))

    # Gives the five number summary
    print(summary(values(WorldPop)))

    # Number of cells that have an NA value
    print(sum(is.na(values(WorldPop))))

    # print(as.raster(WorldPop))

    WorldPop <- replace(WorldPop, is.na(WorldPop), 0) # Delete this line for clear plot. Check!!!

    # WorldPop <- terra::rast(paste0(tifFolder, tifFileName))
    # Use the above line if fully switching over to terra R package completely
    # Note: rasterBasePlot.R was developed with the terra::rast()
    # Error in (function (classes, fdef, mtable) :
    #             unable to find an inherited method for function ‘classify’ for signature ‘"RasterLayer"’

    # print(WorldPop)
    # print(nrow(WorldPop))
    # print(ncol(WorldPop))
    # print(ncell(WorldPop))
    # print(res(WorldPop))
    # print(ext(WorldPop))

    if (rasterAgg == 0 || rasterAgg == 1) {
      Aggregated <- WorldPop
    } else {
      Aggregated <- aggregate(WorldPop, fact = c(rasterAgg, rasterAgg), fun = sum, na.rm = TRUE)
    }

    # print(Susceptible)

    returnList <- list("Susceptible" = WorldPop, "Aggregated" = Aggregated, "nRows" = nrow(WorldPop), "nCols" = ncol(WorldPop), "nCells" = ncell(WorldPop))

    return(returnList)
  }

  #------------------------#
  # Example Function Calls #
  #------------------------#

  # NOTE: The isCropped argument is unused

  # createSusceptibleLayer("Nigeria", 0)

  createRasterStack <- function(selectedCountry, rasterAgg, isCropped = FALSE, level1Names = NULL, susceptibleLayer) {
    inputISO <- countrycode(selectedCountry, origin = "country.name", destination = "iso3c") # Converts country name to ISO Alpha

    Susceptible <- susceptibleLayer$Aggregated



    # print(Susceptible)

    #---------------------------------------#
    # Source 2: From GADM: Level1Identifier #
    #---------------------------------------#

    gadmFileName <- paste0("gadm36_", inputISO, "_1_sp.rds") # name of the .rds file

    gadmFolder <- "gadm/" # .rds files should be stored in local gadm/ folder

    # print(paste0(gadmFolder, gadmFileName))

    Level1Identifier <- readRDS(paste0(gadmFolder, gadmFileName))

    # print(Level1Identifier)
    # print(Level1Identifier$NAME_1) # List of all states/provinces/regions

    if (isCropped) {
      # positions <- which(Level1Identifier$NAME_1 %in% level1Names)  # Determines the position of which indices are TRUE.
      # print(positions)

      Level1Identifier <- Level1Identifier[which(Level1Identifier$NAME_1 %in% level1Names), ]
      # print(Level1Identifier) # It is a SpatialPolygonsDataFrame

      # Level1Identifier <- Level1Identifier[Level1Identifier$NAME_1 %in% level1Names, ]
      # print(Level1Identifier)

      # print(which(Level1Identifier$NAME_1 %in% level1Names)) # Assigns 1, 2, ...
      # print(table(Level1Identifier$NAME_1 %in% level1Names))
    }

    Level1Identifier <- vect(Level1Identifier)
    crs(Level1Identifier) <- crs(Susceptible, proj = TRUE)
    # convert Level1Identifier into a spatVector

    Level1Raster <- crop(Level1Identifier, Susceptible)
    # print(Level1Raster) # It is still a SpatialPolygonsDataFrame

    Level1Raster <- rast(Level1Identifier, resolution = res(Susceptible)[1])
    # print(Level1Raster) # It is now a RasterLayer
    # print(values(Level1Raster))

    Level1Raster <- rasterize(Level1Identifier, Level1Raster)
    # print(Level1Raster) # It is now a RasterLayer
    # print(values(Level1Raster))

    # print(extent(Level1Raster))
    # Level1Raster <- replace(Level1Raster, is.na(Level1Raster), 0)

    # print(table(values(Level1Raster)))
    # print(freq(Level1Raster))

    # Resampling methods
    # "ngb": Nearest-neighbor; assigns the value of the nearest cell
    # "bilinear": Bilinear interpolation; assigns a weighted average of the four nearest cells (the default)

    # Level1Raster <- round(resample(Level1Raster, Susceptible, method = "bilinear"))
    # Level1Raster <-  round(resample(Level1Raster, Susceptible, method = "ngb", fun ='modal'))

    Level1Raster <- resample(Level1Raster, Susceptible, method = "near")

    values(Level1Raster) <- ifelse(values(Level1Raster) > 0, values(Level1Raster), 0) # Refill the rasterLayer with 0, 1, 2, 3, ....
    # print(table(values(Level1Raster)))
    # print(freq(Level1Raster))

    # Level1Raster <- replace(Level1Raster, values(Level1Raster) < 0, 0)
    # Unless you are using method = "ngb" the above line is needed for some countries.
    # The other method is called "bilinear"

    Level1Raster <- replace(Level1Raster, is.na(Level1Raster), 0)
    # print(table(values(Level1Raster)))
    # print(freq(Level1Raster))

    # Background: Aggregating typically an entire column or an entire row or both worth of NAs are added to the Level1Raster
    # NOTE: If rasterAgg = 0 or 1, no NAs are added.
    values(Susceptible) <- ifelse(values(Level1Raster) > 0, values(Susceptible), 0) # Fill the Susceptible Layer with either a 0 or it's actual susceptible.

    Inhabitable <- Vaccinated <- Exposed <- Infected <- Recovered <- Dead <- Susceptible

    values(Vaccinated) <- values(Exposed) <- values(Infected) <- values(Recovered) <- values(Dead) <- 0 # Fill the entire rasterLayer with zeroes
    values(Inhabitable) <- ifelse(values(Susceptible) > 0, 1, 0) # Fill the rasterLayer with either a 0 or 1.

    inhabitableTrim <- trim(Inhabitable, value = 0)
    # print(extent(inhabitableTrim))

    # inhabitableRows <- inhabitableCols <-
    #
    # # print(Inhabitable[1][1])
    # #
    # # for (i in seq_len(nrow(Inhabitable))) {
    # #   for (j in seq_len(ncol(Inhabitable))){
    # #     if (Inhabitable[i][j] == 1) {
    # #       append(inhabitableRows, i)
    # #       append(inhabitableCols, j)
    # #     }
    # #   }
    # # }
    # #
    # # minRow = min(inhabitableRows)
    # # maxRow = max(inhabitableRows)
    # # minCol = min(inhabitableCols)
    # # maxCol = max(inhabitableCols)
    # #
    # print(minRow)
    # print(maxRow)
    # print(minCol)
    # print(maxCol)

    # print("isCropped selected")

    # print(table(as.matrix(Inhabitable)))
    # print(table(as.matrix(Vaccinated)))

    # print(freq(Level1Raster)) # Frequency table of the values of a RasterLayer.
    # print(freq(Inhabitable))

    # print(dim(Level1Raster)); print(dim(Susceptible))
    # print(res(Level1Raster)); print(res(Susceptible))
    # print(origin(Level1Raster)); print(origin(Susceptible))

    # print(extent(Level1Raster))
    # print(extent(Susceptible))
    rasterStack <- c(Susceptible, Vaccinated, Exposed, Infected, Recovered, Dead, Inhabitable, Level1Raster)

    if (isCropped) {
      rasterStack <- crop(rasterStack, inhabitableTrim)
    }

    names(rasterStack) <- c("Susceptible", "Vaccinated", "Exposed", "Infected", "Recovered", "Dead", "Inhabitable", "Level1Raster")

    # print(rasterStack)

    returnList <- list(
      "rasterStack" = rasterStack,
      "Level1Identifier" = Level1Identifier,
      "selectedCountry" = selectedCountry,
      "rasterAgg" = rasterAgg,
      "nRows" = susceptibleLayer$nRows,
      "nCols" = susceptibleLayer$nCols,
      "nCells" = susceptibleLayer$nCells
    )

    return(returnList)
  }

  genQ <- function(nrows, ncols, varCovarFunc, QVar, QCorrLength, nbhd, states_observable = 2) {
    p <- nrows * ncols
    print(p)

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
    #     print(paste("icol =", icol))
    #     irow <- (i-1)%%ncols + 1
    #     print(paste("irow =", irow))
    #     for (j in 1:p) {
    #       #if (!(i %in% uninhabitableCells)) {
    #         jcol <- ceiling(j/ncols)
    #         print(paste("jcol =", jcol))
    #         jrow <- (j-1)%%ncols + 1
    #         print(paste("jrow =", jrow))
    #         d <- sqrt((irow-jrow)^2 + (icol - jcol)^2)
    #         print(d)
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
    #          print(val)
    #          }
    #       #}
    #     }
    #   #}
    # }

    diag(Q) <- ifelse(diag(Q) == 0, QVar, diag(Q))

    print(dim(Q))
    QFull <- Q

    if (states_observable == 2) {
      Q0 <- matrix(0, p, p)

      Qtop <- cbind(QFull, Q0)
      # print(dim(Qtop))

      Qbottom <- cbind(Q0, QFull)
      # print(dim(Qbottom))

      QFull <- rbind(Qtop, Qbottom)
    }

    print(dim(QFull))
    print(QFull[1:5, 1:5])

    return(list("Q" = Q, "QFull" = QFull))
  }

  ##################################################################
  ## Call the above functions to produce the "original" Q matrix. ##
  ##################################################################
  rs <- createRasterStack(
    selectedCountry = COD,
    rasterAgg = rasterAggregationFactor,
    susceptibleLayer = createSusceptibleLayer(COD, rasterAgg = rasterAggregationFactor)
  )

  options("max.print" = 1)

  results.nbhd.zero <- lapply(
    X = variableCovarianceFunctions,
    function(x) {
      genQ(
        varCovarFunc = x,
        nrows = nrow(rs$rasterStack),
        ncols = ncol(rs$rasterStack),
        nbhd = 0,
        QVar = QVar,
        QCorrLength = QCorrLength,
        states_observable = statesObservable
      )
    }
  )

  results.nbhd.one <- lapply(
    X = c("DBD", "Balgovind", "Exponential", "Gaussian", "Spherical"),
    function(x) {
      genQ(
        varCovarFunc = x,
        nrows = nrow(rs$rasterStack),
        ncols = ncol(rs$rasterStack),
        nbhd = 1,
        QVar = 0,
        QCorrLength = 0,
        states_observable = 2
      )
    }
  )

  max.print.old <- getOption("max.print")

  list(
    No.neighborhood = results.nbhd.zero,
    Moore.neighborhood = results.nbhd.one
  )
}

newCodeResults <- {
  subregionsSpatVector <- terra::vect(
    system.file(
      "extdata",
      ## COD: Nord-Kivu and Ituri (Democratic Republic of Congo)
      "subregionsSpatVector",
      package = "spatialEpisim.foundation",
      mustWork = TRUE
    )
  )
  susceptibleSpatRaster <- terra::rast(
    system.file(
      "extdata",
      "susceptibleSpatRaster.tif", # Congo population
      package = "spatialEpisim.foundation",
      mustWork = TRUE
    )
  )
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

test_that("Forecast error covariance matrix, Q, is correct", {
  ## Access oldCodeResults
  ## Access newCodeResults
})
