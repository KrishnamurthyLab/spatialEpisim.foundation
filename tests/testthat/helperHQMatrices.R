subregionsSpatVector <-
  terra::vect(
           system.file(
             "extdata",
             ## COD: Nord-Kivu and Ituri (Democratic Republic of Congo)
             "subregionsSpatVector",
             package = "spatialEpisim.foundation",
             mustWork = TRUE
           )
         )
susceptibleSpatRaster <-
  terra::rast(
           system.file(
             "extdata",
             "susceptibleSpatRaster.tif", # Congo population
             package = "spatialEpisim.foundation",
             mustWork = TRUE
           )
         )
data("healthZonesCongo", package = "spatialEpisim.foundation")
data("initialInfections.fourCities", package = "spatialEpisim.foundation")

suppressPackageStartupMessages({
  require("Matrix")
  require("sp")
  require("raster")
  require("countrycode")
  require("terra")
})

createSusceptibleLayer <- function(selectedCountry, rasterAgg = 0) {
    #----------------------------------------------------------------#
    # Source 1: WorldPop UN-Adjusted Population Count GeoTIFF raster #
    #----------------------------------------------------------------#

    inputISO <- countrycode(selectedCountry, origin = "country.name", destination = "iso3c") # Converts country name to ISO Alpha
    inputISOLower <- tolower(inputISO)

    url <- paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/2020/", inputISO, "/", inputISOLower, "_ppp_2020_1km_Aggregated_UNadj.tif")

    tifFileName <- basename(url) # name of the .tif file

    ## NOTE: this is modified so that the test works, but it's just providing a
    ## different root than what was originally written.
    tifFolder <- test_path("tif/") # .tif files should be stored in local tif/ folder

    if (!file.exists(paste0(tifFolder, tifFileName))) {
      download.file(url, paste0(tifFolder, tifFileName), mode = "wb")
    }

    WorldPop <- rast(paste0(tifFolder, tifFileName))

    # Gives the five number summary
    message(summary(values(WorldPop)))

    # Number of cells that have an NA value
    message(sum(is.na(values(WorldPop))))

    # message(as.raster(WorldPop))

    WorldPop <- replace(WorldPop, is.na(WorldPop), 0) # Delete this line for clear plot. Check!!!

    # WorldPop <- terra::rast(paste0(tifFolder, tifFileName))
    # Use the above line if fully switching over to terra R package completely
    # Note: rasterBasePlot.R was developed with the terra::rast()
    # Error in (function (classes, fdef, mtable) :
    #             unable to find an inherited method for function ‘classify’ for signature ‘"RasterLayer"’

    # message(WorldPop)
    # message(nrow(WorldPop))
    # message(ncol(WorldPop))
    # message(ncell(WorldPop))
    # message(res(WorldPop))
    # message(ext(WorldPop))

    if (rasterAgg == 0 || rasterAgg == 1) {
      Aggregated <- WorldPop
    } else {
      Aggregated <- aggregate(WorldPop, fact = c(rasterAgg, rasterAgg), fun = sum, na.rm = TRUE)
    }

    # message(Susceptible)

    returnList <- list("Susceptible" = WorldPop, "Aggregated" = Aggregated, "nRows" = nrow(WorldPop), "nCols" = ncol(WorldPop), "nCells" = ncell(WorldPop))

    return(returnList)
  }

createRasterStack <- function(selectedCountry, rasterAgg, isCropped = F, level1Names = NULL, susceptibleLayer) {
    inputISO <- countrycode(selectedCountry, origin = "country.name", destination = "iso3c") # Converts country name to ISO Alpha

    Susceptible <- susceptibleLayer$Aggregated



    # message(Susceptible)

    #---------------------------------------#
    # Source 2: From GADM: Level1Identifier #
    #---------------------------------------#

    gadmFileName <- paste0("gadm36_", inputISO, "_1_sp.rds") # name of the .rds file

    gadmFolder <- test_path("gadm/") # .rds files should be stored in local gadm/ folder

    # message(paste0(gadmFolder, gadmFileName))

    Level1Identifier <- readRDS(paste0(gadmFolder, gadmFileName))

    # message(Level1Identifier)
    # message(Level1Identifier$NAME_1) # List of all states/provinces/regions

    if (isCropped) {
      # positions <- which(Level1Identifier$NAME_1 %in% level1Names)  # Determines the position of which indices are TRUE.
      # message(positions)

      Level1Identifier <- Level1Identifier[which(Level1Identifier$NAME_1 %in% level1Names), ]
      # message(Level1Identifier) # It is a SpatialPolygonsDataFrame

      # Level1Identifier <- Level1Identifier[Level1Identifier$NAME_1 %in% level1Names, ]
      # message(Level1Identifier)

      # message(which(Level1Identifier$NAME_1 %in% level1Names)) # Assigns 1, 2, ...
      # message(table(Level1Identifier$NAME_1 %in% level1Names))
    }

    Level1Identifier <- vect(Level1Identifier)
    crs(Level1Identifier) <- crs(Susceptible, proj = TRUE)
    # convert Level1Identifier into a spatVector

    Level1Raster <- crop(Level1Identifier, Susceptible)
    # message(Level1Raster) # It is still a SpatialPolygonsDataFrame

    Level1Raster <- rast(Level1Identifier, resolution = res(Susceptible)[1])
    # message(Level1Raster) # It is now a RasterLayer
    # message(values(Level1Raster))

    Level1Raster <- rasterize(Level1Identifier, Level1Raster)
    # message(Level1Raster) # It is now a RasterLayer
    # message(values(Level1Raster))

    # message(extent(Level1Raster))
    # Level1Raster <- replace(Level1Raster, is.na(Level1Raster), 0)

    # message(table(values(Level1Raster)))
    # message(freq(Level1Raster))

    # Resampling methods
    # "ngb": Nearest-neighbor; assigns the value of the nearest cell
    # "bilinear": Bilinear interpolation; assigns a weighted average of the four nearest cells (the default)

    # Level1Raster <- round(resample(Level1Raster, Susceptible, method = "bilinear"))
    # Level1Raster <-  round(resample(Level1Raster, Susceptible, method = "ngb", fun ='modal'))

    Level1Raster <- resample(Level1Raster, Susceptible, method = "near")

    values(Level1Raster) <- ifelse(values(Level1Raster) > 0, values(Level1Raster), 0) # Refill the rasterLayer with 0, 1, 2, 3, ....
    # message(table(values(Level1Raster)))
    # message(freq(Level1Raster))

    # Level1Raster <- replace(Level1Raster, values(Level1Raster) < 0, 0)
    # Unless you are using method = "ngb" the above line is needed for some countries.
    # The other method is called "bilinear"

    Level1Raster <- replace(Level1Raster, is.na(Level1Raster), 0)
    # message(table(values(Level1Raster)))
    # message(freq(Level1Raster))

    # Background: Aggregating typically an entire column or an entire row or both worth of NAs are added to the Level1Raster
    # NOTE: If rasterAgg = 0 or 1, no NAs are added.
    values(Susceptible) <- ifelse(values(Level1Raster) > 0, values(Susceptible), 0) # Fill the Susceptible Layer with either a 0 or it's actual susceptible.

    Inhabitable <- Vaccinated <- Exposed <- Infected <- Recovered <- Dead <- Susceptible

    values(Vaccinated) <- values(Exposed) <- values(Infected) <- values(Recovered) <- values(Dead) <- 0 # Fill the entire rasterLayer with zeroes
    values(Inhabitable) <- ifelse(values(Susceptible) > 0, 1, 0) # Fill the rasterLayer with either a 0 or 1.

    inhabitableTrim <- trim(Inhabitable, value = 0)
    # message(extent(inhabitableTrim))

    # inhabitableRows <- inhabitableCols <-
    #
    # # message(Inhabitable[1][1])
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
    # message(minRow)
    # message(maxRow)
    # message(minCol)
    # message(maxCol)

    # message("isCropped selected")

    # message(table(as.matrix(Inhabitable)))
    # message(table(as.matrix(Vaccinated)))

    # message(freq(Level1Raster)) # Frequency table of the values of a RasterLayer.
    # message(freq(Inhabitable))

    # message(dim(Level1Raster)); message(dim(Susceptible))
    # message(res(Level1Raster)); message(res(Susceptible))
    # message(origin(Level1Raster)); message(origin(Susceptible))

    # message(extent(Level1Raster))
    # message(extent(Susceptible))
    rasterStack <- c(Susceptible, Vaccinated, Exposed, Infected, Recovered, Dead, Inhabitable, Level1Raster)

    if (isCropped) {
      rasterStack <- crop(rasterStack, inhabitableTrim)
    }

    names(rasterStack) <- c("Susceptible", "Vaccinated", "Exposed", "Infected", "Recovered", "Dead", "Inhabitable", "Level1Raster")

    # message(rasterStack)

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
