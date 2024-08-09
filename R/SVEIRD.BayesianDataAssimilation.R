##' Set the options for this package; primarily set path
##' `spatialEpisim.foundation.datapath` used for storing data used by this
##' package which is downloaded from the Web.
##'
##' The option `spatialEpisim.foundation.datapath` is only set to its default
##' value if it is not already set.
##' @title Cause side-effects when loading this package
##' @author Bryce Carson
##' @param libname The name of the library in which the package is stored.
##' @param pkgname The name of the loaded package.
.onLoad <- function(libname, pkgname) {
  if(is.null(getOption("spatialEpisim.foundation.datapath"))) {
    options("spatialEpisim.foundation.datapath" = path.expand("~/spatialEpisim.foundation.data"))
    if(!dir.exists(getOption("spatialEpisim.foundation.datapath")))
      dir.create(getOption("spatialEpisim.foundation.datapath"),
                 recursive = TRUE,
                 showWarnings = TRUE)
    else
      warning(sprintf("%s already exists, but option `%s` was unset.",
                      getOption("spatialEpisim.foundation.datapath"),
                      "spatialEpisim.foundation.datapath"))
  }
}

##' Download 1 km aggregated UN-adjusted population count spatial data from 2020
##' through the WorldPop servers.
##'
##' This function only downloads the file if necessary (if it doesn't exist at
##' path).
##'
##' ppp is the code used by WorldPop to mean population estimate for population
##' counts; similar data, but which grids population density is coded _pd_
##' rather than _ppp_. 1km means 1 km aggregation of the data. The description
##' WorldPop gives for this data is "individual countries 2000-2020 UN adjusted
##' aggregated to 1km resolution using 100m resolution population count
##' datasets". See: https://hub.worldpop.org/doi/10.5258/SOTON/WP00671.
##' @title Download WorldPop 2020 data
##' @param countryCodeISO3C The uppercase ISO three character code a recognized
##'   country.
##' @param folder The destination folder the downloaded data will be stored in
##' @returns An absolute path where the data was downloaded, or the path at which
##'   the file already existed
##' @author Bryce Carson
##' @export
##' @examples
##' downloadWorldPopData(countryCodeISO3C = "COD", folder = tempdir())
downloadWorldPopData <- function(countryCodeISO3C, folder = getOption("spatialEpisim.foundation.datapath")) {
  ## NOTE: Code borrowed from `httr2:::dots_to_path`; see citation("httr2")
  dotsToPath <- \(...) sub("^([^/])", "/\\1", paste(c(...), collapse = "/"))
  ## Construct the path to the data on data.worldpop.org
  urlPath <-
    c("GIS",
      "Population",
      "Global_2000_2020_1km_UNadj",
      "2020",
      countryCodeISO3C,
      basename = sprintf("%s_ppp_2020_1km_Aggregated_UNadj.tif",
                         tolower(countryCodeISO3C))) %>%
    dotsToPath()
  url <- httr2::url_build(structure(list(scheme = "https",
                                  hostname = "data.worldpop.org",
                                  path = urlPath),
                             class = "httr2_url"))

  ## Download the GeoTIFF file if it doesn't already exist.
  if (!file.exists(file.path(getOption("spatialEpisim.foundation.datapath"), basename(urlPath)))) {
    download.file(url,
                  file.path(getOption("spatialEpisim.foundation.datapath"),
                            basename(urlPath)),
                  mode = "wb")
  }

  file.path(getOption("spatialEpisim.foundation.datapath"),
            basename(urlPath))
}

##' Returns a [terra::SpatVector] for the requested country on demand, either
##' retrieving the data from disk if it has been downloaded before, or
##' downloading the SpatVector to disk at the location determined by
##' [geodata::geodata_path()].
##'
##' The level one boundaries are the least granular (geographically largest)
##' administrative boundaries that countries create to subdivide themselves.
##'
##' The default resolution is one kilometer, meaning points in the vector are no
##' closer to each other than one kilometer, therefore details smaller than one
##' kilometer cannot be differentiated using this data.
##' @title Retrieve a SpatRaster of level 1 boundaries from local or remote disk
##' @param countryCodeISO3C The uppercase ISO three character code a recognized
##'   country.
##' @param folder The path where GADM data should be found or stored; passed on
##'   to geodata::gadma.
##' @returns SpatVector
##' @author Bryce Carson
##' @export
##' @examples
##' library(geodata)
##' geodata_path_backup <- geodata_path()
##' options(geodata_default_path = getOption("spatialEpisim.foundation.datapath"))
##' lvl1AdminBorders("COD")
##' options(geodata_default_path = geodata_path_backup)
##'
##' \dontrun{
##' ## It's recommended to set geodata_default_path in your .Rprofile rather
##' ## than explicitly provide the path.
##' lvl1AdminBorders("COD", file.path("data", "gadm"))
##' }
lvl1AdminBorders <- function(countryCodeISO3C, folder = geodata::geodata_path()) {
  geodata::gadm(
    country = countryCodeISO3C,
    level = 1,
    path = folder,
    version = "3.6",
    resolution = 1
  )
}

##' Create a named RasterLayer object useful for spatiotemporal epidemic
##' compartmental modelling.
##' @title Create a Susceptible-component RasterLayer
##' @param countryCodeISO3C The uppercase ISO three character code a recognized
##'   country.
##' @param folder Passed on to method downloadWorldPopData
##' @returns A SpatRaster of WorldPop population count data with the name
##'   Susceptible, with all NAs replaced by zeros.
##' @author Bryce Carson
##' @author Michael Myer
##' @author Ashok Krishnmaurthy
##' @export
##' @examples
##' getCountryPopulation.SpatRaster("COD")
getCountryPopulation.SpatRaster <- function(countryCodeISO3C, folder = NULL) {
  countryRaster <-
    if(is.null(folder)) {
      downloadWorldPopData(countryCodeISO3C)
    } else {
      downloadWorldPopData(countryCodeISO3C, folder)
    }

  countrySpatRaster <- terra::rast(countryRaster)
  replace(countrySpatRaster, is.na(countrySpatRaster), 0) %>%
    `names<-`("Susceptible")
}

##' Read an RDS file from the GADM folder for a given country's subregion.
##'
##' @title Get a country's GADM verison 3.6 raster
##' @param countryCodeISO3C The uppercase ISO three character code a recognized
##'   country.
##' @param level1Region The subregions of the country to crop to
##' @param folder The folder which should be searched for the GADM data, passed
##'   on to lvl1AdminBorders.
##' @returns SpatVector for the specificed country.
##' @author Bryce Carson
##' @export
##' @examples
##' library(geodata)
##' geodata_path_backup <- geodata_path()
##' options(geodata_default_path = getOption("spatialEpisim.foundation.datapath"))
##' getCountrySubregions.SpatVector("COD", c("Nord-Kivu", "Ituri"))
##' getCountrySubregions.SpatVector("CZE", "Prague")
##' getCountrySubregions.SpatVector("NGA", "Kwara")
##' options(geodata_default_path = geodata_path_backup)
getCountrySubregions.SpatVector <- function(countryCodeISO3C = "COD",
                                            level1Region = c("Nord-Kivu", "Ituri"),
                                            folder) {
  stopifnot(countryCodeISO3C %in% countrycode::codelist$iso3c)
  lvl1AdminBorders(countryCodeISO3C, folder) %>%
    subset(.$NAME_1 %in% level1Region)
}

## NOTE: This function is related to cropping "base" rasters, and plotting them
## using the Haxby colour table. For more information on the Haxby colour table,
## see the entry for Haxby in the r.colors documentation of the GRASS GIS
## software suite: https://grass.osgeo.org/grass83/manuals/r.colors.html; this
## states "relative colors for bathymetry or topography".
##' Plot a raster cropped to a specific region(s)
##'
##' As a side effect, it writes a raster to file.
##' @title Plot a raster cropped to sub-region(s) of a country
##' @author Bryce Carson
##' @author gursDhaliwal
##' @author Ashok Krishnmaurthy
##' @author Michael Myer
##' @author Thomas White
##' @param subregions a SpatRaster object for subregions, created with
##'   getCountrySubregions.SpatVector
##' @param susceptible SpatRaster or SpatVector to be masked by GADM data
##'   of the country
##' @export
##' @keywords internal
##' @examples
##' library(lattice)
##' subregionsSpatVector <- terra::vect(
##'   system.file(
##'     "extdata",
##'     ## COD: Nord-Kivu and Ituri (Democratic Republic of Congo)
##'     "subregionsSpatVector",
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' susceptibleSpatRaster <- terra::rast(
##'   system.file(
##'     "extdata",
##'     "susceptibleSpatRaster.tif", # Congo population
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' masked <- maskAndClassifySusceptibleSpatRaster(subregionsSpatVector,
##'                                                susceptibleSpatRaster)
##' terra::plot(masked)
##'
##' \dontrun{
##' maskAndClassifySusceptibleSpatRaster(getSubregion("CZE", "Prague"),
##'                                      getCountryPopulation.SpatRaster("CZE"))
##'
##' maskAndClassifySusceptibleSpatRaster(getSubregion("NGA", "Lagos"),
##'                                      getCountryPopulation.SpatRaster("NGA"))
##'
##' maskAndClassifySusceptibleSpatRaster(getSubregion("COD", "Ituri"),
##'                                      getCountryPopulation.SpatRaster("COD"))
##'
##' maskAndClassifySusceptibleSpatRaster(getSubregion("COD", c("Nord-Kivu", "Ituri")),
##'                                      getCountryPopulation.SpatRaster("COD"))
##' }
maskAndClassifySusceptibleSpatRaster <- function(subregions, susceptible) {
  terra::crs(subregions) <- terra::crs(susceptible, proj = TRUE)

  ## NOTE: these two values were betwixt the call of rast and classify, which
  ## are only piped to remove the use of an unnecssary assignment whilst these
  ## interjecting values are unused and therefore don't need to be calculated.
  ## Rather than interjecting and elongating the pipe's line count, they are
  ## moved here (before the pipe), so only a relevant comment interjects in the
  ## pipeline.

  ## Mask the susceptible SpatRaster by the subregions SpatRaster
  terra::crop(susceptible, subregions, mask = TRUE) %>%
  ## NOTE: this is refactored, but untested. NOTE: this seems to classify the
  ## values as 1:7, if they have a value between the bin lower limits; that is,
  ## bin the values into classes delineated by the given vector of minimums for
  ## the bins.
  terra::classify(c(0, 10, 25, 50, 100, 250, 1000, 10000)) %>%
    "levels<-"(levels(.)[[1]])
}

## FIXME: refactor the return value so that the Inhabited layer doesn't
## continually need to be ignored through subsetting the layers.
##' @title Create a SpatRaster of SVEIRD model compartment raster data
##' @description Create a list of SpatRaster objects, one for each component in
##'   an SVEIRD epidemic model.
##' @details The SpatRaster objects for the VEIRD components are empty, while
##'   the Inhabited SpatRaster is a binary classification on habitation of land
##'   area.
##' @param subregions a SpatVector object of subregions used to crop the SVEIRD SpatRaster
##' @param susceptible a RasterLayer, pre-aggregated if that is wished.
##' @param aggregationFactor the number of cells in any direction to aggregate
##'   together into one, in the susceptible ratser, after masking with the
##'   subregions vector, before creating the list
##' @returns a SpatRaster, with layers for the SVERID components and an
##'   additional layer classifying the Inhabited status of a cell
##' @author Bryce Carson
##' @author Michael Myer
##' @author Ashok Krishnmaurthy
##' @export
##' @examples
##' subregionsSpatVector <- terra::vect(
##'   system.file(
##'     "extdata",
##'     ## COD: Nord-Kivu and Ituri (Democratic Republic of Congo)
##'     "subregionsSpatVector",
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' susceptibleSpatRaster <- terra::rast(
##'   system.file(
##'     "extdata",
##'     "susceptibleSpatRaster.tif", # Congo population
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' getSVEIRD.SpatRaster(subregionsSpatVector, susceptibleSpatRaster,
##'                      aggregationFactor = 35)
##'
##' ## Omitting the aggregation factor argument will prevent aggregation. An
##' ## aggregation factor of zero or one is meaningless and will produce an
##' ## error.
##' getSVEIRD.SpatRaster(subregionsSpatVector, susceptibleSpatRaster)
getSVEIRD.SpatRaster <- function(subregions, susceptible, aggregationFactor = NULL) {
  susceptible <- maskAndClassifySusceptibleSpatRaster(subregions, susceptible)

  if (!is.null(aggregationFactor)) {
    susceptible <- terra::aggregate(susceptible, fact = aggregationFactor, fun = sum)
  }

  terra::values(susceptible)[terra::values(susceptible) < 0] <- 0
  Inhabited <- susceptible
  terra::values(Inhabited)[terra::values(Inhabited) > 0] <- 1
  terra::values(Inhabited)[terra::values(Inhabited) < 1] <- 0

  ## NOTE: It is faster to use fun = NA, rather than fun = 0, because
  ## castSeedDataQueensNeighbourhood, for example, will not perform calculations,
  ## rather than multiple many cells by zero to no effect; it also produces
  ## clearer plots without values outside the bounds of the borders of the
  ## spatial region of the geographical region represented in the raster.
  empty <- terra::init(susceptible, fun = NA)
  c(susceptible, rep(empty, 5), Inhabited) %>%
    "names<-"(c("Susceptible",
                "Vaccinated",
                "Exposed",
                "Infected",
                "Recovered",
                "Dead",
                "Inhabited"))
}

##' @title Average Euclidean Distance
##' @description Measure the exposure influence upon susceptible individuals at
##'   spatial position (x, y) of the infectious individuals at a spatial
##'   position (u, v).
##' @details The 2020 United Nations-adjusted population count raster data, at
##'   the 30-arc second resolution, consists of hundreds of thousands (and for
##'   some countries millions) of grid cells. This function is used to calculate
##'   the average Euclidean distance from one cell to others, adjusting for
##'   raster aggregation.
##'
##'   The mathematical modelling of human or non-human mobility patterns is
##'   non-trivial. This function is a limited implementation of the effective
##'   area for only human mobility, with distances travelled per day (lambda)
##'   measured in kilometers.
##'
##'   NOTE: Our weight function does not take the same arguments as shown in the
##'   slideshow: (w x y u v); The term "effective area" comes from the fact that
##'   if the kernel is constant on a finite disk (and zero outside it), then the
##'   formula due to Bolker (1999) gives the area of the disk.
##'
##'   See the article titled *A clarification of transmission terms in
##'   host-microparasite models by Begon et al.* (2002).
##'
##'   The count of infected persons, or incidence, raster data provided to the
##'   caller of this function, i.e. [transmissionLikelihoodWeightings], may be
##'   aggregated by a given factor greater than one, `aggregationFactor`; it is
##'   passed regardless.
##' @param lambda the average dialy movement distance of an individual (in
##'   kilometers).
##' @param aggregationFactor the factor of aggregation applied to the raster
##'   data mentioned in the function details (expressed in kilometers).
##' @param epsilon a rounding error correction (a real number; by default,
##'   zero).
##' @returns a matrix of the average Euclidean distances
##' @export
##' @keywords internal
##' @author Bryce Carson
##' @author Thomas White
##' @examples
##' averageEuclideanDistance(lambda = 15) # No raster aggregation
##' averageEuclideanDistance(lambda = 15, aggregationFactor = 35)
##'
##' lattice::levelplot(as.array(averageEuclideanDistance(15)))
##' lattice::levelplot(as.array(averageEuclideanDistance(100, 35)))
averageEuclideanDistance <-
  function(lambda, aggregationFactor = 1, epsilon = 0) {
    radius <-
      if (lambda <= aggregationFactor)
        round((lambda - aggregationFactor) / aggregationFactor + epsilon) + 1
      else
        lambda + aggregationFactor

    stopifnot(radius %in% seq(floor(radius), ceiling(radius)))

    avg.euc.dist <- function(i, j) {
      exp(-sqrt(sum((c(i, j) - c(radius + 1, radius + 1))^2)) / lambda)
    }
    len <- seq_len(1 + radius * 2)
    weights <-
      dplyr::mutate(
               tidyr::expand(tibble::tibble(i = len, j = len), i, j),
               averageEuclideanDistance = purrr::map2_dbl(i, j, avg.euc.dist)
             ) %>%
      dplyr::select(averageEuclideanDistance) %>%
      unlist(use.names = FALSE) %>%
      base::matrix(byrow = TRUE, ncol = sqrt(length(.)))
    stopifnot(dim(weights)[1] == dim(weights)[2])
    stopifnot(dim(weights)[1] %% 2 == 1)
    return(weights)
  }

##' @title Weighted Sums
##' @description Calculate a matrix of weights respecting human mobility
##'   patterns.
##' @details The pattern of human mobility used is described in a slideshow
##'   here:
##'   https://docs.google.com/presentation/d/1_gqcEh4d8yRy22tCZkU0MbGYSsGu3Djh/edit?usp=sharing&ouid=102231457806738400087&rtpof=true&sd=true.
##' @inherit SVEIRD.BayesianDataAssimilation
##' @param infections a matrix of the count of infections per aggregate area in
##'   a raster of terrestrial data.
##' @returns a matrix of weightings for the calculation of the proportionOfSusceptible of
##'   exposed individuals who will become infectious.
##' @author Bryce Carson
##' @author Thomas White
##' @export
##' @keywords internal
##' @examples
##' library(terra)
##' subregionsSpatVector <- vect(
##'   system.file(
##'     "extdata",
##'     ## COD: Nord-Kivu and Ituri (Democratic Republic of Congo)
##'     "subregionsSpatVector",
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' susceptibleSpatRaster <- rast(
##'   system.file(
##'     "extdata",
##'     "susceptibleSpatRaster.tif", # Congo population
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' layers <- getSVEIRD.SpatRaster(subregionsSpatVector,
##'                                susceptibleSpatRaster,
##'                                aggregationFactor = 35)
##' transmissionLikelihoodWeightings(layers$Infected, 15, 35)
transmissionLikelihoodWeightings <-
  function(infections, lambda, aggregationFactor) {
    terra::focal(infections, w = averageEuclideanDistance(lambda, aggregationFactor))
  }

##' @title Linear (Forward) Interpolation Operator matrices for one or two state
##'   vectors
##' @details Create a linear forward interpolation operator matrix with as many
##'   columns as cells in the SpatRaster layers, and as many rows as health
##'   zones for which coordinates are provided (when creating a matrix for use
##'   with one state vector). The matrix is either a trivial matrix as
##'   described, or two partitions in a sparse, block diagonal matrix.
##'
##' When used to creat an interpolation matrix for two state vectors, the result
##'   is a block diagonal matrix of identical partitions, with each partition as
##'   described for one state vector.
##' @param layers The SpatRaster object with SVEIRD compartment layers, and a
##'   layer classifying habitation. Created with the getSVEIRD.SpatRaster
##'   function.
##' @param healthZoneCoordinates a table of values giving the latitude and
##'   longitude coordinates for health zones in the country of interest. See the
##'   description of the healthZoneCoordinates argument in the function
##'   SVEIRD.BayesianDataAssimilation for more details.
##'
##'   \preformatted{
##'     HealthZone    Latitude    Longitude
##'     Alimbongo     -0.365515   29.1911818
##'     Beni          0.49113     29.47306
##'     Biena         0.57923     29.115633
##'     Butembo       0.140692    29.335014
##'     Goma          -1.658271   29.220132
##'     Kalunguta     0.323085    29.354677
##'     Katwa         0.116985    29.411838
##'     Kayna         -0.603936   29.174066
##'     Kyondo        -0.005622   29.408813
##'     Lubero        -0.15633    29.24057
##'     Mabalako      0.461257    29.210687
##'     Manguredjipa  0.353433    28.765479
##'     Masereka      -0.133333   29.333333
##'     Musienene     0.04022     29.26246
##'     Mutwanga      0.32893     29.74576
##'     Nyiragongo    -1.516667   29.25
##'     Oicha         0.698681    29.518834
##'     Pinga         -0.9830504  28.687911
##'     Vuhovi        0.1416      29.4075
##'     Ariwara       3.136873    30.706615
##'     Bunia         1.566667    30.25
##'     Komanda       1.367496    29.774322
##'     Lolwa         1.352969    29.498455
##'     Mambasa       1.359731    29.029226
##'     Mandima       1.35551     29.08173
##'     Nyakunde      1.431271    30.029283
##'     Rwampara      1.4053      30.3449
##'     Tchomia       1.4412      30.4845
##'   }
##' @param compartmentsReported either 1 or 2. Previously identified as
##'   states_observable, this is the count of compartments that are reported on
##'   and which will have data assimilated; if it is 2, the matrix is a block
##'   diagonal matrix.
##' @returns a linear forward interpolation operator matrix used in
##'   SVEIRD.BayesianDataAssimilation simulations.
##' @author Bryce Carson
##' @author Ashok Krishnmaurthy
##' @author Michael Myer
##' @export
##' @keywords internal
##' @examples
##' subregionsSpatVector <- terra::vect(
##'   system.file(
##'     "extdata",
##'     ## COD: Nord-Kivu and Ituri (Democratic Republic of Congo)
##'     "subregionsSpatVector",
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' susceptibleSpatRaster <- terra::rast(
##'   system.file(
##'     "extdata",
##'     "susceptibleSpatRaster.tif", # Congo population
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' data("healthZonesCongo", package = "spatialEpisim.foundation")
##' linearInterpolationOperator(
##'   layers = getSVEIRD.SpatRaster(subregionsSpatVector,
##'                                 susceptibleSpatRaster,
##'                                 aggregationFactor = 35),
##'   healthZoneCoordinates = healthZonesCongo
##' )
##'
##' H.matrix <-
##'   linearInterpolationOperator(
##'     layers = getSVEIRD.SpatRaster(subregionsSpatVector,
##'                                   susceptibleSpatRaster,
##'                                   aggregationFactor = 35),
##'     healthZoneCoordinates = healthZonesCongo,
##'     neighbourhood.order = 1
##'   )
##' plot(as.array(H.matrix[1, ]))
##'
##' linearInterpolationOperator(
##'   layers = getSVEIRD.SpatRaster(subregionsSpatVector,
##'                                 susceptibleSpatRaster,
##'                                 aggregationFactor = 35),
##'   healthZoneCoordinates = healthZonesCongo,
##'   neighbourhood.order = 2
##' )
linearInterpolationOperator <-
  function(layers,
           healthZoneCoordinates,
           neighbourhood.order = 0,
           compartmentsReported = 1) {
    stopifnot(neighbourhood.order %in% c(0, 1, 2))
    ## TODO: add an argument to select one cell (no neighbourhood, a Moore
    ## neighbourhood, or a 2nd-order Chess Queen's neighbourhood). A second order
    ## neighbour can't be calculated with the current algorithm if the number of
    ## columns is less than five.
    if (neighbourhood.order == 2) stopifnot(ncol(layers) >= 5)
    stopifnot(compartmentsReported %in% 1:2) # NOTE: this is defensive programming; I've only implemented these cases.

    queensNeighbours <- function(order, cell, ncols) {
      stopifnot(order %in% 1:2)

      if (order == 1) {
        neighbouringCells <-
          c((cell - ncols - 1) : (cell - ncols + 1),
            cell - 1 , cell + 1,
            (cell + ncols - 1) : (cell + ncols + 1))
        stopifnot(length(neighbouringCells) == 8)
      } else if (order == 2) {
        neighbouringCells <-
          c((cell - ncols * 2 - 2) : (cell - ncols * 2 + 2),
            cell - ncols - 2 , cell - ncols + 2,
            cell - 2 , cell + 2,
            cell + ncols - 2 , cell + ncols + 2,
            (cell + ncols * 2 - 2) : (cell + ncols * 2 + 2))
        stopifnot(length(neighbouringCells) == 16)
      }

      neighbouringCells
    }

    extend.length <- 5
    layers <- terra::extend(layers, extend.length)

    ## NOTE: cells contains the index into the rasters in layers (when converted
    ## to a matrix). MAYBE FIXME: The coordinates are re-ordered as
    ## longitude-latitude, rather than latitude-longitude as they are otherwise
    ## stored; the reason is due to the generation of NAs in the resulting matrix,
    ## otherwise, according to the previous implementation.
    cells <- terra::cellFromXY(layers, terra::as.matrix(healthZoneCoordinates[, 3:2]))
    if (anyDuplicated(cells) > 0)
      warning("Raster aggregation factor is too high to differentiate between two (or more) health zones (they correspond to the same grid cell).")
    if (any(is.na(cells)))
      warning("Ignoring NAs in [cells] object corresponding to coordinates out of bounds of [layers] raster.")

    cells <- cells[!is.na(cells)]

    ## NOTE: preallocate the linear forward interpolation matrix, with
    ## dimensions q * p (health zones by cells in the SpatRaster).
    H.extended <- base::matrix(0, nrow(healthZoneCoordinates), terra::ncell(layers))

    ## NOTE: these are the weightings used for the chess queen zeroth, first,
    ## and second order neighbors. The zeroth order neighbor is the position of
    ## the queen itself.
    neighbour.weights <-
      switch(neighbourhood.order + 1, # the first of ... applies to zero, etc.
             1,
             c(2, 1) * 0.1,
             c(3, 2, 1) * 35^-1)

    for (index in seq_along(cells)) {
      H.extended[index, cells[index]] <- neighbour.weights[1]

      if (neighbourhood.order != 0) {
        neighbour.1st <- queensNeighbours(1, cells[index], terra::ncol(layers)); stopifnot(length(neighbour.1st) == 8)
        H.extended[index, neighbour.1st[neighbour.1st > 0 & neighbour.1st <= terra::ncell(layers)]] <- neighbour.weights[2]
      }

      if (neighbourhood.order == 2) {
        neighbour.2nd <- queensNeighbours(2, cells[index], terra::ncol(layers)); stopifnot(length(neighbour.2nd) == 16)
        if(anyDuplicated(c(neighbour.1st, neighbour.2nd)) > 0)
          simpleError("Duplicate cell indices among neighbours of multiple localities.")
        H.extended[index, neighbour.2nd[neighbour.2nd > 0 & neighbour.2nd <= terra::ncell(layers)]] <- neighbour.weights[3]
      }
    }

    stopifnot(dplyr::near(sum(H.extended), nrow(healthZoneCoordinates)))
    stopifnot(dplyr::near(sum(matrix(H.extended[1, ],
                                     ncol = ncol(layers),
                                     byrow = TRUE)),
                          1))

    if (compartmentsReported == 2) H.extended <- Matrix::bdiag(H.extended, H.extended)

    ## NOTE: the extended areas of the matrix are now dropped to return the matrix
    ## to the expected size for the input.
    interpolationOperatorMatrix <-
      apply(X = H.extended,
            MARGIN = 1, # apply the function to rows
            FUN =
              function(row) {
                m <- matrix(row, byrow = TRUE, ncol = ncol(layers))
                m[(extend.length + 1):(terra::nrow(m) - extend.length),
                (extend.length + 1):(terra::ncol(m) - extend.length)] %>%
                  Matrix::t() %>% # row-major order (byrow)
                  as.vector()
              }) %>% Matrix::t() # rows should be health zones

    return(interpolationOperatorMatrix)
  }

##' @description generates a block diagonal error covariance matrix with exponential decay
##' @details TODO: write the details about the implementation of this function.
##' @title Create a Q-matrix
##' @param layers The SpatRaster object with SVEIRD compartment layers, and a
##'   layer classifying habitation. Created with the getSVEIRD.SpatRaster
##'   function.
##' @param variableCovarianceFunction TODO
##' @param Q.backgroundErrorStandardDeviatio TODO
##' @param Q.characteristicCorrelationLength TODO
##' @param neighbourhood TODO
##' @param compartmentsReported TODO
##' @returns TODO
##' @author Bryce Carson
##' @author Ashok Krishnmaurthy
##' @author Michael Myer
##' @author Thomas White
##' @export
##' @keywords internal
##' @examples
##' subregionsSpatVector <- terra::vect(
##'   system.file(
##'     "extdata",
##'     ## COD: Nord-Kivu and Ituri (Democratic Republic of Congo)
##'     "subregionsSpatVector",
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' susceptibleSpatRaster <- terra::rast(
##'   system.file(
##'     "extdata",
##'     "susceptibleSpatRaster.tif", # Congo population
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' layers <- getSVEIRD.SpatRaster(subregionsSpatVector,
##'                                susceptibleSpatRaster,
##'                                aggregationFactor = 35)
##' Ituri.Q.forecastErrorCov <- Q.forecastErrorCov(layers, "DBD", 2, 0.8, 4, 2)
Q.forecastErrorCov <- function(layers,
                               variableCovarianceFunction,
                               Q.backgroundErrorStandardDeviation,
                               Q.characteristicCorrelationLength,
                               neighbourhood,
                               compartmentsReported = 2) {
  ncols <- terra::ncol(layers)
  numberOfCells <- terra::ncell(layers)
  Q <- matrix(0, numberOfCells, numberOfCells) # (rows * columns)^2 sparse
  rows <- rep(1:terra::nrow(layers), each = ncols) # 111 ... 222 ... 333
  cols <- rep(1:ncols, times = terra::nrow(layers)) # 123 ... 123 ... 123
  point.a <- matrix(rep(rows, length(rows)), nrow = numberOfCells, byrow = TRUE)
  point.b <- matrix(rep(cols, length(cols)), nrow = numberOfCells, byrow = TRUE)
  point.c <- Matrix::t(point.a)
  point.d <- Matrix::t(point.b)
  ## This appears to the the Balgovind form of the correlation function C, as
  ## mentioned by Ashok here:
  ## https://codereview.stackexchange.com/questions/224536. The above lines
  ## appear to be those shared by user "minem":
  ## https://codereview.stackexchange.com/a/224901.
  d <- sqrt((point.a - point.c)^2 + (point.b - point.d)^2)

  varCov.fun <-
    switch(variableCovarianceFunction,
           DBD = function() {
             Q.backgroundErrorStandardDeviation *
               Q.characteristicCorrelationLength^d
           },
           ## NOTE: isotropic Balgovind form of Q; the Balgovind model
           ## parameterizes the isotopic decaying correlation.
           Balgovind = function() {
             Q.backgroundErrorStandardDeviation *
               (1 + (d / Q.characteristicCorrelationLength)) *
               exp(-d / Q.characteristicCorrelationLength)
           },
           Exponential = function() {
             Q.backgroundErrorStandardDeviation *
               exp(-d / Q.characteristicCorrelationLength)
           },
           Guassian = function() {
             Q.backgroundErrorStandardDeviation *
               exp(-d^2 / 2 * Q.characteristicCorrelationLength^2)
           },
           Spherical = function() {
             ## NOTE: the Q.characteristicCorrelationLength actually refers to
             ## the radius of the spherical variance-covariance function.
             varMatrix <- Q.backgroundErrorStandardDeviation *
               ((3 * d) / (2 * Q.characteristicCorrelationLength) -
                d^3 / (2 * Q.characteristicCorrelationLength^3))
             varMatrix[d >= Q.characteristicCorrelationLength] <- 0
             varMatrix
           },
           stop(r"[Provided name of variableCovarianceFunction is invalid.
Valid function names are:
 - DBD
 - Balgovind
 - Exponential
 - Gaussian
 - Spherical]"))

  ## Assign to cells of the sparse matrix where the decay function has a value
  ## less than neighbourhood.
  Q[d < neighbourhood] <- varCov.fun()[d < neighbourhood]
  diag(Q) <- ifelse(Matrix::diag(Q) == 0, Q.backgroundErrorStandardDeviation, Matrix::diag(Q))

  if (compartmentsReported == 2) Q <- Matrix::bdiag(Q, Q)

  return(Q)
}

##' @description Run a SVEIRD compartmental model of an epidemic, optionally
##'   using Bayesian data assimilation.
##' @details TODO: DETAILS of the function.
##' @title SVEIRD compartmental model with optional Bayesian data assimilation
##' @param psi.diagonal TODO
##' @param layers The SpatRaster object with SVEIRD compartment layers, and a
##'   layer classifying habitation. Created with the getSVEIRD.SpatRaster
##'   function.
##' @param startDate The date (in YYYY-MM-DD format) the simulation begins.
##' @param countryCodeISO3C The ISO three character code for a recognized country.
##' @param lambda the average dialy movement distance of an individual (in
##'   kilometers).
##' @param aggregationFactor The number of adjacent cells in any one direction to
##'   aggregate into a single cell. The aggregation factor must be the same as
##'   that used to generate the SpatRaster for layers.
##' @param alpha The rate of vaccination (per day)
##' @param beta The rate of exposure (per day)
##' @param gamma The rate of becoming infectious (per day)
##' @param sigma The rate of recovery (per day)
##' @param delta The fatality rate (per day)
##' @param n.days The number of days the simulation will run for, beginning from
##'   the startDate.
##' @param neighbourhood.order The order of the queen's neighbourhood over which
##'   to distribute the Exposed and Infected intial state variable; zero is no
##'   neighbourhood, and corresponds to only the same cell where the queen
##'   already is. Higher orders follow the simple formula `(2 * x + 1)^2` to
##'   determine the number of cells over which to evenly disperse the exposed
##'   and infected state variables. All other state variables are placed
##'   directly in the grid cell corresponding to the health zone.
##' @param seedData a dataframe, as described in [initialInfections.fourCities].
##' @param simulationIsDeterministic Whether stochasticity is enabled or not; if
##'   the simulation is deterministic then no stochastic processes are used and
##'   the simulation is entirely deterministic. Either `TRUE` or `FALSE`.
##' @param dataAssimilationEnabled Whether Bayesian data assimilation will be
##'   used for state reporting data.
##' @param healthZoneCoordinates The coordinates of health zones in the country
##'   of interest, TODO: describe the use of the data so users understand why
##'   they must include it; e.g. (which will be used to group and summarize the
##'   compartmental model data at the end of the simmulation.)
##'
##'   \preformatted{
##'     HealthZone    Latitude    Longitude
##'     Alimbongo     -0.365515   29.1911818
##'     Beni          0.49113     29.47306
##'     Biena         0.57923     29.115633
##'     Butembo       0.140692    29.335014
##'     Goma          -1.658271   29.220132
##'     Kalunguta     0.323085    29.354677
##'     Katwa         0.116985    29.411838
##'     Kayna         -0.603936   29.174066
##'     Kyondo        -0.005622   29.408813
##'     Lubero        -0.15633    29.24057
##'     Mabalako      0.461257    29.210687
##'     Manguredjipa  0.353433    28.765479
##'     Masereka      -0.133333   29.333333
##'     Musienene     0.04022     29.26246
##'     Mutwanga      0.32893     29.74576
##'     Nyiragongo    -1.516667   29.25
##'     Oicha         0.698681    29.518834
##'     Pinga         -0.9830504  28.687911
##'     Vuhovi        0.1416      29.4075
##'     Ariwara       3.136873    30.706615
##'     Bunia         1.566667    30.25
##'     Komanda       1.367496    29.774322
##'     Lolwa         1.352969    29.498455
##'     Mambasa       1.359731    29.029226
##'     Mandima       1.35551     29.08173
##'     Nyakunde      1.431271    30.029283
##'     Rwampara      1.4053      30.3449
##'     Tchomia       1.4412      30.4845
##'   }
##' @param incidenceData A "situation report" dataframe. The first column
##'   provides the date of the officially reported, observed incidence of the
##'   disease, in ISO format (YYYY-MM-DD). MAYBE TODO: enforce the startDate
##'   parameter to be one week prior to the first observed data?
##'
##'   \preformatted{
##'     Date    Beni  Biena  Butembo  Goma  Kalunguta
##'     05-Aug  3     0      2        0     0
##'     12-Aug  2     0      0        0     0
##'     20-Aug  1     0      0        0     0
##'     26-Aug  5     0      0        0     0
##'     02-Sep  8     0      0        0     1
##'     09-Sep  5     0      2        0     0
##'     16-Sep  5     0      3        0     0
##'     23-Sep  4     0      1        0     0
##'     02-Oct  10    0      1        0     0
##'     07-Oct  14    0      4        0     0
##'     15-Oct  28    0      2        0     1
##'     21-Oct  19    0      3        0     0
##'     28-Oct  28    0      8        0     0
##'     04-Nov  16    0      6        0     1
##'     11-Nov  14    0      0        0     21
##'   }
##' @param deathData Data of the same format as incidenceData, but observations
##'   represent deaths, not infections.
##' @param variableCovarianceFunction Passed directly to [Q.forecastErrorCov()] to generate a
##'   Q matrix.
##' @param Q.backgroundErrorStandardDeviation TODO
##' @param Q.characteristicCorrelationLength TODO
##' @param neighbourhood TODO
##' @param callback a callback function to run, with no arguments, which will be
##'   called every time the main loop of the simulation iterates.
##' @returns a summaryTable dataframe for the simulation, showing changes in the
##'   compartment values over time, the daily values, and cumulative values.
##' @author Bryce Carson
##' @author Ashok Krishnmaurthy
##' @author Michael Myer
##' @export
##' @examples
##' subregionsSpatVector <- terra::vect(
##'   system.file(
##'     "extdata",
##'     ## COD: Nord-Kivu and Ituri (Democratic Republic of Congo)
##'     "subregionsSpatVector",
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' susceptibleSpatRaster <- terra::rast(
##'   system.file(
##'     "extdata",
##'     "susceptibleSpatRaster.tif", # Congo population
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' data("healthZonesCongo", package = "spatialEpisim.foundation")
##' data("initialInfections.fourCities", package = "spatialEpisim.foundation")
##' data("Congo.EbolaIncidence", package = "spatialEpisim.foundation")
##' SVEIRD.BayesianDataAssimilation(
##'   ## Parameters
##'   alpha = 3.5e-5,
##'   beta = 7e-3,
##'   gamma = 1/7,
##'   sigma = 1/36,
##'   delta = 2/36,
##'   lambda = 15,
##'   ## Model runtime
##'   n.days = 31, # a month, permitting three assimilations of observed data
##'   ## Model data
##'   seedData = initialInfections.fourCities,
##'   neighbourhood.order = 1,
##'   layers = getSVEIRD.SpatRaster(subregionsSpatVector,
##'                                 susceptibleSpatRaster,
##'                                 aggregationFactor = 35),
##'   aggregationFactor = 35,
##'   startDate = "2018-08-05",
##'   countryCodeISO3C = "COD",
##'   incidenceData = Congo.EbolaIncidence,
##'   ## Model options
##'   simulationIsDeterministic = TRUE,
##'   dataAssimilationEnabled = TRUE,
##'   healthZoneCoordinates = healthZonesCongo,
##'   variableCovarianceFunction = "DBD",
##'   ## Special parameters
##'   Q.backgroundErrorStandardDeviation = 0.55,
##'   Q.characteristicCorrelationLength = 6.75e-1,
##'   neighbourhood = 3,
##'   psi.diagonal = 1e-3
##' )
SVEIRD.BayesianDataAssimilation <-
  function(## Parameters influencing differential equations # TODO: verify the description of these parameters, i.e., do none of the other parameters influence the differential equations?
           alpha,
           beta,
           gamma,
           sigma,
           delta,

           ## Model runtime
           n.days,

           ## Model data
           seedData,        ## these three arguments influence the progression of infection
           neighbourhood.order = 0,  ## these three arguments influence the progression of infection
           lambda,          ## these three arguments influence the progression of infection
           layers,
           aggregationFactor,
           startDate,
           countryCodeISO3C,
           incidenceData = NULL,
           deathData = NULL,

           ## Model options
           simulationIsDeterministic = TRUE,
           dataAssimilationEnabled = FALSE,
           healthZoneCoordinates,
           variableCovarianceFunction,

           ## Special parameters
           Q.backgroundErrorStandardDeviation,
           Q.characteristicCorrelationLength,
           neighbourhood,
           psi.diagonal,

           ## Monitoring and logging
           callback = `{`) {
    ## MAYBE: is missing a better option? NULL defaults influence the decision...
    compartmentsReported <- sum(!is.null(incidenceData), !is.null(deathData))

    ## Preallocate a zeroed data frame with the following column names, and
    ## store it in a symbol named "summaryTable".
    names <- c(## Population and epidemic compartments (states)
               "N", "S", "V", "E", "I", "R", "D",
               ## Daily values of new vaccinations, exposures, infections,
               ## recoveries, and deaths
               "newV", "newE", "newI", "newR","newD",
               ## Cumulative values of exposed or infected people through the
               ## simulation runtime
               "cumE", "cumI")
    summaryTable <-
      data.frame(matrix(data = 0, ncol = length(names), nrow = n.days)) %>%
      "colnames<-"(names)

    ## NOTE: cast the seed data from the initial infections equitably, in a
    ## Moore Neighborhood of cells.
    seededLayers <- castSeedDataQueensNeighbourhood(seedData, neighbourhood.order, layers)

    adjustedSusceptible <- layers$Susceptible -
      seededLayers$Vaccinated -
      seededLayers$Exposed -
      seededLayers$Infected -
      seededLayers$Recovered -
      seededLayers$Dead

    layers$Vaccinated <- seededLayers$Vaccinated
    layers$Exposed <- seededLayers$Exposed
    layers$Infected <- seededLayers$Infected
    layers$Recovered <- seededLayers$Recovered
    layers$Dead <- seededLayers$Dead

    message(sprintf(r"{Adjusting the susceptible layer to account for seeded epidemic compartments...
Susceptible before seeding = %s
Susceptible after seeding = %s
}",
terra::global(layers$Susceptible, sum, na.rm = TRUE),
terra::global(adjustedSusceptible, sum, na.rm = TRUE)))
    stopifnot(adjustedSusceptible > 1)

    layers$Susceptible <- adjustedSusceptible

    if (dataAssimilationEnabled) {
      matrices.Bayes <-
        setupBayesianDataAssimilation(layers,
                                      healthZoneCoordinates,
                                      compartmentsReported,
                                      variableCovarianceFunction,
                                      Q.backgroundErrorStandardDeviation,
                                      Q.characteristicCorrelationLength,
                                      neighbourhood)
      H <- linearInterpolationMatrix <- matrices.Bayes$H
      HQHt <- matrices.Bayes$HQHt
      QHt <- matrices.Bayes$QHt
    }

    ## NOTE: preallocate the list which will hold a timeseries of SpatRaster
    ## objects.
    layers.timeseries <- vector(mode = "list", length = n.days)

    ## TODO: absolutely do not use this "datarow". There's a much better way.
    ## NOTE: datarow is a sentinel value used to prevent trying to assimilate
    ## more data than exists in the observed data dataframe. Seventy-six is a
    ## magic number, which is actually the number of rows of observed data plus
    ## one (for one-based indexing). The rows are weekly data; it may be easier
    ## to simply check if there is more data to assimilate, and whenever the
    ## simulation date is one day before or the same day as the observed data it
    ## is then assimilated. Using a sentinel value for the row from which we
    ## read the data to assimilate each week is unnecessary, though easy.
    datarow <- 1

    reclassifyNegatives <-
      function(spit) terra::classify(spit, cbind(-Inf, 1, 0), right = FALSE)

    ## TODO: all of the calcualtions within this loop should be spatial
    ## calcualtions on the SpatRasters, and no conversion to vector or matrix
    ## should be performed unless absolutely necessary.
    for (today in seq(n.days)) {
      ## TODO: At this URL, StackOverflow user Roman provides a reprex for a
      ## waitress callback function to generate a progress bar.
      ## https://stackoverflow.com/a/77496657/14211497. DONT change this; the
      ## callback function call here should only be modified to include a
      ## general set of arguments that a callback function may be interested in
      ## using. The arguments should be provided as a list.
      callback() # Run the callback function, or NULL expression.

      compartments <- terra::subset(layers, "Inhabited", negate = TRUE)
      compartments$Dead %<>% "*"(-1)
      compartments %<>% terra::global(sum, na.rm = TRUE)

      ## Set NSVEI counts in the summaryTable table; the date column is calculated later.
      summaryTable[today, 1] <- round(sum(compartments[-6, ], na.rm = TRUE)) # DONT include the dead.
      summaryTable[today, 2] <- round(compartments["Susceptible", ])
      summaryTable[today, 3] <- round(compartments["Vaccinated", ])
      summaryTable[today, 4] <- round(compartments["Exposed", ])
      summaryTable[today, 5] <- round(compartments["Infected", ])
      summaryTable[today, 6] <- round(compartments["Recovered", ])
      summaryTable[today, 7] <- round(compartments["Dead", ])

      ## The population is the sum of the susceptible, vaccinated, exposed,
      ## infected, and recovered compartments.
      numberLiving <- sum(terra::subset(layers, c("Dead", "Inhabited"), negate = TRUE))
      newVaccinated <- alpha * reclassifyNegatives(layers$Susceptible)

      ## Some susceptible people who come in contact with nearby infected are
      ## going to be newly exposed. Also known by the symbols βN⁻¹.
      proportionSusceptible <-
        terra::subst(layers$Susceptible / numberLiving, NaN, 0)

      ## NOTE: also known by the symbol Ĩ
      transmissionLikelihoods <-
        transmissionLikelihoodWeightings(layers$Infected, lambda, aggregationFactor)

      growth <- matrix(as.vector(beta * proportionSusceptible)
                       * as.vector(transmissionLikelihoods),
                       nrow = 20,
                       ncol = 14,
                       byrow = TRUE)

      ## TODO: stochasticity is not properly implemented yet; it was not fully
      ## supported in the previous implementation. It's likely that using
      ## stochasticity now will (still) produce an error.
      newExposed <- if(simulationIsDeterministic) growth else stats::rpois(1, growth)
      ## NOTE: any indices of these objects which were less than one will be the
      ## same indices that are set to zero in the newExposed object.
      indices <- c(proportionSusceptible < 1, transmissionLikelihoods < 1)
      newExposed[unlist(indices[-length(indices)])] <- 0
      dailyExposed <- sum(newExposed)

      newInfected <- reclassifyNegatives(gamma * (layers$Exposed + newExposed))
      dailyInfected <- sum(newInfected)

      ## Some infectious people are going to either recover or die
      infectious <- reclassifyNegatives(layers$Infected + newInfected)
      newRecovered <- sigma * infectious
      dailyRecovered <- sum(newRecovered)

      newDead <- reclassifyNegatives(delta * infectious)
      dailyDead <- sum(newDead)

      ## NOTE: each of these should be a SpatRaster, and thereby the arithmetic
      ## is cell-by-cell.
      layers$Susceptible <- layers$Susceptible - newExposed - newVaccinated;
      layers$Vaccinated  <- layers$Vaccinated  - newVaccinated
      layers$Exposed     <- layers$Exposed     + newExposed  - newInfected
      layers$Infected    <- layers$Infected    + newInfected - newDead - newRecovered
      layers$Recovered   <- layers$Recovered   + newRecovered
      layers$Dead        <- layers$Dead        + newDead

      ## NOTE: assimilate observed data weekly, not more or less frequently,
      ## while there is still data to assimilate. TODO: there are some
      ## alternatives that could be implemented for when to assimilate data; a
      ## check against the current simulaiton date and the reporting dates can
      ## be made, and whenever these align assimilation could occur.
      shouldAssimilateData <- all(dataAssimilationEnabled,
                                  today %% 7 == 0,
                                  datarow <= nrow(incidenceData))
      if (shouldAssimilateData) {
        infectedExposedLayers <- assimilateData(layers,
                                                linearInterpolationMatrix,
                                                incidenceData[datarow, -c(1, 2)],
                                                healthZoneCoordinates,
                                                psi.diagonal,
                                                QHt,
                                                HQHt)
        layers$Infected <- infectedExposedLayers$Infected
        layers$Exposed <- infectedExposedLayers$Exposed
        datarow <- datarow + 1
      }

      layers.timeseries[[today]] <- layers
    }

    ## TODO: the entire column can be calculated one time and added to the
    ## summaryTable data at the end of the function.
    ## summaryTable[today, 1]  <- toString(as.Date(startDate) + n.days(today - 1))
    ## NOTE: NA MAYBE a valid statistical value, it may not be appropriate to
    ## always replace it in our summaryTable table; for some variables it makes sense
    ## that the value is zero (no fatalities, infections, exposures, et cetera),
    ## but why would it be NA (missing) or NaN (not a number)?
    summaryTable[is.na(summaryTable)] <- 0
    return(list(table = tibble::as_tibble(summaryTable),
                rast = layers))
  }

##' Using the provided parameters and SpatRaster, the necessary setup functions
##' and values for Bayesian data assimilation are called, with values used later
##' on returned.
##' @title Setup Bayesian data assimilation
##' @param layers a SpatRaster with the following layers: Susceptible,
##'   Vaccinated, Exposed, Infected, Recovered, and Dead.
##' @param healthZoneCoordinates a dataframe with three columns, a location
##'   name, latitude, and longitude describing the geographical locations of
##'   reporting health zones.
##' @param compartmentsReported the number of compartments or epidemic states
##'   reported on; either one or two. Higher numbers are not supported. One
##'   corresponds only to infection, while two compartments corresponds to
##'   exposure and infection.
##' @param variableCovarianceFunction a function to calculate the error covariance, returned by Q.forecastErrorCov.
##' @param Q.backgroundErrorStandardDeviation the "background" or default amount of error, in standard deviations.
##' @param Q.characteristicCorrelationLength TODO
##' @param neighbourhood TODO
##' @returns the HQHt matrix
##' @author Bryce Carson
setupBayesianDataAssimilation <-
  function(layers,
           healthZoneCoordinates,
           compartmentsReported,
           variableCovarianceFunction,
           Q.backgroundErrorStandardDeviation,
           Q.characteristicCorrelationLength,
           neighbourhood) {
    ## Generate the linear interpolation operator matrix (function works for
    ## two compartments, at most).
    H <- linearInterpolationMatrix <-
      linearInterpolationOperator(layers,
                                  healthZoneCoordinates,
                                  compartmentsReported)

    ## NOTE: Q and H are both block diagonal, sparse matrices (but not of
    ## class sparseMatrix:
    ## <https://stat.ethz.ch/R-manual/R-patched/library/Matrix/html/sparseMatrix-class.html>).
    ## NOTE: create the model error covariance matrix, which, given we are
    ## using an ensemble-type data assimilation process, is time invariant.
    ## Immediately it is used to calculate QHt, and otherwise is unused.
    Q <- forecastErrorCovarianceMatrix <-
      Q.forecastErrorCov(layers,
                         variableCovarianceFunction,
                         Q.backgroundErrorStandardDeviation,
                         Q.characteristicCorrelationLength,
                         neighbourhood,
                         compartmentsReported)
    QHt <- Q %*% Matrix::t(H) # MAYBE TODO: alias these with better names.
    HQHt <- H %*% QHt # MAYBE TODO: alias these with better names.

    ## NOTE: this is based on old, dead code from the previous implementation,
    ## and also based on commented code from a StackOverflow question Ashok
    ## asked in July 2019: https://codereview.stackexchange.com/q/224536. It
    ## probably isn't necessary to retain, but it's here. Ashok can make a
    ## decision about its usage later. stopifnot(sum(eigen(Q)$values) ==
    ## ncell(layers))

    return(list(QHt = QHt, HQHt = HQHt, H = H))
  }

##' Assimilation of data (Bayesian) using optimal statistical inference, a
##' modified Kalman Filter.
##' @title Bayesian data assimilation
##' @param layers The SpatRaster object with layers Susceptible, Vaccinated,
##'   Exposed, Infected, Recovered, and Dead.
##' @param linearInterpolationMatrix a linear forward interpolation operator
##'   matrix with as many columns as cells in the SpatRaster layers, and as many
##'   rows as health zones for which coordinates are provided (when creating a
##'   matrix for use with one state vector). The matrix is either a trivial
##'   matrix as described elsewhere (see [linearInterpolationOperator]), or two
##'   partitions in a sparse, block diagonal matrix.
##' @param incidenceData A "situation report" dataframe. The first column
##'   provides the date of the officially reported, observed incidence of the
##'   disease, in ISO format (YYYY-MM-DD). MAYBE TODO: enforce the startDate
##'   parameter to be one week prior to the first observed data?
##'
##'   \preformatted{
##'     Date    Beni  Biena  Butembo  Goma  Kalunguta
##'     05-Aug  3     0      2        0     0
##'     12-Aug  2     0      0        0     0
##'     20-Aug  1     0      0        0     0
##'     26-Aug  5     0      0        0     0
##'     02-Sep  8     0      0        0     1
##'     09-Sep  5     0      2        0     0
##'     16-Sep  5     0      3        0     0
##'     23-Sep  4     0      1        0     0
##'     02-Oct  10    0      1        0     0
##'     07-Oct  14    0      4        0     0
##'     15-Oct  28    0      2        0     1
##'     21-Oct  19    0      3        0     0
##'     28-Oct  28    0      8        0     0
##'     04-Nov  16    0      6        0     1
##'     11-Nov  14    0      0        0     21
##'   }
##' @param healthZoneCoordinates a table of values giving the latitude and
##'   longitude coordinates for health zones in the country of interest. See the
##'   description of the healthZoneCoordinates argument in the function
##'   SVEIRD.BayesianDataAssimilation for more details.
##'
##'   \preformatted{
##'     HealthZone    Latitude    Longitude
##'     Alimbongo     -0.365515   29.1911818
##'     Beni          0.49113     29.47306
##'     Biena         0.57923     29.115633
##'     Butembo       0.140692    29.335014
##'     Goma          -1.658271   29.220132
##'     Kalunguta     0.323085    29.354677
##'     Katwa         0.116985    29.411838
##'     Kayna         -0.603936   29.174066
##'     Kyondo        -0.005622   29.408813
##'     Lubero        -0.15633    29.24057
##'     Mabalako      0.461257    29.210687
##'     Manguredjipa  0.353433    28.765479
##'     Masereka      -0.133333   29.333333
##'     Musienene     0.04022     29.26246
##'     Mutwanga      0.32893     29.74576
##'     Nyiragongo    -1.516667   29.25
##'     Oicha         0.698681    29.518834
##'     Pinga         -0.9830504  28.687911
##'     Vuhovi        0.1416      29.4075
##'     Ariwara       3.136873    30.706615
##'     Bunia         1.566667    30.25
##'     Komanda       1.367496    29.774322
##'     Lolwa         1.352969    29.498455
##'     Mambasa       1.359731    29.029226
##'     Mandima       1.35551     29.08173
##'     Nyakunde      1.431271    30.029283
##'     Rwampara      1.4053      30.3449
##'     Tchomia       1.4412      30.4845
##'   }
##' @param psi.diagonal TODO
##' @param QHt TODO
##' @param HQHt TODO
##' @returns a list of SpatRasters, the Infected and Exposed SpatRasters
##' @author Bryce Carson
assimilateData <-
  function(layers,
           linearInterpolationMatrix,
           incidenceData,
           healthZoneCoordinates,
           psi.diagonal,
           QHt,
           HQHt) {
    Infected <- terra::as.matrix(layers$Infected, wide = TRUE)
    ## TODO: Ashok is asking Thomas White about the meaning of "rat" in this
    ## context; he is checking with Thomas why 1e-9 was chosen as well.
    rat <- sum(terra::as.matrix(layers$Exposed, wide = TRUE)) / (sum(Infected) + 1e-9) # FIXME: no magic numbers, please.

    Prior <- matrix(Matrix::t(Infected), ncol = 1)

    HForecast <- linearInterpolationMatrix %*% Prior

    ## Create the measurement error covariance matrix.
    Innovation <- as.numeric(incidenceData) - HForecast

    Psi <- as.numeric(incidenceData)
    Psi[Psi < 1] <- psi.diagonal
    Psi <- diag(Psi)

    ## NOTE: The gain matrix, the Kalman filter, determines how the
    ## observational data are assimilated.
    KalmanFilter <- QHt %*% Matrix::solve(HQHt + Psi)
    Posterior <- Prior + KalmanFilter %*% Innovation
    Posterior[Posterior < 0] <- 0

    ## MAYBE TODO: is subsetting even necessary? Is Posterior larger than
    ## "seq(terra::nrow(layers) * terra::ncol(layers))", requiring us to
    ## subset it so that I is not too large? NOTE: when RESTACKING make sure
    ## byrow = TRUE. NOTE: what is restacking? Why did I choose this word
    ## when I wrote the first part of this comment a week ago? NOTE: take
    ## the subset of Posterior which is the same size as `layers`? Using single
    ## element subsetting assumes row-major ordering.
    I <- matrix(Posterior[seq(terra::nrow(layers) * terra::ncol(layers))],
                nrow = terra::nrow(layers),
                ncol = terra::ncol(layers),
                byrow = TRUE)

    ## NOTE: if an area is uninhabitable replace its value with zero; it makes
    ## more sense to instead use NA values to prevent calculating values for
    ## uninhabitable areas.
    infectious <- terra::mask(terra::"crs<-"(terra::"ext<-"(terra::rast(I), terra::ext(layers)), terra::crs(layers)),
                              layers$Inhabited,
                              maskvalues = 0,
                              updatevalue = 0)

    ## MAYBE FIXME: how, exaclty, does the number of compartments reported
    ## impact the assimilation of the data? What is the influence of the
    ## number of compartments reported on the overriding of Infected and
    ## Exposed, if only Infected data is observed? What if Infected and
    ## Exposed data are observed and assimilated, how is RAT used then?
    ## Should any changes be made in that case from the usual algorithm?
    exposures <- rat * layers$Infected

    return(list(Infected = infectious, Exposed = exposures))
  }

##' The Moore neighbourhood around the locations given in the seed data is
##' affected with the initial values recorded in the provided seed data. The
##' infected and exposed compartments are cast in a Moore neighbourhood of cells
##' about the geographical coordinates recorded locations, while other
##' compartments are cast directly into the cell corresponding exazctly to the
##' location, not a Moore neighbourhood.
##'
##' The Moore neighbourhood is calculated using a simple arithmetical algorithm.
##' @title Seed the initial infections, and other compartments, in a spatial
##'   simluation
##' @param seedData a dataframe like the following example; the compartment
##'   columns are the initial values.
##'
##' \preformatted{
##'   Location Latitude Longitude Vaccinated Exposed Infected Recovered Dead
##'   Beni     0.49113  29.47306  0          24      12       0         4
##'   Butembo  0.140692 29.335014 0          0       0        0         0
##'   Mabalako 0.461257 29.210687 0          0       0        0         0
##'   Mandima  1.35551  29.08173  0          0       0        0         0
##' }
##' @param neighbourhood.order The order of the queen's neighbourhood over which
##'   to distribute the Exposed and Infected intial state variable; zero is no
##'   neighbourhood, and corresponds to only the same cell where the queen
##'   already is. Higher orders follow the simple formula `(2 * x + 1)^2` to
##'   determine the number of cells over which to evenly disperse the exposed
##'   and infected state variables. All other state variables are placed
##'   directly in the grid cell corresponding to the health zone.
##' @param layers The SpatRaster object with layers Susceptible, Vaccinated,
##'   Exposed, Infected, Recovered, and Dead.
##' @author Bryce Carson
##' @returns the input layers SpatRaster, modified with the seeded data.
##' @export
##' @examples
##' subregionsSpatVector <- terra::vect(
##'   system.file(
##'     "extdata",
##'     ## COD: Nord-Kivu and Ituri (Democratic Republic of Congo)
##'     "subregionsSpatVector",
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' susceptibleSpatRaster <- terra::rast(
##'   system.file(
##'     "extdata",
##'     "susceptibleSpatRaster.tif", # Congo population
##'     package = "spatialEpisim.foundation",
##'     mustWork = TRUE
##'   )
##' )
##' layers <- getSVEIRD.SpatRaster(subregionsSpatVector,
##'                                susceptibleSpatRaster,
##'                                aggregationFactor = 35)
##' data(initialInfections.fourCities, package = "spatialEpisim.foundation")
##' plot(castSeedDataQueensNeighbourhood(initialInfections.fourCities, 0, layers))
##' castSeedDataQueensNeighbourhood(initialInfections.fourCities, 1, layers)
##' plot(castSeedDataQueensNeighbourhood(initialInfections.fourCities, 3, layers))
castSeedDataQueensNeighbourhood <- function(seedData, neighbourhood.order, layers) {
  stopifnot(neighbourhood.order %in% seq(floor(neighbourhood.order), ceiling(neighbourhood.order)))

  ## NOTE: evenly spread the count of exposed and infected persons across a
  ## Chess Queen's neighbourhood of a given order; this value is assigned to
  ## each cell of a neighbourhood of the same order in the SpatRaster, centered
  ## about the proper SpatRaster cell.
  neighbourhoodQuotient <- function(x) x / (2 * neighbourhood.order + 1)^2
  seedData.equitable <-
    dplyr::group_by(seedData, Location) %>%
    dplyr::summarize(dplyr::across(c("InitialExposed", "InitialInfections"),
                                   neighbourhoodQuotient)) %>%
    dplyr::right_join(seedData, by = dplyr::join_by(Location))

  for (seedingLocation in seedData.equitable$Location) {
    ## Get row and column numbers from the latitude and longitude for this
    ## health region.
    data <- dplyr::filter(seedData.equitable, Location == seedingLocation)
    row <- terra::rowFromY(layers, data$lat)
    col <- terra::colFromX(layers, data$lon)

    if (!any(is.na(c(row, col)))) {
      rowRange <- seq(from = row - neighbourhood.order, to = row + neighbourhood.order)
      columnRange <- seq(from = col - neighbourhood.order, to = col + neighbourhood.order)

      layers$Vaccinated[row, col]            <- data$InitialVaccinated
      layers$Exposed[rowRange, columnRange]  <- data$InitialExposed.x
      layers$Infected[rowRange, columnRange] <- data$InitialInfections.x
      layers$Recovered[row, col]             <- data$InitialRecovered
      layers$Dead[row, col]                  <- data$InitialDead
    }
  }

  ## FIXME: this results in layer$Susceptible having negative values, which is
  ## totally unrealistic!
  seeds <- terra::subset(layers, "Susceptible", negate = TRUE) * -1
  seededLayers <- c(terra::subset(layers, "Susceptible"), seeds)
  layers$Susceptible <- sum(seededLayers, na.rm = TRUE)

  return(layers)
}
