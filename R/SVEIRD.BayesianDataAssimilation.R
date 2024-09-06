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

##' Retrieve a SpatRaster useful for spatiotemporal epidemic compartmental
##' modelling.
##' @title Retrieve a population count SpatRaster
##' @param countryCodeISO3C The uppercase ISO three character code a recognized
##'   country.
##' @param folder Passed on to method downloadWorldPopData
##' @returns A SpatRaster of WorldPop population count data with the layer name
##'   "Population", with all NAs replaced by zeros.
##' @author Bryce Carson
##' @author Michael Myer
##' @author Ashok Krishnmaurthy
##' @export
##' @examples
##' getCountryPopulation.SpatRaster("COD")
getCountryPopulation.SpatRaster <- function(countryCodeISO3C, folder = NULL) {
  countrySpatRaster <-
    terra::rast(downloadWorldPopData(countryCodeISO3C,
                                     if (!is.null(folder)) folder))
  ## MAYBE TODO: replacing NAs with zeros has a performance penalty throughout
  ## the application, but ensures that the "algorithm" is closer to what it was
  ## before. It is worthwhile attempting to produce the same results without
  ## replacing NAs.
  replace(countrySpatRaster, is.na(countrySpatRaster), 0) %>%
    `names<-`("Population")
}

##' Acquire a SpatVector for a given country, optionally limited to subregion(s) thereof.
##'
##' @title Get a country's GADM SpatVector, from GADM verison 3.6
##' @param countryCodeISO3C The uppercase ISO three character code of a
##'   recognized country.
##' @param level1Region Subregions of the country to limit the extent of the
##'   SpatVector to.
##' @param folder The folder which should be searched for the GADM data, passed
##'   on to lvl1AdminBorders.
##' @returns SpatVector for the specificed country, optionally cropped to the
##'   provinces/states specified in the `level1Region` argument.
##' @author Bryce Carson
##' @export
##' @examples
##' library(geodata)
##' geodata_path_backup <- geodata_path()
##' options(geodata_default_path = getOption("spatialEpisim.foundation.datapath"))
##' getCountrySubregions.SpatVector("COD")
##' getCountrySubregions.SpatVector("COD", c("Nord-Kivu", "Ituri"))
##' getCountrySubregions.SpatVector("CZE", "Prague")
##' getCountrySubregions.SpatVector("NGA", "Kwara")
##' options(geodata_default_path = geodata_path_backup)
getCountrySubregions.SpatVector <- function(countryCodeISO3C = "COD",
                                            level1Region,
                                            folder) {
  stopifnot(countryCodeISO3C %in% countrycode::codelist$iso3c)

  if (missing(folder)) folder <- geodata::geodata_path()
  stopifnot(dir.exists(folder))

  countryVector <- geodata::gadm(countryCodeISO3C, path = folder, version = "3.6")

  if (!missing(level1Region))
    countryVector <- terra::subset(countryVector, countryVector$NAME_1 == level1Region)

  return(countryVector)
}

##' @title Create a SpatRaster of SVEIRD model compartment raster data
##' @description Create a SpatRaster objects, with a layer for each component in
##'   an SVEIRD epidemic model.
##' @details The SpatRaster objects for the VEIRD components are empty, the
##'   input SpatRaster is taken as the Suscptible layer, the only layer with
##'   non-zero values (if any existed before).
##' @param subregions a SpatVector object of subregions used to crop the
##'   population SpatRaster before creating the other layers in the raster
##'   object to represent the other compartments in a SVEIRD epidemic model. If
##'   it is missing no cropping is performed.
##' @param population a SpatRaster of population count data; it must be
##'   pre-aggregated if any aggregation is to be used during the simulation;
##'   none is performed by this function or downstream functions in the overall
##'   implementation of [SVEIRD.BayesianDataAssimilation] simulations.
##' @param aggregationFactor the number of cells in any direction to aggregate
##'   together into one, in the input raster, after masking with the subregions
##'   vector.
##' @returns a SpatRaster, with layers for the SVERID components
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
##'                      aggregationFactor = 10)
##'
##' ## Omitting the aggregation factor argument will prevent aggregation. An
##' ## aggregation factor of zero or one is meaningless and will produce an
##' ## error.
##' getSVEIRD.SpatRaster(subregionsSpatVector, susceptibleSpatRaster)
getSVEIRD.SpatRaster <- function(subregions, population, aggregationFactor = NULL) {
  terra::crs(subregions) <- terra::crs(population, proj = TRUE)
  population <- terra::crop(population, subregions, mask = TRUE)

  if (!is.null(aggregationFactor)) {
    population <- terra::aggregate(population, fact = aggregationFactor, fun = "sum", na.rm = TRUE)
  }

  ## NOTE: It is faster to use fun = NA, rather than fun = 0, because
  ## castSeedDataQueensNeighbourhood, for example, will not perform calculations,
  ## rather than multiply many cells by zero to no effect; it also produces
  ## clearer plots without values outside the bounds of the borders of the
  ## spatial region of the geographical region represented in the raster.
  zeroed <- terra::classify(population, cbind(0, Inf, 0))
  c(reclassifyBelowUpperBound(population, upper = 0), rep(zeroed, 5)) %>%
    "names<-"(c("Susceptible", "Vaccinated", "Exposed", "Infected", "Recovered",
                "Dead"))
}

##' The susceptible matrix is reclassified so positive values become one and
##' negative values become zero. Values which are zero remain zero. This results
##' in a binary classification of inhabitance when given a SpatRaster with
##' population count or susceptible data; regardless, whatever input SpatRaster
##' is provided, a binary classification of the values in the raster are
##' returned in a SpatRaster of equal dimension and with the same coordinate
##' reference system.
##'
##' The (re)-classification function built-in to terra is much faster than an
##' ad-hoc implementation of the same behaviour, of course, so the built-in
##' method is used (though it requries the reader to learn a little more of the
##' raster library; a good, a necessary, penalty).
##' @title Create a binary reclassification of the input SpatRaster
##' @param inputRaster A single-layer SpatRaster
##' @returns a binary reclassification of the input SpatRaster
##' @author Bryce Carson
classify.binary <- function(inputRaster) {
  stopifnot(terra::nlyr(inputRaster) == 1)
  ## Unit: milliseconds
  ##                            expr      min       lq     mean   median       uq     max neval
  ##  classify.binary(congo) 562.2161 781.7317 853.3972 802.9192 976.4928 1227.77   100
  ##
  ## terra::values(inputRaster)[terra::values(inputRaster) > 0] <- 1
  ## terra::values(inputRaster)[terra::values(inputRaster) < 1] <- 0

  ## Unit: milliseconds
  ##                            expr      min       lq     mean   median       uq      max neval
  ##  classify.binary(congo) 136.2108 147.4001 158.2525 154.9175 166.6182 236.5314   100
  terra::classify(inputRaster,
                  rbind(cbind(-Inf, 0,   0),
                        cbind(0,    Inf, 1)),
                  right = FALSE,
                  include.lowest = FALSE)
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
##' averageEuclideanDistance(lambda = 15, aggregationFactor = 10)
##'
##' lattice::levelplot(as.array(averageEuclideanDistance(15)))
##' lattice::levelplot(as.array(averageEuclideanDistance(100, 35)))
averageEuclideanDistance <-
  function(lambda, aggregationFactor = 1, epsilon = 0) {
    stopifnot(aggregationFactor > 0)
    radius <-
      if (lambda <= aggregationFactor)
        1
      else
        round((lambda - aggregationFactor) / aggregationFactor + epsilon) + 1

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

    ## NOTE: it is required that the matrix returned is square, and that its
    ## number of rows or columns are odd-numbered.
    stopifnot(dim(weights)[1] == dim(weights)[2])
    stopifnot(dim(weights)[1] %% 2 == 1)
    return(weights)
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
##'
##' The dimensions of the returned matrix are nrow(healthZoneCoordinates) rows
##'   by ncell(layers) columns.
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
##'                                 aggregationFactor = 10),
##'   healthZoneCoordinates = healthZonesCongo
##' )
##'
##' H.matrix <-
##'   linearInterpolationOperator(
##'     layers = getSVEIRD.SpatRaster(subregionsSpatVector,
##'                                   susceptibleSpatRaster,
##'                                   aggregationFactor = 10),
##'     healthZoneCoordinates = healthZonesCongo,
##'     neighbourhood.order = 1
##'   )
##' plot(as.array(H.matrix[1, ]))
##'
##' linearInterpolationOperator(
##'   layers = getSVEIRD.SpatRaster(subregionsSpatVector,
##'                                 susceptibleSpatRaster,
##'                                 aggregationFactor = 10),
##'   healthZoneCoordinates = healthZonesCongo,
##'   neighbourhood.order = 2
##' )
linearInterpolationOperator <- function(layers,
                                        healthZoneCoordinates,
                                        neighbourhood.order = 2,
                                        compartmentsReported = 1) {
  stopifnot(neighbourhood.order %in% c(0, 1, 2))
  if (neighbourhood.order == 2)
    stopifnot(terra::ncell(layers) >= 5 && terra::nrow(layers) >= 5) # required for 2nd order
  stopifnot(compartmentsReported %in% 1:2)

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
  ## to a matrix). Extract coordinates as cbind(lon, lat); it's stored as
  ## cbind(lat, lon).
  cells <- terra::cellFromXY(layers, terra::as.matrix(healthZoneCoordinates[, 3:2]))

  if (any(duplicated(cells))) {
    overaggregatedHealthZones <- tibble::tibble("Health Zone" = healthZoneCoordinates[, 1], Cell = cells) %>%
      dplyr::group_by(Cell) %>%
      dplyr::filter(dplyr::n() != 1)

    warning(sprintf("[Linear interpolation operator] Raster aggregation factor is too high to differentiate between two (or more) health zones (they correspond to the same grid cell).\n%s",
                    ## Based on the helpful answer by Richie Cotton on SO:
                    ## https://stackoverflow.com/a/26083626/14211497, which the following is
                    ## based on.
                    paste(capture.output(print(overaggregatedHealthZones)), collapse = "\n")))
  }
  if (any(is.na(cells)))
    warning("Ignoring NAs in [cells] object corresponding to coordinates out of bounds of [layers] SpatRaster.")

  cells <- cells[!is.na(cells)]

  ## NOTE: preallocate the linear forward interpolation matrix, with
  ## dimensions q * p (health zones by cells in the SpatRaster).
  H <- base::matrix(0, nrow(healthZoneCoordinates), terra::ncell(layers))

  ## NOTE: these are the weightings used for the chess queen zeroth, first,
  ## and second order neighbours. The zeroth order neighbor is the position of
  ## the queen itself.
  neighbour.weights <-
    switch(neighbourhood.order + 1, # the first of ... applies to zero, etc.
           1,
           c(2, 1) * 0.1,
           c(3, 2, 1) * 35^-1)

  for (index in seq_along(cells)) {
    H[index, cells[index]] <- neighbour.weights[1]

    if (neighbourhood.order != 0) {
      neighbour.1st <- queensNeighbours(1, cells[index], terra::ncol(layers))
      H[index, neighbour.1st] <- neighbour.weights[2]
    }

    if (neighbourhood.order == 2) {
      neighbour.2nd <- queensNeighbours(2, cells[index], terra::ncol(layers))
      if(anyDuplicated(c(neighbour.1st, neighbour.2nd)) > 0)
        simpleError("Duplicate cell indices among neighbours of multiple localities.")
      H[index, neighbour.2nd] <- neighbour.weights[3]
    }
  }

  stopifnot(dplyr::near(sum(H), nrow(healthZoneCoordinates)))
  stopifnot(dplyr::near(sum(matrix(H[1, ],
                                   ncol = ncol(layers),
                                   byrow = TRUE)),
                        1))

  dropExtendedArea <- function(row) {
    m <- matrix(row, byrow = TRUE, ncol = ncol(layers))
    m[(extend.length + 1):(base::nrow(m) - extend.length),
    (extend.length + 1):(base::ncol(m) - extend.length)] %>%
      Matrix::t() %>% # row-major order (byrow)
      as.vector()
  }
  H <- Matrix::t(apply(X = H, MARGIN = 1, FUN = dropExtendedArea))

  if (compartmentsReported == 2) H <- as.matrix(Matrix::bdiag(H, H))

  stopifnot(sum(.rowSums(H, m = nrow(H), n = ncol(H))) == (compartmentsReported * nrow(healthZoneCoordinates)))

  return(H)
}

##' @description generates a block diagonal error covariance matrix with exponential decay
##' @title Create a Forecast Error Covariance Matrix
##' @param layers The SpatRaster object with SVEIRD compartment layers, and a
##'   layer classifying habitation. Created with the getSVEIRD.SpatRaster
##'   function.
##' @param variableCovarianceFunction a covariance function used to determine
##'   the covariance between two variables.
##' @param forecastError.cov.sdBackground The "background" standard deviation of
##'   the covariances of the forecast covariance error matrix.
##' @param forecastError.cor.length "correlation length (i.e. the average size
##'   of the fluctuations)," as stated by J. Murray on the Physics
##'   StackExchange, <https://physics.stackexchange.com/a/671317>.
##' @param neighbourhood.Bayes TODO
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
##'                                aggregationFactor = 10)
##' Ituri.forecastError.cov <- forecastError.cov(layers,
##'                                              variableCovarianceFunction = "DBD",
##'                                              forecastError.cov.sdBackground = 2,
##'                                              forecastError.cor.length = 0.8,
##'                                              neighbourhood = 1,
##'                                              compartmentsReported = 2)
forecastError.cov <- function(layers,
                              variableCovarianceFunction,
                              forecastError.cov.sdBackground,
                              forecastError.cor.length,
                              neighbourhood.Bayes,
                              compartmentsReported = 2) {
  validFunctions <- c("DBD", "Balgovind", "Exponential", "Gaussian", "Spherical")
  if (!any(variableCovarianceFunction == validFunctions))
    stop(sprintf("%s is not a valid variable covariance function.", variableCovarianceFunction))

  columns <- terra::ncol(layers)
  rows <- terra::nrow(layers)
  block <- \(n) t(outer(seq(n), seq(n), "-"))
  ## I believe this is the called "the Balgovind form of the correlation
  ## [decay?] function C", as mentioned by Ashok here:
  ## https://codereview.stackexchange.com/questions/224536. Regardless of what
  ## it is called, this is a very effeciently implementation of it, with much
  ## less arithmetic than what was occurring before. It should be an order of
  ## magnitude faster, at least. MAYBE: "isotopic decaying correlation"?
  decay <- sqrt((matrix(1, columns, columns) %x% block(rows))^2 +
                (matrix(1, rows, rows) %x% block(columns))^2)

  ## Assign to cells of the sparse matrix where the decay function has a value
  ## less than neighbourhood.
  forecastErrorCovariance <- matrix(0, terra::ncell(layers), terra::ncell(layers))

  varCov.fun <-
    switch(variableCovarianceFunction,
           DBD = function() {
             forecastError.cov.sdBackground *
               forecastError.cor.length^decay
           },
           ## NOTE: isotropic Balgovind form of forecastErrorCovariance; the
           ## Balgovind model parameterizes the isotopic decaying correlation.
           Balgovind = function() {
             forecastError.cov.sdBackground *
               (1 + (decay / forecastError.cor.length)) *
               exp(-decay / forecastError.cor.length)
           },
           Exponential = function() {
             forecastError.cov.sdBackground *
               exp(-decay / forecastError.cor.length)
           },
           Gaussian = function() {
             forecastError.cov.sdBackground *
               exp(-decay^2 / 2 * forecastError.cor.length^2)
           },
           Spherical = function() {
             ## NOTE: the forecastError.cor.length actually refers to
             ## the radius of the spherical variance-covariance function.
             varMatrix <- forecastError.cov.sdBackground *
               ((3 * decay) / (2 * forecastError.cor.length) -
                decay^3 / (2 * forecastError.cor.length^3))
             varMatrix[decay >= forecastError.cor.length] <- 0
             varMatrix
           })

  forecastErrorCovariance[decay < neighbourhood.Bayes] <-
    varCov.fun()[decay < neighbourhood.Bayes]
  Matrix::diag(forecastErrorCovariance) <-
    ifelse(Matrix::diag(forecastErrorCovariance) == 0,
           forecastError.cov.sdBackground,
           Matrix::diag(forecastErrorCovariance))

  if (compartmentsReported == 2)
    forecastErrorCovariance <-
      Matrix::bdiag(forecastErrorCovariance, forecastErrorCovariance)

  return(forecastErrorCovariance)
}

##' @description Run a SVEIRD compartmental model of an epidemic, optionally using Bayesian data assimilation.
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
##'   directly in the grid cell corresponding to the health zone. At the moment
##'   orders higher than one are prohibited, so the only valid values for this
##'   argument are zero and one.
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
##'   disease, in YYYY-MM-DD format.
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
##' @param variableCovarianceFunction Passed directly to [forecastError.cov()] to generate a
##'   forecastErrorCovariance matrix.
##' @param forecastError.cov.sdBackground the "background" or default amount of error, in standard deviations.
##' @param forecastError.cor.length TODO
##' @param neighbourhood.Bayes The neighbourhood used in Bayesian data
##'   assimilation; this is different from that used in the casting of seed data
##'   about a neighbourhood of cells.
##' @param callback a callback function to run, with no arguments, which will be
##'   called every time the main loop of the simulation iterates, or a list of
##'   callback functions which run before, during, and after the loop, and which
##'   have the following structure. For each component, if the component is a
##'   list it must have members fun and args, where fun is a function symbol and
##'   args is a list of named arguments to the function; if it is not a list,
##'   the component of the list (before, during, or after) must be a function.
##'   See the examples.
##' @returns a list with components table, a tibble, and timeseries, a
##'   SpatRasterDataset representing timeseries of different variables.
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
##' rasterAggregationFactor = 9
##' simulation.days <- 31
##' if (requireNamespace("cli", quietly = TRUE)) {
##'   callback <- list(before = list(fun = cli::cli_progress_bar,
##'                                  args = list(name = "Simulating epidemic (SEI-type)",
##'                                              total = simulation.days)),
##'                    during = cli::cli_progress_update,
##'                    after = cli::cli_progress_done)
##' } else {
##'   callback <- "{" # does nothing
##' }
##' aggregatedCongoSubregions <-
##'   getSVEIRD.SpatRaster(subregionsSpatVector,
##'                        susceptibleSpatRaster,
##'                        aggregationFactor = rasterAggregationFactor)
##' SVEIRD.BayesianDataAssimilation(
##'   ## Parameters
##'   alpha = 3.5e-5,
##'   beta = 7e-3,
##'   gamma = 1/7,
##'   sigma = 1/36,
##'   delta = 2/36,
##'   lambda = 18,
##'   ## Model runtime
##'   n.days = simulation.days,
##'   ## Model data
##'   seedData = initialInfections.fourCities,
##'   neighbourhood.order = 1,
##'   layers = aggregatedCongoSubregions,
##'   aggregationFactor = rasterAggregationFactor,
##'   startDate = "2018-08-01",
##'   countryCodeISO3C = "COD",
##'   ## Model options
##'   simulationIsDeterministic = TRUE,
##'   dataAssimilationEnabled = FALSE,
##'   callback = callback
##' )
##'
##' ## A large, lengthy simulation with Bayesian data assimilation (runtime
##' ## approximately four minutes).
##' simulation.days <- 28
##' SVEIRD.BayesianDataAssimilation(
##'   ## Parameters
##'   alpha = 3.5e-5,
##'   beta = 7e-3,
##'   gamma = 1/7,
##'   sigma = 1/36,
##'   delta = 2/36,
##'   lambda = 18,
##'   ## Model runtime
##'   n.days = simulation.days,
##'   ## Model data
##'   seedData = initialInfections.fourCities,
##'   neighbourhood.order = 1,
##'   layers = aggregatedCongoSubregions,
##'   aggregationFactor = rasterAggregationFactor,
##'   startDate = "2018-08-01",
##'   countryCodeISO3C = "COD",
##'   incidenceData = head(Congo.EbolaIncidence, n = 4),
##'   ## Model options
##'   simulationIsDeterministic = TRUE,
##'   dataAssimilationEnabled = TRUE,
##'   healthZoneCoordinates = healthZonesCongo,
##'   variableCovarianceFunction = "DBD",
##'   ## Special parameters
##'   forecastError.cov.sdBackground = 0.55,
##'   forecastError.cor.length = 6.75e-1,
##'   neighbourhood.Bayes = 3,
##'   psi.diagonal = 1e-3
##' )
SVEIRD.BayesianDataAssimilation <-
  function(## Parameters influencing differential equations
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

           ## Bayesian data assimilation
           incidenceData,
           deathData,

           ## Model options
           simulationIsDeterministic = TRUE,
           dataAssimilationEnabled = FALSE,

           ## Bayesian data assimilation
           healthZoneCoordinates,
           variableCovarianceFunction,
           forecastError.cov.sdBackground,
           forecastError.cor.length,
           neighbourhood.Bayes,
           psi.diagonal,

           ## Monitoring and logging
           callback) {
    startDate <- lubridate::ymd(startDate)

    if (dataAssimilationEnabled) {
      ## TODO: the following two calls can be made into one call if some
      ## higher-order programming is used.
      if (!missing(incidenceData)) {
        if (!(startDate <= dplyr::first(incidenceData$Date <- lubridate::ymd(incidenceData$Date))))
          stop(sprintf("%s is not prior to or equal to %s.", startDate, dplyr::first(incidenceData$Date)))
      }
      if (!missing(deathData)) {
        if (!(startDate <= dplyr::first(deathData$Date <- lubridate::ymd(deathData$Date))))
          stop(sprintf("%s is not prior to or equal to %s.", startDate, dplyr::first(deathData$Date)))
      }
    } else {
      ## When data assimilation is not enabled, neither incidence nor death data
      ## should be provided.
      stopifnot(missing(incidenceData) && missing(deathData))
    }

    ## NOTE: Preallocate a zeroed data frame with the following column names, and
    ## store it in a symbol named "summaryTable".
    names <- c(
      ## Population and epidemic compartments (states)
      "N", # referred to as the count of "living" later on.
      "S", "V", "E", "I", "R", "D",
      ## Daily values of new vaccinations, exposures, infections,
      ## recoveries, and deaths
      "newV", "newE", "newI", "newR","newD",
      ## Cumulative values of exposed or infected people through the
      ## simulation runtime
      "cumE", "cumI"
    )
    summaryTable <-
      data.frame(matrix(data = 0, ncol = length(names), nrow = n.days)) %>%
      "colnames<-"(names)
    summaryTable <- cbind(tibble::tibble(Date = startDate + lubridate::days(0:(n.days - 1))),
                          summaryTable)

    layers %<>% castSeedDataQueensNeighbourhood(seedData, neighbourhood.order)
    timeseries <- terra::sds("names<-"(as.list(layers), names(layers)))

    if (dataAssimilationEnabled) {
      if (!missing(incidenceData) && !all(incidenceData$Date %in% summaryTable$Date)) {
        stop(sprintf("Some of the incidence data to be assimilated exists outside of the temporal bounds of the simulation.\nfirst(incidenceData$Date)\t%s\t\tfirst(summaryTable$Date)\t%s\nlast(incidenceData$Date)\t%s\t\tlast(summaryTable$Date)\t%s\n\nTry an earlier startDate argument, or increase the n.days argument.\n", dplyr::first(incidenceData$Date), dplyr::first(summaryTable$Date), dplyr::last(incidenceData$Date), dplyr::last(summaryTable$Date)))

        ## TODO: improve the error messaging to be clearer, offering an exact
        ## difference (in days) to the user; even better, adjust the startDate
        ## or n.days automatically with a warning.
        ## if (dplyr::first(incidenceData$Date) < dplyr::first(summaryTable$Date)) {
        ## } else if (dplyr::last(incidenceData$Date) < dplyr::last(summaryTable$Date)) {
        ## }
      }
      if (!missing(deathData) && !all(deathData$Date %in% summaryTable$Date)) {
        stop(sprintf("Some of the death data to be assimilated exists outside of the temporal bounds of the simulation.\nfirst(deathData$Date)\t%s\t\tfirst(summaryTable$Date)\t%s\nlast(deathData$Date)\t%s\t\tlast(summaryTable$Date)\t%s\n\nTry an earlier startDate argument, or increase the n.days argument.\n", dplyr::first(deathData$Date), dplyr::first(summaryTable$Date), dplyr::last(deathData$Date), dplyr::last(summaryTable$Date)))

        ## TODO: improve the error messaging to be clearer, offering an exact
        ## difference (in days) to the user; even better, adjust the startDate
        ## or n.days automatically with a warning.
        ## if (dplyr::first(deathData$Date) < dplyr::first(summaryTable$Date)) {
        ## } else if (dplyr::last(deathData$Date) < dplyr::last(summaryTable$Date)) {
        ## }
      }
      ## Create values needed for Bayesian data assimilation
      compartmentsReported <- sum(!missing(incidenceData), !missing(deathData))
      H <- linearInterpolationMatrix <-
        linearInterpolationOperator(layers,
                                    healthZoneCoordinates,
                                    compartmentsReported)
      ## NOTE: the model error covariance matrix is time-invariant.
      Q <- forecastErrorCovariance <-
        forecastError.cov(layers,
                          variableCovarianceFunction,
                          forecastError.cov.sdBackground,
                          forecastError.cor.length,
                          neighbourhood.Bayes,
                          compartmentsReported)
      QHt <- Matrix::tcrossprod(Q, H)
      HQHt <- H %*% QHt
    }

    ## NOTE: execute the "before" callback.
    if (!missing(callback) && hasName(callback, "before")) {
      if (all(c("fun", "args") %in% names(callback$before))) {
        if (is.function(callback$before$fun) && is.list(callback$before$args))
          do.call(callback$before$fun, args = callback$before$args)
        else
          simpleError(r"[The "before" callback is malformed.]")
      } else if (is.function(callback$before)) {
        callback$before()
      }
    }

    for (today in seq(n.days)) {
      ## NOTE: execute the "during" callback.
      if (!missing(callback)) {
        if (hasName(callback, "during")) {
          if (all(c("fun", "args") %in% names(callback$during))) {
            if (is.function(callback$during$fun) && is.list(callback$during$args))
              do.call(callback$during$fun, args = callback$during$args)
            else
              simpleError("The `during' callback is malformed.")
          } else if (is.function(callback$during)) {
            callback$during()
          }
        } else if (is.function(callback)) {
          callback()
        }
      }

      living <- terra::subset(layers, "Dead", negate = TRUE) %>%
        terra::app("sum", na.rm = TRUE) %>%
        "names<-"("Living")

      ## NOTE: set the previous timesteps compartment count values in the
      ## summary table before calculating values for the current timestep.
      summaryTable[today, "N"] <- round(terra::global(living, "sum", na.rm = TRUE))
      summaryTable[today, c("S", "V", "E", "I", "R", "D")] <-
        round(sums <- t(terra::global(layers, "sum", na.rm = TRUE)))

      newVaccinated <- alpha * reclassifyBelowUpperBound(layers$Susceptible, upper = 1)

      proportionSusceptible <- layers$Susceptible / living # alike total-mass action in Episim.

      ## NOTE: Calculate a matrix of weights respecting human mobility patterns.
      transmissionLikelihoods <-
        terra::focal(x = layers$Infected,
                     w = averageEuclideanDistance(lambda, aggregationFactor),
                     na.policy = "omit",
                     fillvalue = 0,
                     fun = "sum",
                     na.rm = TRUE)

      uniqueInfectionLikelihoods <- length(unique(as.vector(transmissionLikelihoods)))
      if (is.nan(sums[, "Infected"]))
        warning("The result of (global) sum total of the Infected compartment is NaN according to terra, so the check that the number of unique infection likelihoods is greater than one is being skipped this iteration. See issue #13.")
      if (all(!is.nan(sums[, "Infected"]),
              as.numeric(sums[, "Infected"]) >= 0,
              !(uniqueInfectionLikelihoods > 1)))
        stop("The number of unique likelihoods of transmission is not more than one, indicating an issue generating the transmissionLikelihoods matrix.")

      growth <- beta * proportionSusceptible * transmissionLikelihoods
      ## MAYBE FIXME: there should not be two layers returned; the masking is not as intended.
      newExposed <- terra::mask(if(simulationIsDeterministic) growth else stats::rpois(1, growth),
                                c(proportionSusceptible < 1, transmissionLikelihoods < 1),
                                maskvalues = TRUE,
                                updatevalue = 0)

      newInfected <- reclassifyBelowUpperBound(gamma * sum(c(layers$Exposed, newExposed)), upper = 1)

      ## Some infectious people are going to either recover or die
      infectious <- reclassifyBelowUpperBound(sum(c(layers$Infected, newInfected)), upper = 1)
      newRecovered <- sigma * infectious

      newDead <- reclassifyBelowUpperBound(delta * infectious, upper = 1)

      layers$Susceptible <- sum(c(layers$Susceptible, -1 * newExposed, -1 * newVaccinated), na.rm = TRUE)
      layers$Vaccinated  <- sum(c(layers$Vaccinated, newVaccinated), na.rm = TRUE)
      layers$Exposed     <- sum(c(layers$Exposed, newExposed, -1 * newInfected), na.rm = TRUE)
      layers$Infected    <- sum(c(layers$Infected, newInfected, -1 * newDead, -1 * newRecovered), na.rm = TRUE)
      layers$Recovered   <- sum(c(layers$Recovered, newRecovered), na.rm = TRUE)
      layers$Dead        <- sum(c(layers$Dead, newDead), na.rm = TRUE)


      if (dataAssimilationEnabled && summaryTable[today, "Date"] %in% incidenceData$Date) {
        todaysIncidenceData <- (dplyr::filter(incidenceData, Date == summaryTable[today, "Date"]))[, -c(1, 2)]
        infectedExposedLayers <-
          assimilateData(layers,
                         linearInterpolationMatrix,
                         todaysIncidenceData,
                         healthZoneCoordinates,
                         psi.diagonal,
                         QHt,
                         HQHt)
        layers$Infected <- infectedExposedLayers$Infected
        layers$Exposed <- infectedExposedLayers$Exposed
      }

      for(layer in names(layers)) {
        i <- which(layer == names(layers))
        terra::add(timeseries[i]) <- layers[layer]
      }
    }

    ## NOTE: execute the "after" callback.
    if (!missing(callback) && hasName(callback, "after")) {
      if (all(c("fun", "args") %in% names(callback$after))) {
        if (is.function(callback$after$fun) && is.list(callback$after$args))
          do.call(callback$after$fun, args = callback$after$args)
        else
          simpleError("The \"after\" callback is malformed.")
      } else if (is.function(callback$after)) {
        callback$after()
      }
    }

    summaryTable %<>%
      dplyr::mutate(newV = dplyr::lead(V) - V,
                    newE = dplyr::lead(E) - E,
                    newI = dplyr::lead(I) - I,
                    newR = dplyr::lead(R) - R,
                    newD = dplyr::lead(D) - D,
                    cumE = dplyr::first(E) + cumsum(dplyr::case_when(newE > 0 ~ newE, .default = 0)),
                    cumI = dplyr::first(I) + cumsum(dplyr::case_when(newI > 0 ~ newI, .default = 0)))

    stopifnot(unique(terra::nlyr(timeseries)) == n.days + 1)
    terra::time(timeseries) <- lubridate::date(startDate) +
      seq(from = 0, length.out = unique(terra::nlyr(timeseries)))

    return(list(table = tibble::as_tibble(summaryTable),
                timeseries = timeseries))
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
##' @param prevalenceData A "situation report" dataframe. The first column
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
##' @param psi.diagonal A replacement value for elements of the Psi matrix'
##'   diagonal which are zero.
##' @param QHt the [tcrossprod] of the Q and H matrices.
##' @param HQHt the product of the H matrix and QHt.
##' @returns a SpatRaster with Infected and Exposed compartment layers.
##' @author Bryce Carson
assimilateData <-
  function(layers,
           linearInterpolationMatrix,
           prevalenceData,
           healthZoneCoordinates,
           psi.diagonal,
           QHt,
           HQHt) {
    Infected <- terra::as.matrix(layers$Infected, wide = TRUE)

    Prior <- matrix(Matrix::t(Infected), ncol = 1)
    ## FIXME DONE: all elements in the Forecast matrix are NaN. TODO:
    ## investigate why the Prior contains NaNs.
    if (any(is.nan(Prior))) {
      warning("Prior contains NaNs, replacing with zeroes.")
      Prior[is.nan(Prior)] <- 0
    }
    Forecast <- linearInterpolationMatrix %*% Prior
    ## Create the measurement error covariance matrix.
    Innovation <- as.numeric(prevalenceData) - Forecast
    Innovation <- as.numeric(prevalenceData) - (linearInterpolationMatrix %*% Prior)

    Psi <- as.numeric(prevalenceData)
    Psi[Psi == 0] <- psi.diagonal
    Psi <- diag(Psi)

    ## NOTE: I did not have great documentation left for myself regarding this.
    ## Why was I concerned about the dimensions? Why was I concerend whether it
    ## was a direct sum or Kronecker sum (or not)?
    if (all(!is.null(getOption("spatialEpisim.foundation.debugMessages")),
            getOption("spatialEpisim.foundation.debugMessages"),
            !identical(dim(HQHt), dim(Psi))))
      message("HQHT + diag(Psi) is either a direct sum or a Kronecker sum, because the dimensions are inequal.")

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
    if (!identical(dim(Posterior), dim(layers)))
      message("Posterior and `layers` don't have the same dimensions, so subsetting Posterior is indeed required.")

    updatedSpatRaster <-
      matrix(Posterior[seq(terra::nrow(layers) * terra::ncol(layers))],
             nrow = terra::nrow(layers),
             ncol = terra::ncol(layers),
             byrow = TRUE) %>%
      terra::rast() %>%
      terra::"ext<-"(terra::ext(layers)) %>%
      terra::"crs<-"(terra::crs(layers)) %>%
      terra::mask(classify.binary(layers$Susceptible), # Inhabited
                  maskvalues = 0,
                  updatevalue = 0) %>%
      "names<-"("Infected")

    if (all(is.nan(unique(terra::values(updatedSpatRaster))))) {
      eachElementOfInfectedIsNaN <- "All values in the Infected compartment update are NaN; that's incorrect!"
      if (interactive()) {
        warning(eachElementOfInfectedIsNaN)
        browser()
      } else  {
        stop(eachElementOfInfectedIsNaN)
      }
    }

    ## NOTE: reclassifying NAs as zeroes will cause the geographical appearance
    ## of the plotted raster to be lost, because values outside of the actual
    ## geographical bounds (but within the bounds of the SpatRaster) will be set
    ## to zero, so the raster will be more of a "field" or "plane" than a raster
    ## of a geographic area with raster values for some variable.
    ## rcl <- rbind(cbind(NA, 0), cbind(NaN, 0))

    ## NOTE: reclassifying NaNs is the proper way to handle division by zero.
    rcl <- cbind(NaN, 0)

    ## DONE: in previous commits there was a question here; my answer to that
    ## question is, "Exposed data should be assimilated before preserving the
    ## ratio of the forecast and the analysis of the exposed and infected
    ## compartments; the forecast of exposures should be replaced with the
    ## observed data prior to preserving the ratio of exposed to infected".
    return("names<-"(c(updatedSpatRaster,
                       terra::classify((layers$Exposed * updatedSpatRaster) / layers$Infected,
                                       rcl)),
                     c("Infected", "Exposed")))
  }

##' The Moore neighbourhood around the locations given in the seed data is
##' affected with the initial values recorded in the provided seed data. The
##' infected and exposed compartments are cast in a Moore neighbourhood of cells
##' about the geographical coordinates recorded locations, while other
##' compartments are cast directly into the cell corresponding exactly to the
##' location, not a Moore neighbourhood.
##'
##' The Moore neighbourhood is calculated using a simple arithmetical algorithm.
##' @title Seed the initial infections, and other compartments, in a spatial
##'   simulation
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
##'                                aggregationFactor = 10)
##' data(initialInfections.fourCities, package = "spatialEpisim.foundation")
##' terra::plot(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 0)$Susceptible)
##' terra::plot(castSeedDataQueensNeighbourhood(layers, initialInfections.fourCities, 1)$Susceptible)
castSeedDataQueensNeighbourhood <-
  function(layers, seedData, neighbourhood.order) {
    ## This function works for higher orders, but we limit it to zeroth or first order.
    stopifnot(any(neighbourhood.order == c(0, 1)))
    sumBeforeSeeding <- sum(terra::global(layers, "sum", na.rm = TRUE))

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
        rowRange <- seq(from = row - neighbourhood.order,
                        to = row + neighbourhood.order)
        columnRange <- seq(from = col - neighbourhood.order,
                           to = col + neighbourhood.order)

        layers$Vaccinated[row, col]            <- data$InitialVaccinated
        layers$Exposed[rowRange, columnRange]  <- data$InitialExposed.x
        layers$Infected[rowRange, columnRange] <- data$InitialInfections.x
        layers$Recovered[row, col]             <- data$InitialRecovered
        layers$Dead[row, col]                  <- data$InitialDead
      }
    }

    layers$Susceptible <- layers$Susceptible -
      layers$Vaccinated -
      layers$Exposed -
      layers$Infected -
      layers$Recovered -
      layers$Dead

    stopifnot(terra::global(layers$Susceptible, min, na.rm = TRUE) >= 0)
    stopifnot(sum(terra::global(layers, "sum", na.rm = TRUE)) == sumBeforeSeeding)
    return(layers)
  }

##' Change all numbers in a raster which are less than some upper value to zero.
##'
##' With an upper of zero, negatives will be reclassified to zero, whilst with
##' an upper bound of some positive number, any number less than that will be
##' reclassified to zero.
##' @title Reclassify Below Upper Bound
##' @param raster a SpatRaster to reclassify
##' @param upper the upper bound, not included in the range of values classified
##' @returns the reclassified SpatRaster
##' @author Bryce Carson
reclassifyBelowUpperBound <- function(raster, upper) {
  terra::classify(raster, cbind(-Inf, upper, 0), right = FALSE, include.lowest = FALSE)
}
