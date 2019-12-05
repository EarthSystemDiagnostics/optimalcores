##
## aim:
## functions to select ECHAM5-wiso data and data from a specific target site
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores
##

##' Select climate model data
##'
##' Select the data from a certain climate model simulation.
##' @param data length one character vector with the short name of the climate
##' model simulation to load from the data folder. Defaults to \code{'past1000'}
##' for the currently available ECHAM5-MPIOM-wiso past millennium simulation.
##' @return A list of 10 elements; each list element is a two-dimensional array
##' (time, position) with the time series for the respective model variable. The
##' variables are:
##' * \code{t2m}: annual mean 2 metre temperature
##' * \code{t2m.pw}: annual mean 2 metre temperature precipitation-weighted
##' * \code{oxy}: annual mean d18O of precipitation
##' * \code{oxy.pw}: annual mean d18O of precipitation precipitation-weighted
##' * \code{prc}: annual precipitation
##' * \code{lnd.t2m}: as above but only for continental Antarctic sites
##' * \code{lnd.t2m.pw}: as above but only for continental Antarctic sites
##' * \code{lnd.oxy}: as above but only for continental Antarctic sites
##' * \code{lnd.oxy.pw}: as above but only for continental Antarctic sites
##' * \code{lnd.prc}: as above but only for continental Antarctic sites.
##' The first five list elements are \code{"pField"} objects, the latter five
##' are \code{"pTs"} objects.
##' @seealso \code{\link[ecustools]{pField}}, \code{\link[ecustools]{pTs}}
##' @author Thomas Münch
selectData <- function(data = "past1000") {

  # Load data set

  if (data == "past1000") {
    filePath <- file.path("data",
                         "echam5-mpiom-wiso-past1000_ts-fields-annual.rda")
  } else if (data == "eraint") {
    filePath <- file.path("data",
                         "echam5-wiso-era_ts-fields-annual.rda")
  }  else if (data == "piControl") {
    filePath <- file.path("data",
                         "echam6-mpiom-wiso-piControl_ts-fields-annual.rda")
  } else {
    stop("Unknown data set requested.", call. = FALSE)
  }

  if (file.exists(filePath)) {
    load(filePath)
  } else {
    stop(sprintf("No '%s' data found in data/.", data), call. = FALSE)
  }
  
  # Output data

  res <- list()

  res$t2m    <- tsfld.ann$t2m
  res$t2m.pw <- tsfld.ann$t2m.pw
  res$oxy    <- tsfld.ann$oxy
  res$oxy.pw <- tsfld.ann$oxy.pw
  res$prc    <- tsfld.ann$prc

  # filter out sites only containing NA (ocean sites)
  i.sel <- which(!is.na(colSums(tsfld.ann$t2m)))

  res$lnd.t2m    <- res$t2m[, i.sel]
  res$lnd.t2m.pw <- res$t2m.pw[, i.sel]
  res$lnd.oxy    <- res$oxy[, i.sel]
  res$lnd.oxy.pw <- res$oxy.pw[, i.sel]
  res$lnd.prc    <- tsfld.ann$prc[, i.sel]

  return(res)

}

##' Select a target site
##'
##' Select a target site (i.e. a certain grid cell) from a spatial field based
##' on latidude and longitude coordinates.
##' @param field a \code{\link[ecustools]{pField}} object with the spatial field
##' of a certain variable.
##' @param site a length-one character vector giving a shorthand ID for the
##' requested target site; see details for a list of implemented sites; defaults
##' to the EDML site at Kohnen Station. Alternatively, set \code{site} to
##' \code{NULL} and directly specify target latitude and longitude coordinates.
##' @param lat0 requested latitude of the target site (only used if
##' \code{site = NULL}).
##' @param lon0 requested longitude of the target site (only used if
##' \code{site = NULL}).
##' @details
##' Implemented shorthand IDs for target sites are:
##' * "edml":   EDML deep drilling site at Kohnen Station
##' * "vostok": Vostok station
##' * "domef":  Dome F deep drilling site
##' * "domec":  Dome C deep drilling site
##' * "talos":  Talos Dome deep drilling site
##' * "wdc":    WAIS Divide deep drilling site.
##'
##' The actual target site is set to the nearest grid cell of the spatial
##' \code{field} based on the distance minimization approach in
##' \code{\link[ecustools]{SelPoint}}.
##'
##' Note that the output field of distances relative to the target site only
##' refers to continental Antarctic sites.
##' @return A named list of the following four elements:
##' * "dat":  the time series at the specified target site.
##' * "lat0": latitude of actual grid cell used closest to the requested site.
##' * "lin0": longitude of actual grid cell used closest to the requested site.
##' * "dist": 2D array with distances relative to the target site of the grid
##'           cells in \code{field}.
##' @author Thomas Münch
setTarget <- function(field, site = "edml", lat0 = NULL, lon0 = NULL) {

  # Set target site

  if (!is.null(site)) {
    
    if (site == "edml") {
      lat0 <- -75
      lon0 <- 0
    } else if (site == "vostok") {
      lat0 <- -78.47
      lon0 <- 106.83
    } else if (site == "domef") {
      lat0 <- -77.32
      lon0 <- 39.7
    } else if (site == "domec") {
      lat0 <- -74.7
      lon0 <- 124.2
    } else if (site == "talos") {
      lat0 <- -72.8
      lon0 <- 159.1
    } else if (site == "wdc") {
      lat0 <- -79.46
      lon0 <- 360 - 112.09
    } else {
      stop("Unknown target site.")
    }

  } else {

    if (length(lat0) == 0 & length(lon0) == 0)
      stop("Specifiy both 'lat0' and 'lon0'.")

  }

  # Select target and output data

  res <- list()
  
  res$dat <- ecustools::SelPoint(field, lat = lat0, lon = lon0)

  res$lat0 <- ecustools::GetLat(res$dat)
  res$lon0 <- ecustools::GetLon(res$dat)

  if (res$lon0 >= 180) res$lon0 <- res$lon0 - 360

  # distances relative to target site of only land sites
  i.sel <- which(apply(field, 2, function(c){!all(is.na(c))}))

  # check whether subset field is pField or pTs
  if (length(i.sel == ncol(field))) {
    # is pField
    res$dist <-
      ecustools::GetDistanceField(field,
                                  lat = lat0, lon = lon0,
                                  get.nearest = TRUE)
    
  } else {
    # is pTs
    res$dist <-
      ecustools::GetDistance(lat0 = lat0, lon0 = lon0,
                             lat = ecustools::GetLat(field[, i.sel]),
                             lon = ecustools::GetLon(field[, i.sel]),
                             get.nearest = TRUE)

  }

  return(res)

}

