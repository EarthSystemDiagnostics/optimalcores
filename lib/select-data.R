##
## aim:
## functions to select ECHAM5-wiso data and data from a specific target site
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores
##

##' Select climate model data
##'
##' Select the data from a certain climate model simulation.
##'
##' @param data length one character vector with the short name of the climate
##'   model simulation to load from the data folder. Defaults to
##'   \code{'past1000'} for the currently available ECHAM5-MPIOM-wiso past
##'   millennium simulation.
##' @param threshold factor to define an outlier threshold applied to the time
##'   series from each model grid cell; i.e., time series values which lie
##'   above, or below, the time series' interquartile range times this factor
##'   are defined as outliers and are set to \code{NA}. Note that the outlier
##'   removal is only performed on the oxygen isotope and precipitation-weighted
##'   oxygen isotope fields.
##' @return A list of 10 elements; each list element is a two-dimensional array
##'   (time, position) with the time series for the respective model
##'   variable. The variables are:
##'   * \code{t2m}: annual mean 2 metre temperature
##'   * \code{t2m.pw}: annual mean 2 metre temperature precipitation-weighted
##'   * \code{oxy}: annual mean d18O of precipitation
##'   * \code{oxy.pw}: annual mean d18O of precipitation precipitation-weighted
##'   * \code{prc}: annual precipitation
##'   * \code{lnd.t2m}: as above but only for continental Antarctic sites
##'   * \code{lnd.t2m.pw}: as above but only for continental Antarctic sites
##'   * \code{lnd.oxy}: as above but only for continental Antarctic sites
##'   * \code{lnd.oxy.pw}: as above but only for continental Antarctic sites
##'   * \code{lnd.prc}: as above but only for continental Antarctic sites.
##'   The first five list elements are \code{"pField"} objects, the latter five
##'   are \code{"pTs"} objects.
##' @seealso \code{\link[pfields]{pField}}, \code{\link[pfields]{pTs}}
##' @author Thomas Münch
selectData <- function(data = "past1000", threshold = 4) {

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
  res$oxy    <- ApplyTime(tsfld.ann$oxy, removeOutlier,
                          threshold = threshold)
  res$oxy.pw <- ApplyTime(tsfld.ann$oxy.pw, removeOutlier,
                          threshold = threshold)
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

##' Locate outliers
##'
##' Find the positions in a time series at which the values exceed a certain
##' threshold.
##'
##' @param x a time series to analyse.
##' @param threshold factor to define the outlier threshold; i.e., values which
##'   lie above, or below, the data's interquartile range times this factor are
##'   defined as outliers.
##' @return the index positions in \code{x} of the outliers or \code{integer(0)}
##'   if no outliers are found.
##' @author Thomas Münch
locateOutlier <- function(x, threshold) {

  lower.whisker <- boxplot(x, range = threshold, plot = FALSE)$stats[1, 1]
  upper.whisker <- boxplot(x, range = threshold, plot = FALSE)$stats[5, 1]

  c(which(x < lower.whisker), which(x > upper.whisker))

}

##' Remove outliers
##'
##' Find the positions in a time series at which the values exceed a certain
##' threshold and set these values to \code{NA}.
##'
##' @param x a time series to analyse.
##' @param threshold factor to define the outlier threshold; i.e., values which
##'   lie above, or below, the data's interquartile range times this factor are
##'   defined as outliers.
##' @return the input time series with outliers removed.
##' @author Thomas Münch
removeOutlier <- function(x, threshold) {

  if (length(i <- locateOutlier(x, threshold))) x[i] <- NA
  return(x)
}


##' Select a target site
##'
##' Select a target site (i.e. a certain grid cell) from a spatial field based
##' on latidude and longitude coordinates.
##' @param field a \code{\link[pfields]{pField}} object with the spatial field
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
##' \code{\link[pfields]{SelPoint}}.
##'
##' Note that the output field of distances relative to the target site only
##' refers to continental Antarctic sites.
##' @return A named list of the following four elements:
##' * "dat":  the time series at the specified target site.
##' * "lat0": latitude of actual grid cell used closest to the requested site.
##' * "lon0": longitude of actual grid cell used closest to the requested site.
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
  
  res$dat <- pfields::SelPoint(field, lat = lat0, lon = lon0)

  res$lat0 <- pfields::GetLat(res$dat)
  res$lon0 <- pfields::GetLon(res$dat)

  if (res$lon0 >= 180) res$lon0 <- res$lon0 - 360

  # distances relative to target site of only land sites
  i.sel <- which(apply(field, 2, function(c){!all(is.na(c))}))

  # check whether subset field is pField or pTs
  if (length(i.sel) == ncol(field)) {
    # is pField
    res$dist <-
      pfields::GetDistanceField(field,
                                lat = lat0, lon = lon0,
                                get.nearest = TRUE)
    
  } else {
    # is pTs
    res$dist <-
      geostools::GetDistance(lat0 = lat0, lon0 = lon0,
                             lat = pfields::GetLat(field[, i.sel]),
                             lon = pfields::GetLon(field[, i.sel]),
                             get.nearest = TRUE)

  }

  return(res)

}

##' Set target region
##'
##' Select those grid cells from a \code{"pfield"} or \code{"pTs"} object that
##' lie within a defined latitude-longitude region.
##' @param field a \code{"pfield"} or \code{"pTs"} object; see
##' \code{\link[pfields]{pField}} and \code{\link[pfields]{pTs}}.
##' @param min.lat the minimum latitude of the region.
##' @param max.lat the maximum latitude of the region.
##' @param min.lon the minimum longitude of the region.
##' @param max.lon the maximum longitude of the region.
##' @param verbose logical; if \code{TRUE} (default), prints a message giving
##'   the distance in km between the midpoint of the defined region and its
##'   outer corners.
##' @return a data frame with three columns giving the indices
##'   (\code{field.indices}) of the grid cells of \code{field} that lie within
##'   the defined region and their latitude (\code{lat}) and longitudes
##'   (\code{lon}).
##' @author Thomas Münch
setTargetRegion <- function(field,
                            min.lat = -80, max.lat = -70,
                            min.lon = -17.5, max.lon = 17.5,
                            verbose = TRUE) {

  if (pfields::is.pField(field)) {
    coord.field <- pfields::GetLatLonField(field)
    lats <- coord.field$lat2d
    lons <- coord.field$lon2d
  } else if (pfields::is.pTs(field)) {
    lats <- pfields::GetLat(field)
    lons <- pfields::GetLon(field)
  } else {
    stop("'field' whether pField nor pTs object.")
  }

  lons[lons > 180] <- lons[lons > 180] - 360

  i <- which(lons >= min.lon & lons <= max.lon)

  lons <- lons[i]
  lats <- lats[i]

  j <- which(lats >= min.lat & lats <= max.lat)

  lons <- lons[j]
  lats <- lats[j]

  # distances from midpoint to x/y border of defined region
  if (verbose) {

    mid.lat <- (min.lat + max.lat) / 2
    mid.lon <- (min.lon + max.lon) / 2
    x <- geostools::GetDistance(mid.lat, mid.lon, mid.lat, max(lons))
    y <- geostools::GetDistance(mid.lat, mid.lon, min(lats), mid.lon)

    midpoint.border.distances <- c(x = x, y = y)
    cat(sprintf("Midpoint-border distances:\nx = %4.0f, y = %4.0f.\n",
                midpoint.border.distances["x"],
                midpoint.border.distances["y"]))
  }

  lons[lons < 0] <- lons[lons < 0] + 360
  coordinates <- data.frame(field.indices = i[j], lat = lats, lon = lons)

  return(coordinates)

}

