##
## aim:
## functions for the random picking of cores from a climate model data field
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores
##

##' Pick cores
##'
##' Pick N cores (grid cells) randomly from a given region of a climate model
##' field and correlate their average time series with a target time series.
##'
##' @param N integer; the number of grid cells to pick.
##' @param target a numeric vector with the target time series.
##' @param field a spatial field (matrix) of class \code{"pField"} or
##'   \code{"pTs"} with the climate model data from which grid cells are to be
##'   picked.
##' @param picking.sites a vector of the grid indices which define the region of
##'   \code{field} from which grid cells are to be picked.
##' @param nmc integer; number of random picks.
##' @param replace logical; if \code{TRUE} radnom picking is performed with
##'   replacement. Defaults to \code{FALSE} (no replacement).
##' @param COSTFUN a cost function to determine a metric for the agreement
##'   between the averaged cores and the target; defaults to the Pearson
##'   correlation coefficient (\code{\link[stats]{cor}}).
##' @param ... further arguments passed on to \code{COSTFUN}.
##' @param maximize logical; should the determined metric for each picking be
##'   sorted in descending (\code{TRUE}, the default) or ascending order upon
##'   output?
##' @param return.all logical; should all picking results be returned
##'   (\code{TRUE}, the default) or only the optimal one (according to the
##'   setting of \code{maximize})?
##' @param verbose logical; switch on/off printing of info messages.
##' @param .parallel logical; whether to parallelize the correlation
##'   computations of the individual sampling possibilities. Defaults to
##'   \code{FALSE}.
##' @param mc.cores integer; the number of cores to use for a parallel
##'   computation, i.e. at most how many child processes will be run
##'   simultaneously. The default \code{NULL} means to use the value from
##'   \code{parallel::detectCores()}.
##' @return a list of three elements:
##'   * "N": integer; the number of picked cores.
##'   * "sample": a list of \code{nmc} data frames (if \code{return.all = TRUE})
##'     with the coordinates of all \code{N} picked cores, or a single data
##'     frame (if \code{return.all = FALSE}) with the coordinates of only the
##'     optimal set of \code{N} cores.
##'   * "metric": a numeric vector of length \code{nmc} (if \code{return.all =
##'     TRUE}) with the metrics from \code{COSTFUN} for each picking, or a
##'     length-one vector with the metric of the optimal picking (if
##'     \code{return.all = FALSE}).
##' @author Thomas Münch
doPicking <- function(N = 1, target, field, picking.sites,
                      nmc = 1000, replace = FALSE,
                      COSTFUN = cor, ..., maximize = TRUE,
                      return.all = TRUE, verbose = TRUE,
                      .parallel = FALSE, mc.cores = NULL) {

  if (pfields::is.pField(field)) {

    coord.field <- pfields::GetLatLonField(field)
    lats <- coord.field$lat2d
    lons <- coord.field$lon2d

  } else if (pfields::is.pTs(field)) {

    lats <- pfields::GetLat(field)
    lons <- pfields::GetLon(field)

  } else {

    stop("Unknown class of input field.")

  }

  class(field) <- attr(field, "oclass")

  sav <- list(
    N = N,
    sample = list(),
    metric = vector()
  )

  samplePicks <- function(picked.sites, target, field, lats, lons, ...) {

    y <- rowMeans(field[, picked.sites, drop = FALSE], na.rm = TRUE)

    list(
      data.frame(
        index = picked.sites,
        lat = lats[picked.sites],
        lon = lons[picked.sites]
      ),

      COSTFUN(target, y, ...)
    )

  }

  if (N == 1) {

    if (verbose) cat("Sampling all possibilities.\n")

    res <- lapply(picking.sites, samplePicks, target, field, lats, lons, ...)
    sav$sample <- lapply(res, function(x) {x[[1]]})
    sav$metric <- sapply(res, function(x) {x[[2]]})

  } else {

    nPossibilities <- choose(length(picking.sites), N)

    if (nPossibilities > nmc) {

      if (verbose) {
        cat("Too many possibilities; resorting to Monte Carling sampling.\n")
      }

      picking.sites.sampled <- lapply(seq(nmc), function(x) {
        sample(picking.sites, size = N, replace = replace)})

    } else {

      if (verbose) cat("Sampling all possibilities.\n")

      picking.sites.sampled <- arrangements::combinations(picking.sites, k = N)
      picking.sites.sampled <- split(picking.sites.sampled,
                                     row(picking.sites.sampled))
    }

    if (.parallel) {

      require(parallel)
      if (!length(mc.cores)) mc.cores <- parallel::detectCores()

      res <- parallel::mclapply(picking.sites.sampled, samplePicks,
                                target, field, lats, lons, ...,
                                mc.cores = 3)

    } else {

      res <- lapply(picking.sites.sampled, samplePicks,
                    target, field, lats, lons, ...)
    }

    sav$sample <- lapply(res, function(x) {x[[1]]})
    sav$metric <- sapply(res, function(x) {x[[2]]})

  }

  decreasing <- ifelse(maximize, TRUE, FALSE)
  rank <- sort.int(sav$metric, decreasing = decreasing, index.return = TRUE)

  sav$sample <- sav$sample[rank$ix]
  sav$metric <- sav$metric[rank$ix]

  if (!return.all) {
    sav$sample <- sav$sample[[1]]
    sav$metric <- sav$metric[1]
  }

  return(sav)

}

##' Pick a set of N cores from a circle around a target
##'
##' Pick a set of N cores (grid cells) randomly from a circular region of a
##' climate model field around a target site and correlate, for each N, their
##' average time series with the time series of the target variable at the
##' target site (i.e. centre of the circle).
##'
##' @param N an integer vector with a set of N cores to pick.
##' @param target a character string with the name of the target site, as
##'   used by \code{\link{setTarget}}, which defines the centre of the circle
##'   from which the cores are picked; set to \code{NULL} to directly specify
##'   target coordinates.
##' @param lat0 numeric; (optional) latitude of the target site if no
##'   \code{target} name is specified.
##' @param lon0 numeric; (optional) longitude of the target site if no
##'   \code{target} name is specified.
##' @param radius numeric; the radius (in km) of the circle around the target
##'   site from which the cores are to be picked.
##' @param target.field a \code{"pField"} object with the spatial field (matrix)
##'   of a climate model variable from which the target site time series is
##'   extracted.
##' @param study.field a spatial field (matrix) of class \code{"pField"} or
##'   \code{"pTs"} with the climate model data from which the cores (grid cells)
##'   are to be picked.
##' @param nmc integer; number of random picks to perform for each N.
##' @param replace logical; if \code{TRUE} radnom picking is performed with
##'   replacement. Defaults to \code{FALSE} (no replacement).
##' @param COSTFUN a cost function to determine a metric for the agreement
##'   between the averaged cores and the target; defaults to the Pearson
##'   correlation coefficient (\code{\link[stats]{cor}}).
##' @param ... further arguments passed on to \code{COSTFUN}.
##' @param maximize logical; should the determined metric for each picking be
##'   sorted in descending (\code{TRUE}, the default) or ascending order upon
##'   output?
##' @param return.all logical; should all picking results (number of \code{nmc}
##'   for each N) be returned (\code{TRUE}, the default) or only the optimal one
##'   (according to the setting of \code{maximize})?
##' @param verbose logical; switch on/off printing of progress message.
##' @param .parallel logical; whether to parallelize the correlation
##'   computations of the individual sampling possibilities. Defaults to
##'   \code{FALSE}.
##' @param mc.cores integer; the number of cores to use for a parallel
##'   computation, i.e. at most how many child processes will be run
##'   simultaneously. The default \code{NULL} means to use the value from
##'   \code{parallel::detectCores()}.
##' @return a list of four elements:
##'   * "correlation.map": a \code{"pField"} or \code{"pTs"} object with the
##'     correlation between every grid cell in \code{field} and the
##'     \code{target} time series.
##'   * "target": a data frame with the name, latitude and longitude of the
##'     target site grid cell.
##'   * "radius": the radius (in km) of the circle around the target
##'     site from which the cores were picked.
##'   * "picking": a list the same length as \code{N}; each list element is the
##'     output of \code{\link{doPicking}}.
##' @author Thomas Münch
pickNCores <- function(N = 1, target = "edml", lat0 = NULL, lon0 = NULL,
                       radius = 2000, target.field, study.field,
                       nmc = 1000, replace = FALSE,
                       COSTFUN = cor, ..., maximize = TRUE,
                       return.all = TRUE, verbose = TRUE,
                       .parallel = FALSE, mc.cores = NULL) {

  target.site <- setTarget(target.field, site = target,
                           lat0 = lat0, lon0 = lon0)

  correlation.map <- pfields::cor.pTs(target.site$dat, study.field)

  target.coord <- data.frame(
    name = if (is.null(target)) "site i" else target,
    lat  = target.site$lat0,
    lon  = target.site$lon0
  )

  picking.sites <- which(target.site$dist <= radius)

  res <- list(
    correlation.map = correlation.map,
    target          = target.coord,
    radius          = radius,
    picking         = list()
  )

  for (i in 1 : length(N)) {

    if (verbose) cat(sprintf("Number of cores: %i\n", N[i]))
    res$picking[[i]] <- doPicking(N = N[i], target = target.site$dat,
                                  field = study.field,
                                  picking.sites = picking.sites,
                                  nmc = nmc, replace = replace,
                                  COSTFUN = COSTFUN, ...,
                                  maximize = maximize, return.all = return.all,
                                  verbose = verbose, .parallel = .parallel,
                                  mc.cores = mc.cores)

  }

  return(res)

}
      
