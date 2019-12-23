##' Analyse a region
##'
##' Analyse all sites of a climate field within a given latitude-longitude
##' region for the average correlation of the mean of N time series of a climate
##' variable with the time series of a target climate field at a target
##' site. The N time series are sampled from consecutive ring bins.
##'
##' If N = 1, the time series are sampled from the individual rings. For N > 1,
##' the time series are sampled from rings iterating through all possibilities
##' of combining N rings.
##' @param region a data frame of coordinate indices and latitude-longitude
##'   values defining a subset region of the climate field in
##'   \code{target.field}.
##' @param target.field a \code{"pField"} object from which the grid cells
##'   defined in \code{region} are to be selected as target sites.
##' @param study.field a \code{"pTs"} object with a climate field for which the
##'   correlations with the target sites are to be calculated.
##' @param N integer; the number of grid cells to average before computing the
##'   correlation to the target.
##' @param max.dist the maximum distance of the outer ring (in km).
##' @param delta.d the ring width (in km).
##' @param nmc number of Monte Carlo iterations when sampling more than two
##' sites, i.e. for \code{N} > 2.
##' @param verbose logical; if \code{TRUE}, print a progess message giving
##'   the number of the currently analysed target site of the \code{region} and
##' its latitude and longitude; defaults to \code{FALSE}.
##' @return The return value depends on the number of averaged sites: if {N} =
##'   1, it is the output of \code{\link{processSingles}}, else the output of
##'   \code{\link{processNCores}}.
##' @author Thomas Münch
analyseTargetRegion <- function(region, target.field, study.field, N = 1,
                                max.dist = 2000, delta.d = 250, nmc = 100,
                                verbose = FALSE) {
    
    res <- list()
    for (i in 1 : nrow(region)) {

        target.site <- setTarget(field = target.field, site = NULL,
                                 lat0 = region$lat[i],
                                 lon0 = region$lon[i])

        if (verbose) {
            cat(sprintf("Run %i/%i:\n", i, nrow(region)))
            cat(sprintf("lat = %2.2f\n", target.site$lat0))
            cat(sprintf("lon = %2.2f\n", target.site$lon0))
            cat("\n")
        }

        if (N == 1) {
            tmp <- SampleOneFromRings(max.dist = max.dist, delta.d = delta.d,
                                      field = study.field,
                                      target.lst = target.site)
        } else if (N == 2) {
            tmp <- SampleTwoFromRings(field = study.field,
                                      target.lst = target.site)
        } else {
            tmp <- SampleNFromRings(N = N, field = study.field,
                                    target.lst = target)
        }

        if (N == 1) {
            res[[i]] <- ProcessSingles(tmp)
        } else {
            res[[i]] <- ProcessNCores(tmp)
        }

    }

    return(res)
    
}
