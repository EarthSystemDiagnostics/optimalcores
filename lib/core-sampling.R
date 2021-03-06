##
## aim:
## functions for the sampling of core combinations from climate model data
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores
##

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS

##' Convert a list into a matrix
##'
##' @param lst a list of length 'n', where each list element is a numeric vector.
##' @return an \code{n x m} matrix , where 'n' is the length of \code{lst} and
##'   'm' is the length of the longest element in \code{lst}.
##' @author Thomas Münch
list2mat <- function(lst) {
    
  t(sapply(lst, "[", seq(max(sapply(lst, length)))))
}

##' Grid cells within a ring bin
##'
##' Determine all grid cells that lie within a ring bin of defined inner and
##' outer radius with respect to a target site.
##'
##' @param dist.start numeric; the distance in km to a target site of the inner
##'   ring's radius.
##' @param delta.d numeric; the constant width of the ring bins in km.
##' @param distance.field numeric vector of distances (in km) relative to a
##'   target site of the grid cells which are used to define the ring bins
##'   around the target.
##' @return numeric vector of indices of \code{distance.field} which lie in
##'   between \code{dist.start} and \code{dist.start + delta.d}.
##' @author Thomas Münch
getRingGrids <- function(dist.start = 0, delta.d = 250, distance.field) {

  dist.min <- dist.start
  dist.max <- dist.start + delta.d
  
  idx <- which(distance.field >= dist.min &
               distance.field <= dist.max)

  return(idx)

}

##' Number of grid cells within radial bins
##'
##' Determine the number of grid cells that lie within radial distance bins of
##' specified inner and outer radius given a field of grid cell distances
##' relative to a target site.
##'
##' @param distance.field numeric vector of distances of the grid cells in an
##'   investigated model grid relative to a target site.
##' @param start inner radius of the innermost (first) radial distance bin.
##' @param end inner radius of the outermost (last) radial distance bin.
##' @param binsize radial width of the requested distance bins.
##' @return integer vector the same length as \code{start} with the number of
##'   grid cells that lie within each radial distance bin, starting with the bin
##'   from \code{start} to \code{start + binsize} to the bin from \code{end} to
##'   \code{end + binsize}.
##' @author Thomas Münch
getNumberOfCells <- function(distance.field, start, end, binsize) {

  inner <- seq(start, end, binsize)
  outer <- inner + binsize

  cell.number <- integer()
  for (i in 1 : length(inner)) {
    cell.number <- c(cell.number,
      length(which(distance.field >= inner[i] & distance.field < outer[i])))
  }

  return(cell.number)

}

##' Produce symmetric data matrix
##'
##' Produce a symmetric matrix of data sampled from the same bins in two
##' dimensions.
##'
##' @param input a list of two elements:
##'   * "bins": numeric vector with the set of possible bins from which the data
##'     could sampled.
##'   * "samples": a list or data frame of three elements/columns of equal
##'     length 'n', where the first (second) element are the actual bins from
##'     which the data has been sampled in the first (second) dimension, and the
##'     third element is the data.
##' @return a symmetric 'n' x 'n' matrix filled with the data.
##' @author Thomas Münch
data2matrix <- function(input) {

  res <- matrix(nrow = length(input$bins), ncol = length(input$bins))
  for (i in 1 : nrow(input$samples)) {

    m <- match(input$samples[[1]][i], input$bins)
    n <- match(input$samples[[2]][i], input$bins)

    res[m, n] <- input$samples[[3]][i]
  }

  # make matrix symmetric
  res[lower.tri(res)] <- t(res)[lower.tri(res)]

  # put bins as scale attribute
  attr(res, "scale") <- input$bins

  return(res)

}

# ------------------------------------------------------------------------------
# SAMPLING FUNCTIONS

##' Sample cores
##'
##' Sample all possible combinations for selecting either a single core (grid
##' cell) from a set of grid cells of a given climate field, or for selecting two
##' cores from the same or from two different sets of grid cells. For each core
##' combination, compute the average time series of the grid cells and the
##' correlation of this average time series with a target site time series.
##'
##' @param grid.index1 numeric vector with a set of grid cells of \code{field}
##'   from which the cores are sampled.
##' @param grid.index2 numeric vector with an optional second set of grid cells
##'   of \code{field} in case a combination of two cores shall be
##'   sampled. Defaults to \code{NULL}, i.e. sampling one or two cores only from
##'   the one set in \code{grid.index1}.
##' @param N integer; the number of cores to sample, either one or two.
##' @param field a \code{"pField"} or \code{"pTs"} object with a climate
##'   field from which the grid cell time series are sampled.
##' @param target numeric vector with the target site reference time series to
##'   which the sampled cores are correlated.
##' @return an \code{n x m} matrix with the correlation (first column) and the
##'   grid cell indices ('m' - 1 columns) for each sampled core combination,
##'   where 'n' is the total number of sampled combinations and 'm' - 1 is the
##'   number of sampled cores.
##' @author Thomas Münch
sampleCores <- function(grid.index1, grid.index2 = NULL, N = 1,
                        field, target) {

  if (is.null(grid.index2)) {

    i.grids <- combn(grid.index1, N)

    grid.cor <- apply(i.grids, 2, function(i) {

      if (N == 1) {
        x <- field[, i[1]]
      } else {
        x <- rowMeans(field[, c(i[1], i[2])])
      }
      c(cor(x, target, use = "pairwise"), i)
    })
    
  } else {

    N <- 2

    i.grids <- expand.grid(grid.index1, grid.index2)
    grid.cor <- apply(i.grids, 1, function(i) {

      x <- rowMeans(field[, c(i[1], i[2])])
      c(cor(x, target, use = "pairwise"), i)
    })

  }

  grid.cor <- t(grid.cor)
  colnames(grid.cor) <- c("cor", paste0("index", 1 : N))
  
  return(grid.cor)

}

##' Sample single core
##'
##' Sample all time series of single cores (grid cells) from a climate field
##' which fall in consecutive rings around a target site, and compute the
##' correlation of these cores with the target site reference time series.
##' 
##' @param max.dist numeric; the maximum ring bin distance (in km) to study,
##'   i.e. the last ring bin will span from \code{max.dist} to
##'   \code{max.dist + delta.d}.
##' @param delta.d numeric; the constant width of the ring bins in km.
##' @param field a \code{"pField"} or \code{"pTs"} object with a climate
##'   field from which the single core time series are sampled. The grid
##'   structure of the field must match the structure of \code{distance.field}.
##' @param target numeric vector with the target site reference time series to
##'   which the sampled cores are correlated.
##' @param distance.field numeric vector of the distances (in km) of the grid
##'   cells in \code{field} relative to the grid cell of the target site. The
##'   spatial structure of these distances must follow the structure of
##'   \code{field}.
##' @return a list of four elements:
##'   * "ring.distances": numeric vector with the mid-point distances in km
##'     of the sampled ring bins from the target.
##'   * "ring.combinations": an \code{n x 1} matrix with the indices of the ring
##'     bins from which the single cores have been sampled, where 'n' is the
##'     number of sampled ring bins.
##'   * "grid.indices": an \code{n x m x 1} matrix with the sampled grid cell
##'     indices, where 'n' is the number of sampled ring bins and 'm' is the
##'     highest number of grid cells found across all ring bins. For those bins
##'     with a lower number of grid cells, the remaining matrix cells are filled
##'     with \code{NA}.
##'   * "correlations": an \code{n x m} matrix with the sample correlations,
##'     where 'n' is the number of sampled ring bins and 'm' is the highest
##'     number of grid cells found across all ring bins. For those bins with a
##'     lower number of grid cells, the remaining matrix cells are filled with
##'     \code{NA}. 
##' @author Thomas Münch
sampleOneFromRings <- function(max.dist = 2000, delta.d = 250,
                               field, target, distance.field) {

  class(field) <- attr(field, "oclass")

  if (ncol(field) != length(distance.field)) {
    stop("Number of grid points in 'field' ",
         "does not match length of 'distance.field'. ",
         "Check your input.")
  }

  # define ring bins
  ring.dist <- seq(0, max.dist, delta.d)

  # get all grid cells within the ring segments
  sites <- sapply(ring.dist, getRingGrids, distance.field = distance.field,
                  delta.d = delta.d)
  
  # get correlations of grid cells with target
  tmp <- lapply(sites, sampleCores, field = field, target = target)

  # output

  sitesAsMatrix <- list2mat(sites)
  sitesAsArray <- array(sitesAsMatrix, dim = c(dim(sitesAsMatrix), 1))
#  sites.array[, , 1] <- sites.mat

  res <- list()
  
  res$ring.distances <- ring.dist + delta.d / 2
  res$ring.combinations <-  matrix(seq(ring.dist), ncol = 1)
  res$grid.indices <- sitesAsArray
  res$correlations <- list2mat(lapply(tmp, function(x) {x[, "cor"]}))

  return(res)

}

##' Sample two cores
##'
##' Sample all combinations of two cores (grid cells) from a climate field
##' which fall in either the same or in two different consecutive rings around a
##' target site, and compute the correlation of the average of these cores with
##' the target site reference time series. All possibilities of pairwise
##' combinations of two ring bins are used, including sampling two cores from
##' the same bin.
##'
##' @param max.dist numeric; the maximum bin ring distance (in km) to study,
##'   i.e. the last ring bin will span from \code{max.dist} to
##'   \code{max.dist + delta.d}.
##' @param delta.d numeric; the constant width of the ring bins in km.
##' @param field a \code{"pField"} or \code{"pTs"} object with a climate
##'   field from which the cores are sampled. The grid structure of the field
##'   must match the structure of \code{distance.field}.
##' @param target numeric vector with the target site reference time series to
##'   which the average of the sampled cores are correlated.
##' @param distance.field numeric vector of the distances (in km) of the grid
##'   cells in \code{field} relative to the grid cell of the target site. The
##'   spatial structure of these distances must follow the structure of
##'   \code{field}.
##' @param .fix.central logical; if \code{TRUE}, ring sampling only includes
##'   combining the central ring with itself and all other rings; else (the
##'   default) all unique pairwise ring combinations are sampled.
##' @return a list of four elements:
##'   * "ring.distances": numeric vector with the mid-point distances in km
##'     of the sampled ring bins from the target.
##'   * "ring.combinations": an \code{n x p} matrix with the indices of the ring
##'     bins from which the two cores have been sampled, where 'n' is the
##'     number of all possibilities for combining two ring bins and 'p' is the
##'     number of sampled cores (i.e. 2).
##'   * "grid.indices": an \code{n x m x p} matrix with the sampled grid cell
##'     indices, where 'n' is the number of all possibilities for combining two
##'     ring bins, 'm' is the highest number of grid cell combinations found
##'     across all ring bin combinations, and 'p' is the number of sampled cores
##'     (i.e. 2). For those ring bin combinations where the number of grid cell
##'     combinations is smaller than 'm', the remaining matrix cells are filled
##'     with \code{NA}.
##'   * "correlations": an \code{n x m} matrix with the sample correlations,
##'     where 'n' is the number of all possibilities for combining two ring
##'     bins and 'm' is the highest number of grid cell combinations found
##'     across all ring bin combinations. For those ring bin combinations where
##'     the number of grid cell combinations is smaller than 'm', the remaining
##'     matrix cells are filled with \code{NA}. 
##' @author Thomas Münch
sampleTwoFromRings <- function(max.dist = 2000, delta.d = 250,
                               field, target, distance.field,
                               .fix.central = FALSE) {

  class(field) <- attr(field, "oclass")

  if (ncol(field) != length(distance.field)) {
    stop("Number of grid points in 'field' ",
         "does not match length of 'distance.field'. ",
         "Check your input.")
  }

  # define ring bins
  ring.dist <- seq(0, max.dist, delta.d)

  # get all grid cells within the ring segments
  sites <- sapply(ring.dist, getRingGrids, distance.field = distance.field,
                  delta.d = delta.d)

  # get all unique pairwise combinations of ring segments
  if (.fix.central) {

    # all combinations for fixing the central ring
    ring.comb <- cbind(1, 1 : length(ring.dist))

  } else {

    # all possible combinations
    ring.comb <- arrangements::combinations(1 : length(ring.dist),
                                            k = 2, replace = TRUE)
  }
  
  # get correlation of grid cell pairs with target
  tmp <- apply(ring.comb, 1, function(i) {
    
    if (diff(i) == 0) {
      
      # sample two cores from the same ring            
      sampleCores(grid.index1 = sites[[i[1]]], N = 2,
                  field = field, target = target)

    } else {
      
      # sample two cores from two different rings
      sampleCores(grid.index1 = sites[[i[1]]],
                  grid.index2 = sites[[i[2]]],
                  N = 2, field = field, target = target)

    }
  })

  # output

  res <- list()

  res$ring.distances <- ring.dist + delta.d / 2
  res$ring.combinations <- ring.comb

  res$grid.indices <- abind::abind(
    list2mat(lapply(tmp, function(x) {x[, "index1"]})),
    list2mat(lapply(tmp, function(x) {x[, "index2"]})),
    along = 3)
  attr(res$grid.indices, "dimnames") <- NULL
  
  res$correlations <- list2mat(lapply(tmp, function(x) {x[, "cor"]}))

  return(res)

}

##' Sample N cores
##'
##' Monte Carlo sample a given number of possibilities of combining N cores
##' (grid cells) of a climate field from N ring bins, compute the average of
##' these N cores, and compute the correlation with a target site reference time
##' series. This is performed for all possible combinations of N ring bins,
##' including sampling from the same ring bin.
##' 
##' @param N integer; the number of cores to sample.
##' @param nmc integer; the number of times N cores are sampled from
##'   \code{field} for each ring bin combination. Each time, one core is
##'   sampled randomly from the grid cells of one of the N ring bins.
##' @param max.dist numeric; the maximum ring bin distance (in km) to study,
##'   i.e. the last ring bin will span from \code{max.dist} to
##'   \code{max.dist + delta.d}.
##' @param delta.d numeric; the constant width of the ring bins in km.
##' @param field a \code{"pField"} or \code{"pTs"} object with a climate
##'   field from which the cores are sampled. The grid structure of the field
##'   must match the structure of \code{distance.field}.
##' @param target numeric vector with the target site reference time series to
##'   which the average of the sampled cores are correlated.
##' @param distance.field numeric vector of the distances (in km) of the grid
##'   cells in \code{field} relative to the grid cell of the target site. The
##'   spatial structure of these distances must follow the structure of
##'   \code{field}.
##' @param ngroups specify an integer number of groups into which the ring bin
##'   combinations are subdivided, which forces the ring sampling to be executed
##'   successively in a loop over these groups. This setup can be necessary in
##'   order to limit RAM demand needed for running a large number of cores
##'   (\code{N} > 5) in combination with a large number of Monte Carlo
##'   iterations. Defaults to \code{NULL} which uses no grouping.
##' @param default.ring.combination an optional matrix with \code{N} columns
##'   where each row specifies a set of \code{N} indices for a certain ring bin
##'   combination. If this parameter is provided, only grid cells from these
##'   given ring bin combinations are sampled. Defaults to \code{NULL} which
##'   samples all available ring bin combinations.
##' @param .parallel logical; whether to parallelize the correlation
##'   computations of the individual ring bin combinations. Defaults to
##'   \code{TRUE}.
##' @param mc.cores integer; the number of cores to use for the parallel
##'   computation, i.e. at most how many child processes will be run
##'   simultaneously. The default \code{NULL} means to use the value from
##'   \code{parallel::detectCores()}.
##' @return a list with four elements:
##'   * "ring.distances": numeric vector with the mid-point distances in km
##'     of the sampled ring bins from the target.
##'   * "ring.combinations": an \code{n x p} matrix with the indices of the ring
##'     bins from which the cores have been sampled, where 'n' is the number of
##'     all possibilities for combining \code{N} ring bins and 'p' is the number
##'     of sampled cores (i.e. \code{N}).
##'   * "grid.indices": an \code{n x m x p} matrix with the sampled grid cell
##'     indices, where 'n' is the number of all possibilities for combining
##'     \code{N} ring bins, 'm' is the number of Monte Carlo samples
##'     (\code{nmc}), and 'p' is the number of sampled cores (i.e. \code{N}).
##'   * "correlations": an \code{n x m} matrix with the sample correlations,
##'     where 'n' is the number of all possibilities for combining \code{N} ring
##'     bins and 'm' is the number of Monte Carlo samples (\code{nmc}).
##' @author Thomas Münch
sampleNFromRings <- function(N = 2, nmc = 100, max.dist = 2000, delta.d = 250,
                             field, target, distance.field,
                             ngroups = NULL, default.ring.combination = NULL,
                             .parallel = TRUE, mc.cores = NULL) {

  class(field) <- attr(field, "oclass")

  if (ncol(field) != length(distance.field)) {
    stop("Number of grid points in 'field' ",
         "does not match length of 'distance.field'. ",
         "Check your input.")
  }

  # define helper functions

  sampleRingMonteCarlo <- function(combi, sites, nmc) {

    sapply(combi, function(x, sites) {
      sample(sites[[x]], nmc, replace = TRUE)}, sites)
  }

  sampleCorrelation <- function(x, grid.indices, field, target) {

    core.means <- apply(grid.indices[x, , ], 1, function(x, field) {
      rowMeans(field[, x])}, field = field)
    correlations <- c(cor(core.means, target, use = "pairwise"))

    return(correlations)

  }

  if (!length(ngroups)) ngroups <- 1

  # define ring bins
  ring.dist <- seq(0, max.dist, delta.d)

  # grid cells within each ring bin
  sites <- sapply(ring.dist, getRingGrids, distance.field = distance.field)

  if (!is.null(default.ring.combination)) {
    # use only given ring bin combinations
    if (ncol(default.ring.combination) != N) {
      stop("Number of default ring bins does not match core number to sample.")
    }
    ring.comb.full <- default.ring.combination
  } else {
    # use all possible combinations of ring bins for N cores
    ring.comb.full <- arrangements::combinations(x = 1 : length(ring.dist),
                                                 k = N, replace = TRUE)
  }

  # index vector of the number of ring bin combinations
  x <- seq(nrow(ring.comb.full))

  # split in nearly equal groups of indices
  groups <- split(x, sort(x %% ngroups))

  # loop over ring bin combination groups
  cor.lst <- list()
  for (i in 1 : length(groups)) {

    if (ngroups != 1) print(i)

    # combinations for this index group
    ring.comb <- ring.comb.full[groups[[i]], ]

    # fix matrix if needed
    if (is.null(dim(ring.comb))) ring.comb <- rbind(ring.comb)

    # obtain nmc Monte Carlo grid cell sets for each ring bin combination
    grid.indices <- lapply(split(ring.comb, row(ring.comb)),
                           sampleRingMonteCarlo,
                           sites = sites, nmc = nmc)
    # convert to array
    grid.indices <- aperm(aperm(
      abind::abind(grid.indices, along = 3),
      perm = c(3, 2, 1)), perm = c(1, 3, 2))
    attr(grid.indices, "dimnames") <- NULL

    # obtain correlation with target for each grid cell set
    if (.parallel) {

      require(parallel)
      if (!length(mc.cores)) mc.cores <- parallel::detectCores()

      correlations <- parallel::mclapply(seq(nrow(ring.comb)), sampleCorrelation,
                                         grid.indices = grid.indices,
                                         field = field, target = target,
                                         mc.cores = mc.cores)

    } else {

      correlations <- lapply(seq(nrow(ring.comb)), sampleCorrelation,
                             grid.indices = grid.indices,
                             field = field, target = target)
    }

    if (ngroups != 1) grid.indices <- NULL

    cor.lst[[i]] <- t(simplify2array(correlations))

  }

  correlations <- do.call(rbind, cor.lst)

  return(list(ring.distances = ring.dist + delta.d / 2,
              ring.combinations = ring.comb.full,
              grid.indices = grid.indices,
              correlations = correlations))

}

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
##'
##' @param region a data frame of coordinate indices and latitude-longitude
##'   values defining a subset region of the climate field in
##'   \code{target.field}.
##' @param target.field a \code{"pField"} object from which the grid cells
##'   defined in \code{region} are to be selected as target sites.
##' @param study.field a \code{"pField"} or \code{"pTs"} object with a climate
##'   field for which the correlations with the target sites are to be
##'   calculated. Its grid structure must match the structure of the distance
##'   field obtained from selecting a target site from the \code{target.field}.
##' @param N integer; the number of grid cells to average before computing the
##'   correlation to the target.
##' @param max.dist the inner radius of the outermost ring (in km).
##' @param delta.d the ring width (in km).
##' @param nmc number of Monte Carlo iterations when sampling more than two
##' sites, i.e. for \code{N} > 2.
##' @param ngroups [only for N > 2]: specify an integer number of groups into
##'   which the ring bin combinations are subdivided, which forces the ring
##'   sampling to be executed successively in a loop over these groups. This
##'   setup can be necessary in order to limit RAM demand needed for running a
##'   large number of cores (\code{N} > 5) in combination with a large number of
##'   Monte Carlo iterations. Defaults to \code{NULL} which uses no grouping.
##' @param default.ring.combination [only for N > 2]: an optional matrix with
##'   \code{N} columns where each row specifies a set of \code{N} indices for a
##'   certain ring bin combination. If this parameter is provided, only grid
##'   cells from these given ring bin combinations are sampled. Defaults to
##'   \code{NULL} which samples all available ring bin combinations.
##' @param .parallel [only for N > 2]: logical; whether to parallelize the
##'   correlation computations of the individual ring bin combinations. Defaults
##'   to \code{TRUE}.
##' @param mc.cores [only for N > 2]: integer; the number of cores to use for
##'   the parallel computation, i.e. at most how many child processes will be
##'   run simultaneously. The default \code{NULL} means to use the value from
##'   \code{parallel::detectCores()}.
##' @param verbose logical; if \code{TRUE}, print a progess message giving
##'   the number of the currently analysed target site of the \code{region} and
##'   its latitude and longitude; defaults to \code{FALSE}.
##' @return a list of the same length as the number of sites in \code{region};
##'   each list element contains the output of \code{processCores} being applied
##'   on the results of the respective regional site.
##' @author Thomas Münch
analyseTargetRegion <- function(region, target.field, study.field, N = 1,
                                max.dist = 2000, delta.d = 250, nmc = 100,
                                ngroups = NULL, default.ring.combination = NULL,
                                .parallel = TRUE, mc.cores = NULL,
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
      tmp <- sampleOneFromRings(max.dist = max.dist, delta.d = delta.d,
                                field = study.field,
                                target = target.site$dat,
                                distance.field = target.site$dist)
    } else if (N == 2) {
      tmp <- sampleTwoFromRings(max.dist = max.dist, delta.d = delta.d,
                                field = study.field,
                                target = target.site$dat,
                                distance.field = target.site$dist)
    } else {
      tmp <- sampleNFromRings(max.dist = max.dist, delta.d = delta.d,
                              N = N, nmc = nmc,
                              field = study.field,
                              target = target.site$dat,
                              distance.field = target.site$dist,
                              ngroups = ngroups, .parallel = .parallel,
                              mc.cores = mc.cores,
                              default.ring.combination =
                                default.ring.combination)
    }

    res[[i]] <- processCores(tmp)

  }

  return(res)

}

# ------------------------------------------------------------------------------
# PROCESSING FUNCTIONS

##' Process sampling results
##'
##' Process the results from sampling 1 to N cores from ring bins defined around
##' a target site.
##'
##' @param input the output of the functions 'sampleOneFromRings()',
##'   'sampleTwoFromRings()' or 'sampleNFromRings()'.
##' @param probs numeric; single value of the probability applied to calculate
##'   the lower (\code{probs}) and upper (\code{1 - probs}) quantile of each
##'   mean sample correlation.
##' @param n.optim an integer value for the number of maximum mean sample
##'   correlations to return. Set to \code{NULL} (the default) to return an
##'   upper quantile of correlations instead.
##' @param upper.quantile numeric; value for an upper quantile which
##'   defines an optimal set of all mean sample correlations; ignored if
##'   \code{n.optim} is not \code{NULL}.
##' @return a list of six elements:
##'   * "input": a copy of the \code{input}.
##'   * "N": number of cores averaged for sample correlations in \code{input}.
##'   * "ring.distances.sampled": a data frame with the mid-point distances of
##'     the ring bins from which each of the N cores has been sampled.
##'   * "correlation": a data frame of the sample correlations (mean value,
##'     lower quantile, upper quantile, and standard deviation) across all
##'     individual sampled core combinations for each ring bin combination.
##'   * "distances": a data frame with the mean and the standard deviation of
##'     the distance from the target of the ring bins sampled for each ring
##'     bin combination.
##'   * "optimal.rings": a list of four elements:
##'     (1) a subset of the "correlation" component with the
##'         \code{n.optim} highest or the upper quantile, as defined by
##'         \code{upper.quantile}, of the mean correlations.
##'     (2) the corresponding subset of the "distances" component.
##'     (3) the corresponding subset of the "ring.combinations" component.
##'     (4) "counts": a matrix with the numbers each ring bin has been sampled
##'         for each optimal mean correlation.
##' @author Thomas Münch
processCores <- function(input, probs = 1/3,
                         n.optim = NULL, upper.quantile = 0.95) {

  res <- list()

  #----------
  # return input
  res$input <- input

  #----------
  # number of sampled cores
  res$N <- ncol(input$ring.combinations)

  #----------
  # sampled ring distances
  res$ring.distances.sampled <- as.data.frame(
    matrix(input$ring.distances[input$ring.combinations], ncol = res$N))
  colnames(res$ring.distances.sampled) <- paste0("core", 1 : res$N)

  #----------
  # sample correlation average and variability
  res$correlation <- data.frame(

    # mean correlation and variability
    mean = apply(input$correlations, 1, mean, na.rm = TRUE),
    q1   = apply(input$correlations, 1, quantile,
                 probs = probs, na.rm = TRUE),
    q2   = apply(input$correlations, 1, quantile,
                 probs = 1 - probs, na.rm = TRUE),
    sd   = apply(input$correlations, 1, sd, na.rm = TRUE)
  )

  #----------
  # ring distance mean and variability
  res$distances <- data.frame(
    
    mean = apply(res$ring.distances.sampled, 1, mean),
    sd   = apply(res$ring.distances.sampled, 1, sd)
  )

  #----------
  # best n or upper quantile combinations

  res <- getOptimalCorrelations(res, n.optim = n.optim,
                                upper.quantile = upper.quantile)

  return(res)

}

##' Obtain optimal correlations
##'
##' Obtain the optimal correlations from a set of sampled ring bin
##' combinations. The optimal correlations are determined either from selecting
##' the first N highest correlation values, or from selecting an upper quantile
##' range of correlation values.
##'
##' @param data a list with the input data following the structure defined by
##'   \code{processCores}.
##' @param n.optim an integer value for the number of maximum correlations to
##'   return. Set to \code{NULL} (the default) to return an upper quantile of
##'   correlations instead.
##' @param upper.quantile numeric; the lower threshold of an upper quantile
##'   defining the optimal correlations; ignored if \code{n.optim} is not
##'   \code{NULL}.
##' @return a list of four elements:
##'   (1) a subset of the "correlation" component in \code{data} with the
##'       \code{n.optim} highest or the upper quantile, as defined by
##'       \code{upper.quantile}, of the mean correlations.
##'   (2) the corresponding subset of the "distances" component in \code{data}.
##'   (3) the corresponding subset of the "ring.combinations" component in
##'       \code{data}.
##'   (4) "counts": a matrix with the numbers each ring bin has been sampled
##'       for each optimal mean correlation.
##' @author Thomas Münch
getOptimalCorrelations <- function(data, n.optim = NULL,
                                   upper.quantile = 0.95) {

  if (!length(n.optim)) {
    corQuantile <- quantile(data$correlation$mean, probs = upper.quantile)
    i <- which(data$correlation$mean >= corQuantile)
  } else if (length(n.optim) == 1) {
    rank <- sort.int(data$correlation$mean, decreasing = TRUE,
                     index.return = TRUE)
    i <- rank$ix[1 : n.optim]
  } else {
    stop("Set single number of optimal cores to output.", call. = FALSE)
  }

  data$optimal.rings <- list(

  correlation  = data$correlation[i, ],
  distances    = data$distances[i, ],
  combinations = data$input$ring.combinations[i, , drop = FALSE]
  )

  data$optimal.rings$counts = t(
    apply(data$optimal.rings$combinations, 1, function(x) {
      hist(x, 1 : (length(data$input$ring.distances) + 1),
           plot = FALSE, right = FALSE)$counts})
  )

  return(data)

}

##' Process regional mean correlation structure
##'
##' Compute the average correlation across results obtained from sampling ring
##' bins for several target sites in a region.
##'
##' @param data a list with the results from running \code{analyseTargetRegion}.
##' @return a list of five elements:
##'   * "input": a 2-element list with the midpoint distances of the ring bins
##'     which were sampled and with the sampled ring bin combinations for each
##'     "N".
##'   * "N": the number of grid cells ("ice cores") averaged; identical to the
##'     number of combined rings for each possibility.
##'   * "ring.distances.sampled": data frame with the midpoint distances of the
##'     sampled and combined ring bins for each N (ice core). The number of
##'     distances (rows) for each N corresponds to the maximum number of
##'     possible ring bin combinations.
##'   * "distances": data frame with the mean and SD across the columns of
##'     "ring.distances.sampled".
##'   * "correlation": data frame with the mean and SD of the expected
##'     correlation for each ring bin combination from averaging across all
##'     target sites in the target region.
##' @author Thomas Münch
processRegionalMean <- function(data) {

  res <- list(

    input = list(
      ring.distances    = data[[1]]$input$ring.distances,
      ring.combinations = data[[1]]$input$ring.combinations
    ),

    N = data[[1]]$N,

    ring.distances.sampled = data[[1]]$ring.distances.sampled,

    distances = data[[1]]$distances,

    correlation = data.frame(
      mean = data %>%
        sapply(function(x) {x$correlation$mean}) %>%
        apply(1, mean),
      sd = data %>%
        sapply(function(x) {x$correlation$mean}) %>%
        apply(1, sd)
    )
  )

  return(res)

}

##' Prepare matrix conversion of double core results
##'
##' Prepare the matrix conversion of the results of sampling two cores from
##  consecutive ring bins, i.e. the output of this function is tailored to the
##  input expected by \code{data2matrix}.
##'
##' @param data the output of \code{processCores} or \code{processRegionalMean}
##'   run with the data from \code{sampleTwoFromRings} or from
##'   \code{analyseTargetRegion} with the parameter \code{N} set to 2.
##' @return a list of two elements:
##'   * "bins": numeric vector with the midpoint distances of the possible ring
##'     bins from which grid cells could be sampled.
##'   * "samples": a data frame with, for each ring bin combination, the
##'     midpoint distances of the sampled ring bins and the corresponding
##'     correlation with the target time series.
##' @author Thomas Münch
prepareMatrixConversion <- function(data) {

  list(
    bins    = data$input$ring.distances,
    samples = data.frame(
      data$ring.distances.sampled, cor = data$correlation$mean)
  )
}

##' Arrange ring bin occurrences
##'
##' Produce a matrix with the occurrences of sampled ring bins; more
##' specifically, their mid-point distances from the target site.
##'
##' @param ring.counts the \code{counts} component of the output of
##'   \code{processCores}.
##' @param ring.distances numeric vector of all ring bin mid-point distances.
##' @param dx horizontal offset in the same units as \code{ring.distances} to
##'   use when more than one core has been sampled from the same ring bin.
##' @return an 'n' x 'm' matrix with the sampled mid-point distances, where 'n'
##'   is the number of samples and 'm' a number equal or larger than the length
##'   of \code{ring.distances}, depending on the number of cores which have been
##'   sampled from the same bin.
##' @author Thomas Münch
arrangeRingOccurrences <- function(ring.counts, ring.distances, dx = 50) {

  arrangeRingGroup <- function(n, x, dx) {

    if (n == 0) return(NA)

    x <- x + seq(0, (n - 1) * dx, dx)
    x <- x - ((n - 1) / 2) * dx

    x
  }

  ring.occurrences <- apply(ring.counts, 1, function(x) {

    unlist(sapply(1 : length(x), function(i) {
      arrangeRingGroup(n = x[i], x = ring.distances[i], dx = dx)
    }))
  })

  if (is.list(ring.occurrences)) {
    list2mat(ring.occurrences)
  } else if (is.matrix(ring.occurrences)) {
    t(ring.occurrences)
  } else {
    stop("Something's wrong in 'arrangeRingOccurrences()': ",
         "variable 'ring.occurrences' is neither list nor matrix.",
         call. = FALSE)
  }

}
