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
  attr(res, "scale") <- bins

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
##' the target site reference time series. All possibilities of combining two
##' ring bins are used, including sampling two cores from the same bin.
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

  # get all unique combinations of ring segments
  ring.comb <- arrangements::combinations(1 : length(ring.dist),
                                          k = 2, replace = TRUE)
  
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
                             field, target, distance.field) {

  class(field) <- attr(field, "oclass")

  if (ncol(field) != length(distance.field)) {
    stop("Number of grid points in 'field' ",
         "does not match length of 'distance.field'. ",
         "Check your input.")
  }

  # define ring bins
  ring.dist <- seq(0, max.dist, delta.d)
  # grid cells within each ring bin
  sites <- sapply(ring.dist, getRingGrids, distance.field = distance.field)

  # all possible combinations of ring bins for N cores
  ring.comb <- arrangements::combinations(x = 1 : length(ring.dist),
                                          k = N, replace = TRUE)

  correlations <- array(dim = c(nrow(ring.comb), nmc))
  grid.indices <- array(dim = c(nrow(ring.comb), nmc, N))

  for (i in 1 : nrow(ring.comb)) {
    for (j in 1 : nmc) {

      # run nmc Monte Carlo samples for given ring combination

      i.grids <- sapply(ring.comb[i, ], function(x) {
        sample(sites[[x]], 1)})

      correlations[i, j] <-
        cor(rowMeans(field[, i.grids]), target, use = "pairwise")

      grid.indices[i, j, ] <- i.grids

    }
  }

  return(list(ring.distances = ring.dist + delta.d / 2,
              ring.combinations = ring.comb,
              grid.indices = grid.indices,
              correlations = correlations))

}

# ------------------------------------------------------------------------------
# PROCESSING FUNCTION

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
##' @param upper.quantile numeric; value for an upper quantile which
##'   defines an optimal set of all mean sample correlations.
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
##'     (1) a subset of "correlation" with the upper quantile as defined by
##'         \code{upper.quantile} of the mean correlations.
##'     (2) the corresponding subset of "distances".
##'     (3) the corresponding subset of the "ring.combinations" in the input.
##'     (4) "counts": a matrix with the numbers each ring bin has been sampled
##'         for each mean correlation in the the upper quantile.
##' @author Thomas Münch
processCores <- function(input, probs = 1/3, upper.quantile = 0.95) {

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
    
    mean = apply(res$ring.distances, 1, mean),
    sd   = apply(res$ring.distances, 1, sd)
  )

  #----------
  # upper quantile combinations
  corQuantile <- quantile(res$correlation$mean, probs = upper.quantile)
  i <- which(res$correlation$mean >= corQuantile)
  
  res$optimal.rings <- list()

  res$optimal.rings$correlation <- res$correlation[i, ]
  res$optimal.rings$distances <- res$distances[i, ]

  res$optimal.rings$combinations <- input$ring.combinations[i, , drop = FALSE]

  res$optimal.rings$counts <- t(
    apply(res$optimal.rings$combinations, 1, function(x) {
      hist(x, 1 : (length(input$ring.distances) + 1),
           plot = FALSE, right = FALSE)$counts})
  )

  return(res)

}

##' Process regional mean correlation structure
##'
##' Compute the average correlation across results obtained from sampling ring
##' bins for several target sites in a region.
##'
##' @param input a list with the results from running \code{analyseTargetRegion}.
##' @return a list of two elements:
##'   * "bins": numeric vector with the midpoint distances of the possible ring
##'     bins from which grid cells could be sampled.
##'   * "samples": a data frame with, for each ring bin combination, the
##'     midpoint distances of the sampled ring bins and the average correlation
##'     across all analysed target sites.
##' @author Thomas Münch
processRegionalMean <- function(input) {

  res <- list(

       bins = input[[1]]$input$ring.distances,
    samples = input[[1]]$ring.distances.sampled
  )

  res$samples$cor = input %>%
    sapply(function(x) {x$correlation$mean}) %>%
    apply(1, mean)

  return(res)

}
