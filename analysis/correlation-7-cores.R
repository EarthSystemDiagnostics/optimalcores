#!/usr/bin/env Rscript

##
## aim:
## script to analyse the sampling structure of rings around the target for
## averaging 10 isotope cores, assessed with expectation value across ring
## combinations using Monte Carlo sampling of grid cells.
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

# make script to also run on the terminal
if (interactive()) {

  # interactive usage
  source("init.R")

} else {

  # script usage
  pwd <- Sys.getenv()[["PWD"]]
  setwd(pwd)
  source("../setup.R")
  source("init.R")
}

library(parallel)

# save the correlation data?
SAVE <- TRUE

# ------------------------------------------------------------------------------
# Custom version of 'sampleNFromRings' to reduce RAM demand

sampleNFromRings.custom <- function(N = 2, nmc = 100,
                                    max.dist = 2000, delta.d = 250,
                                    field, target, distance.field,
                                    .parallel = TRUE, mc.cores = NULL,
                                    ngroups = 10) {

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

  # define ring bins
  ring.dist <- seq(0, max.dist, delta.d)

  # grid cells within each ring bin
  sites <- sapply(ring.dist, getRingGrids, distance.field = distance.field)

  # all possible combinations of ring bins for N cores
  ring.comb.full <- arrangements::combinations(x = 1 : length(ring.dist),
                                               k = N, replace = TRUE)

  # index vector of the number of ring bin combinations
  x <- seq(nrow(ring.comb.full))

  # split in nearly equal groups of indices
  groups <- split(x, sort(x %% ngroups))

  # loop over ring bin combination groups
  cor.lst <- list()
  for (i in 1 : length(groups)) {

    print(i)

    # combinations for this index group
    ring.comb <- ring.comb.full[groups[[i]], ]

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

    grid.indices <- NULL

    cor.lst[[i]] <- t(simplify2array(correlations))

  }

  correlations <- do.call(rbind, cor.lst)

  return(list(ring.distances = ring.dist + delta.d / 2,
              ring.combinations = ring.comb.full,
              grid.indices = grid.indices,
              correlations = correlations))

}

# ------------------------------------------------------------------------------
# Select data and set paramters

model  <- selectData()
field  <- model$lnd.oxy.pw

n.optim <- 1
nmc     <- 10^5

mc.cores <- ceiling(0.75 * parallel::detectCores())

N <- 7

# ------------------------------------------------------------------------------
# Run sampling for EDML site with N = 7

cat("\n")
cat(as.character(Sys.time()), "\n")
cat("Running EDML...\n")

target <- setTarget(model$t2m, site = "edml")

system.time(
  edml <- sampleNFromRings.custom(N = N, nmc = nmc, field = field,
                                  target = target$dat,
                                  distance.field = target$dist,
                                  mc.cores = mc.cores) %>%
    processCores(n.optim = n.optim)
)

# ------------------------------------------------------------------------------
# Run sampling for Vostok site with N = 7

cat("Running Vostok...\n")

target <- setTarget(model$t2m, site = "vostok")

system.time(
  vost <- sampleNFromRings.custom(N = N, nmc = nmc, field = field,
                                  target = target$dat,
                                  distance.field = target$dist,
                                  mc.cores = mc.cores) %>%
    processCores(n.optim = n.optim)
)

# ------------------------------------------------------------------------------
# Save output

if (SAVE) {

  # strip off input data to save disk space
  edml$input <- NULL
  vost$input <- NULL

  saved <- list(
    edml = edml,
    vost = vost
  )
  attr(saved, "version") <- Sys.Date()

  saveRDS(saved, file = "analysis/ring_occurrences_edml_vostok_N=7.rds")

}

cat(as.character(Sys.time()), "\n")
cat("done.\n")
