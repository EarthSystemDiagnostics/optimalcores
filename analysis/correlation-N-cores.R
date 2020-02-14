#!/usr/bin/env Rscript

##
## aim:
## script to analyse the sampling structure of rings around the target for
## averaging N isotope cores, assessed with expectation value across ring
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
# Select data and set paramters

model  <- selectData()
field  <- model$lnd.oxy.pw

n.optim <- 5
nmc     <- 10^5

mc.cores <- ceiling(0.5 * parallel::detectCores())

# ------------------------------------------------------------------------------
# Run sampling for EDML site with N = 3, 5

cat("\n")
cat(as.character(Sys.time()), "\n")
cat("Running EDML...\n")

target <- setTarget(model$t2m, site = "edml")

system.time(
  edml.n3 <- sampleNFromRings(N = 3, nmc = nmc, field = field,
                              target = target$dat,
                              distance.field = target$dist,
                              mc.cores = mc.cores) %>%
    processCores(n.optim = n.optim)
)

system.time(
  edml.n5 <- sampleNFromRings(N = 5, nmc = nmc, field = field,
                              target = target$dat,
                              distance.field = target$dist,
                              mc.cores = mc.cores) %>%
    processCores(n.optim = n.optim)
)

# ------------------------------------------------------------------------------
# Run sampling for Vostok site with N = 3, 5

cat("Running Vostok...\n")

target <- setTarget(model$t2m, site = "vostok")

system.time(
  vost.n3 <- sampleNFromRings(N = 3, nmc = nmc, field = field,
                              target = target$dat,
                              distance.field = target$dist,
                              mc.cores = mc.cores) %>%
    processCores(n.optim = n.optim)
)

system.time(
  vost.n5 <- sampleNFromRings(N = 5, nmc = nmc, field = field,
                              target = target$dat,
                              distance.field = target$dist,
                              mc.cores = mc.cores) %>%
    processCores(n.optim = n.optim)
)

# ------------------------------------------------------------------------------
# Save output

if (SAVE) {

  # strip off input data to save disk space
  edml.n3$input <- NULL
  edml.n5$input <- NULL
  vost.n3$input <- NULL
  vost.n5$input <- NULL

  saved <- list(
    edml = list(N3 = edml.n3, N5 = edml.n5),
    vost = list(N3 = vost.n3, N5 = vost.n5)
  )
  attr(saved, "version") <- Sys.Date()

  saveRDS(saved, file = "analysis/ring_occurrences_edml_vostok_N=3_N=5.rds")

}

cat(as.character(Sys.time()), "\n")
cat("done.\n")
