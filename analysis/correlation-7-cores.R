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
  edml <- sampleNFromRings(N = N, nmc = nmc, field = field,
                           target = target$dat,
                           distance.field = target$dist,
                           mc.cores = mc.cores, ngroups = 10) %>%
    processCores(n.optim = n.optim)
)

# ------------------------------------------------------------------------------
# Run sampling for Vostok site with N = 7

cat("Running Vostok...\n")

target <- setTarget(model$t2m, site = "vostok")

system.time(
  vost <- sampleNFromRings(N = N, nmc = nmc, field = field,
                           target = target$dat,
                           distance.field = target$dist,
                           mc.cores = mc.cores, ngroups = 10) %>%
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
