#!/usr/bin/env Rscript

##
## aim:
## script to analyse the sampling structure of rings around the target for
## averaging N isotope cores, assessed with expectation value across ring
## combinations using Monte Carlo sampling of grid cells.
##
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Münch, Werner and Laepple, Clim. Past, 2021
##
## Terminal usage:
## ./correlation-N-Cores.R <N> <noptim>
## where <N> is the number of isotope cores to average and <noptim> is the
## number of best optimal solutions to retain (usually 1, 5 for paper Fig. 7).
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
  cat("\n")

  cmd.arg <- commandArgs(trailingOnly = TRUE)
  if (!length(cmd.arg)) {
    stop("Please specify command line options for running the simulation.")
  }
  if (length(cmd.arg) != 2) {
    stop("Two command line options expected: ",
         "1. number of cores, 2. number of best solutions to retain.")
  }
  N <- as.integer(cmd.arg[1])
  n.optim <- as.integer(cmd.arg[2])
  cat("Following command line options found:\n")
  cat(sprintf("Number of cores: %i\n", N))
  cat(sprintf("Number of best solutions: %i\n", n.optim))
}

library(parallel)

# save the correlation data?
SAVE <- TRUE

# ------------------------------------------------------------------------------
# Select data and set parameters

model <- selectData()
field <- model$lnd.oxy.pw

if (interactive()) {
  N <- 3
  n.optim <- 5
}

nmc <- 10^5

mc.cores <- ceiling(0.75 * parallel::detectCores())

ngroups <- NULL
if (N > 5) ngroups <- 10
if (N > 7) ngroups <- 50

file <- sprintf("analysis/ring_occurrences_edml_vostok_N=%i.rds", N)

# ------------------------------------------------------------------------------
# Run sampling for EDML site

cat("\n")
cat(as.character(Sys.time()), "\n\n")
cat("Running EDML...\n")

target <- setTarget(model$t2m, site = "edml")

system.time(
  edml <- sampleNFromRings(N = N, nmc = nmc, field = field,
                           target = target$dat,
                           distance.field = target$dist,
                           mc.cores = mc.cores, ngroups = ngroups) %>%
    processCores(n.optim = n.optim)
)

# strip off input data to save disk space and save intermediate results
if (SAVE) {

  edml$input <- NULL

  saveRDS(edml, file = "analysis/tmp_edml.rds")

  rm(edml)

}
cat("\n")

# ------------------------------------------------------------------------------
# Run sampling for Vostok site

cat("Running Vostok...\n")

target <- setTarget(model$t2m, site = "vostok")

system.time(
  vost <- sampleNFromRings(N = N, nmc = nmc, field = field,
                           target = target$dat,
                           distance.field = target$dist,
                           mc.cores = mc.cores, ngroups = ngroups) %>%
    processCores(n.optim = n.optim)
)

# strip off input data to save disk space and save intermediate results
if (SAVE) {

  vost$input <- NULL

  saveRDS(vost, file = "analysis/tmp_vost.rds")

  rm(vost)

}
cat("\n")

# ------------------------------------------------------------------------------
# Save output

if (SAVE) {

  saved <- list(
    edml = readRDS(file = "analysis/tmp_edml.rds"),
    vost = readRDS(file = "analysis/tmp_vost.rds")
  )
  attr(saved, "version") <- Sys.Date()

  saveRDS(saved, file = file)

}

cat("done.\n\n")
cat(as.character(Sys.time()), "\n\n")
