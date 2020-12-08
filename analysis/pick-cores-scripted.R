#!/usr/bin/env Rscript

##
## aim:
## script to randomly pick cores (grid cells) from a climate model field to find
## the optimal sites to reconstruct a target time series. Intended to be used in
## batch mode from the terminal only!
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

# script usage
pwd <- Sys.getenv()[["PWD"]]
setwd(pwd)
source("../setup.R")
source("init.R")

library(parallel)

# save the correlation data?
SAVE <- TRUE

# ------------------------------------------------------------------------------
# Select data and set paramters

model  <- selectData()

N <- c(1, 3, 5)
nmc <- 1e7

.parallel <- TRUE
mc.cores <- ceiling(0.75 * parallel::detectCores())

# ------------------------------------------------------------------------------
# Do the picking for EDML and Vostok

cat("\n")
cat(as.character(Sys.time()), "\n")
cat("Running EDML...\n")

edml.picks <- pickNCores(N = N, target = "edml",
                         target.field = model$t2m,
                         study.field = model$lnd.oxy.pw,
                         nmc = nmc, return.all = FALSE,
                         .parallel = .parallel, mc.cores = mc.cores)

cat("Running Vostok...\n")

vostok.picks <- pickNCores(N = N, target = "vostok",
                           target.field = model$t2m,
                           study.field = model$lnd.oxy.pw,
                           nmc = nmc, return.all = FALSE,
                           .parallel = .parallel, mc.cores = mc.cores)

if (SAVE) {

  saved <- list(
    edml.picks   = edml.picks,
    vostok.picks = vostok.picks
  )
  attr(saved, "version") <- Sys.Date()

  saveRDS(saved, file = "analysis/picking.rds")

}

cat(as.character(Sys.time()), "\n")
cat("done.\n")
