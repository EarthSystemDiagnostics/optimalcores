#!/usr/bin/env Rscript

##
## aim:
## script to analyse the correlation structure with temperature for averaging
## two isotope cores assessed with expectation value across ring combinations.
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

# make script to also run on the terminal
if (interactive()) {

  # interactive usage
  source("init.R")
  cmd.arg <- "dml"

} else {

  # script usage
  pwd <- Sys.getenv()[["PWD"]]
  setwd(pwd)
  source("../setup.R")
  source("init.R")

  cmd.arg <- commandArgs(trailingOnly = TRUE)
  if (!length(cmd.arg)) {
    stop("Please specify a command line option for running the simulation.")
  }
}

# save the correlation data?
SAVE <- TRUE

# ------------------------------------------------------------------------------
# Select data and define region of target sites

model  <- selectData()

if (cmd.arg == "dml") {

  region <- setTargetRegion(field = model$lnd.t2m, verbose = FALSE)

} else if (cmd.arg == "vostok") {

  region <- setTargetRegion(field = model$lnd.t2m,
                            min.lat = -83.5, max.lat = -73.5,
                            min.lon = 89.5, max.lon = 124.5, verbose = FALSE)

} else {

  stop(sprintf("Unknown command line option '%s' to select analysis region.",
               cmd.arg))
}

# ------------------------------------------------------------------------------
# Expected correlation depending on distance for averaging two cores

cat("\n")
cat(as.character(Sys.time()), "\n")
cat("Running t2m...\n")

t2m <- region %>%
  analyseTargetRegion(target.field = model$t2m,
                      study.field = model$lnd.t2m, N = 2) %>%
  processRegionalMean() %>%
  prepareMatrixConversion() %>%
  data2matrix()

cat("Running t2m.pw...\n")

t2m.pw <- region %>%
  analyseTargetRegion(target.field = model$t2m,
                      study.field = model$lnd.t2m.pw, N = 2) %>%
  processRegionalMean() %>%
  prepareMatrixConversion() %>%
  data2matrix()

cat("Running oxy...\n")

oxy <- region %>%
  analyseTargetRegion(target.field = model$t2m,
                      study.field = model$lnd.oxy, N = 2) %>%
  processRegionalMean() %>%
  prepareMatrixConversion() %>%
  data2matrix()

cat("Running oxy.pw...\n")

oxy.pw <- region %>%
  analyseTargetRegion(target.field = model$t2m,
                      study.field = model$lnd.oxy.pw, N = 2) %>%
  processRegionalMean() %>%
  prepareMatrixConversion() %>%
  data2matrix()

if (SAVE) {

  saved <- list(
    t2m    = t2m,
    t2m.pw = t2m.pw,
    oxy    = oxy,
    oxy.pw = oxy.pw
  )
  attr(saved, "version") <- Sys.Date()

  saveRDS(saved,
          file = sprintf("analysis/ring_correlation_%s_N=2.rds", cmd.arg))

}

cat(as.character(Sys.time()), "\n")
cat("done.\n")

# ------------------------------------------------------------------------------
# Plotting

color.palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdYlBu")))
distances <- attr(t2m, "scale")
label <- c(expression("(" * bold("a") * ") " * italic("T")["2m"]),
           expression("(" * bold("b") * ") " * italic("T")["2m"]^{"(pw)"}),
           expression("(" * bold("x") * ") " * delta^{18} * "O"),
           expression("(" * bold("c") * ") " * delta^{18} * "O"^{"(pw)"}))

filebase <- sprintf("echam5_mpiom_wiso_two_core_correlation_%s_", cmd.arg)

if (cmd.arg == "dml") {
  zlim.oxy = c(0.2, 0.5)
  zlim.oxy.pw = c(0.15, 0.35)
} else {
  zlim.oxy = c(0.1, 0.55)
  zlim.oxy.pw = c(0.15, 0.45)
}

Quartz(dpi = 300, file = file.path(
  SAVEPATH, "main", paste0(filebase, "01_t2m.png")))
plotCorrelationContours(t2m, distances, color.palette, zlim = c(0.2, 1),
                        label = label[1])
dev.off()

Quartz(dpi = 300, file = file.path(
  SAVEPATH, "main", paste0(filebase, "02_t2m.pw.png")))
plotCorrelationContours(t2m.pw, distances, color.palette, zlim = c(0.2, 0.6),
                        label = label[2])
dev.off()

Quartz(dpi = 300, file = file.path(
  SAVEPATH, "main", paste0(filebase, "03_oxy.png")))
plotCorrelationContours(oxy, distances, color.palette, zlim = zlim.oxy,
                        label = label[3])
dev.off()

Quartz(dpi = 300, file = file.path(
  SAVEPATH, "main", paste0(filebase, "04_oxy.pw.png")))
plotCorrelationContours(oxy.pw, distances, color.palette, zlim = zlim.oxy.pw,
                        label = label[4])
dev.off()

