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

} else {

  # script usage
  pwd <- Sys.getenv()[["PWD"]]
  setwd(pwd)
  source("../setup.R")
  source("init.R")
}

# save the correlation data?
SAVE <- TRUE

# ------------------------------------------------------------------------------
# Select data and define region of target sites

model  <- selectData()
dml    <- setTargetRegion(field = model$lnd.t2m, verbose = FALSE)

# ------------------------------------------------------------------------------
# Expected correlation depending on distance for averaging two cores

cat("\n")
cat(as.character(Sys.time()), "\n")
cat("Running t2m...\n")

t2m <- dml %>%
  analyseTargetRegion(target.field = model$t2m,
                      study.field = model$lnd.t2m, N = 2) %>%
  processRegionalMean() %>%
  prepareMatrixConversion() %>%
  data2matrix()

cat("Running t2m.pw...\n")

t2m.pw <- dml %>%
  analyseTargetRegion(target.field = model$t2m,
                      study.field = model$lnd.t2m.pw, N = 2) %>%
  processRegionalMean() %>%
  prepareMatrixConversion() %>%
  data2matrix()

cat("Running oxy...\n")

oxy <- dml %>%
  analyseTargetRegion(target.field = model$t2m,
                      study.field = model$lnd.oxy, N = 2) %>%
  processRegionalMean() %>%
  prepareMatrixConversion() %>%
  data2matrix()

cat("Running oxy.pw...\n")

oxy.pw <- dml %>%
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

  saveRDS(saved, file = "analysis/ring-correlation_N=2.rds")

}

cat(as.character(Sys.time()), "\n")
cat("done.\n")

# ------------------------------------------------------------------------------
# Plotting

color.palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdYlBu")))
distances <- attr(t2m, "scale")
label <- c(expression(bold("(a) ") * italic("T")["2m"]),
           expression(bold("(b) ") * italic("T")["2m"]^{"(pw)"}),
           expression(bold("(x) ") * delta^{18} * "O"),
           expression(bold("(c) ") * delta^{18} * "O"^{"(pw)"}))

Quartz(dpi = 300, file = file.path(
  SAVEPATH, "main", "echam5_mpiom_wiso_two_core_correlation_edml_01_t2m.png"))
plotCorrelationContours(t2m, distances, color.palette, zlim = c(0.2, 1),
                        label = label[1])
dev.off()

Quartz(dpi = 300, file = file.path(
  SAVEPATH, "main", "echam5_mpiom_wiso_two_core_correlation_edml_02_t2m.pw.png"))
plotCorrelationContours(t2m.pw, distances, color.palette, zlim = c(0.2, 0.55),
                        label = label[2])
dev.off()

Quartz(dpi = 300, file = file.path(
  SAVEPATH, "main", "echam5_mpiom_wiso_two_core_correlation_edml_03_oxy.png"))
plotCorrelationContours(oxy, distances, color.palette, zlim = c(0.2, 0.5),
                        label = label[3])
dev.off()

Quartz(dpi = 300, file = file.path(
  SAVEPATH, "main", "echam5_mpiom_wiso_two_core_correlation_edml_04_oxy.pw.png"))
plotCorrelationContours(oxy.pw, distances, color.palette, zlim = c(0.15, 0.35),
                        label = label[4])
dev.off()

