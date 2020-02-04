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

# save the correlation data?
SAVE <- TRUE

# ------------------------------------------------------------------------------
# Select data and define region of target sites

model  <- selectData()
dml    <- setTargetRegion(field = model$lnd.t2m, verbose = FALSE)

# ------------------------------------------------------------------------------
# Preliminary code to run sampling for EDML and Vostok sites with N = 3, 5

target <- setTarget(model$t2m, site = "edml")
n.optim <- 5

system.time(
  edml.n3 <- sampleNFromRings(N = 3, nmc = 1000,
                              field = model$lnd.oxy.pw,
                              target = target$dat,
                              distance.field = target$dist) %>%
    processCores(n.optim = n.optim)
)

system.time(
  edml.n5 <- sampleNFromRings(N = 5, nmc = 100,
                              field = model$lnd.oxy.pw,
                              target = target$dat,
                              distance.field = target$dist) %>%
    processCores(n.optim = n.optim)
)

target <- setTarget(model$t2m, site = "vostok")

system.time(
  vostok.n3 <- sampleNFromRings(N = 3, nmc = 1000,
                                field = model$lnd.oxy.pw,
                                target = target$dat,
                                distance.field = target$dist) %>%
    processCores(n.optim = n.optim)
)

system.time(
  vostok.n5 <- sampleNFromRings(N = 5, nmc = 100,
                                field = model$lnd.oxy.pw,
                                target = target$dat,
                                distance.field = target$dist) %>%
    processCores(n.optim = n.optim)
)

# ------------------------------------------------------------------------------
# Example code to run analysis across a region

system.time(
  test <- analyseTargetRegion(region = dml[1 : 3, ],
                              target.field = model$t2m,
                              study.field = model$lnd.t2m,
                              N = 3, nmc = 100,
                              verbose = TRUE) %>%
    processRegionalMean() %>%
    getOptimalCorrelations(n.optim = n.optim)
)

# ------------------------------------------------------------------------------
# Plot ring occurrences

adj  <- 0.985
padj <- 0.7

Quartz(width = 14, height = 6,
       file = file.path(SAVEPATH, "main", "echam5_mpiom_wiso_binning.pdf"))
op <- par(LoadGraphicsPar(mfcol = c(2, 2),
                          mar = c(0, 0, 0, 0),
                          oma = c(5, 7.5, 2.5, 0.5)))

plotRingOccurrences(edml.n5, xlab = "", xaxt = "n", ylab = "")
abline(h = 0.5, v = 2250,lwd = 5)
mtext("(a)", side = 3, line = -1, las = 0, adj = adj, padj = padj,
      cex = par()$cex.lab, font = 2)

plotRingOccurrences(edml.n3, ylab = "")
abline(v = 2250,lwd = 5)
mtext("(c)", side = 3, line = -1, las = 0, adj = adj, padj = padj,
      cex = par()$cex.lab, font = 2)

plotRingOccurrences(vostok.n5, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
abline(h = 0.5,lwd = 5)
mtext("(b)", side = 3, line = -1, las = 0, adj = adj, padj = padj,
      cex = par()$cex.lab, font = 2)

plotRingOccurrences(vostok.n3, ylab = "", yaxt = "n")
mtext("(d)", side = 3, line = -1, las = 0, adj = adj, padj = padj,
      cex = par()$cex.lab, font = 2)

mtext("Rank", side = 2, line = 3.5, cex = par()$cex.lab, las = 0,
      outer = TRUE, at = 0.5)
mtext("N = 3", side = 2, line = 5.5, cex = par()$cex.lab, las = 0,
      outer = TRUE, at = 0.25)
mtext("N = 5", side = 2, line = 5.5, cex = par()$cex.lab, las = 0,
      outer = TRUE, at = 0.75)

mtext("DML", side = 3, line = 0.25, cex = par()$cex.lab,
      outer = TRUE, at = 0.03)
mtext("Vostok", side = 3, line = 0.25, cex = par()$cex.lab,
      outer = TRUE, at = 0.54)

dev.off()
par(op)
