##
## aim:
## script to plot the occurrences of sampled rings, i.e. the number of cores put
## into each ring; for different N and different target sites.
##
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Münch, Werner and Laepple, Clim. Past, 2021
##

source("init.R")

# ------------------------------------------------------------------------------
# Settings and load data

adj  <- 0.985
padj <- 0.5

dat.N3 <- readRDS("analysis/ring_occurrences_edml_vostok_N=3.rds")
dat.N5 <- readRDS("analysis/ring_occurrences_edml_vostok_N=5.rds")

# ------------------------------------------------------------------------------
# Plotting

op <- grfxtools::Quartz(width = 14, height = 6,
                        file = file.path(SAVEPATH, "main", "fig_07.pdf"),
                        mfcol = c(2, 2), mar = c(0, 0, 0, 0),
                        oma = c(5, 7.5, 2.5, 0.5))

plotRingOccurrences(dat.N3$edml, xlab = "", xaxt = "n", ylab = "")
abline(h = 0.5, v = 2250, lwd = 5)
mtext("a", side = 3, line = -1, las = 0, adj = adj, padj = padj,
      cex = par()$cex.lab, font = 2)

plotRingOccurrences(dat.N5$edml,
                    xlab = "Distance from target site (km)", ylab = "")
abline(v = 2250, lwd = 5)
mtext("c", side = 3, line = -1, las = 0, adj = adj, padj = padj,
      cex = par()$cex.lab, font = 2)

plotRingOccurrences(dat.N3$vost, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
abline(h = 0.5, lwd = 5)
mtext("b", side = 3, line = -1, las = 0, adj = adj, padj = padj,
      cex = par()$cex.lab, font = 2)

plotRingOccurrences(dat.N5$vost,
                    xlab = "Distance from target site (km)", ylab = "", yaxt = "n")
mtext("d", side = 3, line = -1, las = 0, adj = adj, padj = padj,
      cex = par()$cex.lab, font = 2)

mtext("Rank", side = 2, line = 3.5, cex = par()$cex.lab, las = 0,
      outer = TRUE, at = 0.5)
mtext("N = 3", side = 2, line = 5.5, cex = par()$cex.lab, las = 0,
      outer = TRUE, at = 0.75)
mtext("N = 5", side = 2, line = 5.5, cex = par()$cex.lab, las = 0,
      outer = TRUE, at = 0.25)

mtext("EDML", side = 3, line = 0.25, cex = par()$cex.lab,
      outer = TRUE, at = 0.03)
mtext("Vostok", side = 3, line = 0.25, cex = par()$cex.lab,
      outer = TRUE, at = 0.54)

par(op)
dev.off()
