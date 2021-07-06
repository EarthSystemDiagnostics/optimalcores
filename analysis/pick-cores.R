##
## aim:
## script to randomly pick cores (grid cells) from a climate model field to find
## the optimal sites to reconstruct a target time series.
##
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Münch, Werner and Laepple, Clim. Past, 2021
##

source("init.R")

# save the picking data?
SAVE <- TRUE

# ------------------------------------------------------------------------------
# Select data

model  <- selectData()

# ------------------------------------------------------------------------------
# Test number of sufficient Monte Carlo runs

nmc <- c(10, 100, 1e3, 1e4, 5e4, 1e5, 5e5, 1e6)

test.nmc <- list()
for (i in 1 : length(nmc)) {

  print(nmc[i])
  test.nmc[[i]] <- pickNCores(N = 3, target = "edml",
                              target.field = model$t2m,
                              study.field = model$lnd.oxy.pw,
                              nmc = nmc[i], return.all = FALSE,
                              verbose = FALSE)
}

optimal.correlation <- test.nmc %>%
  sapply(function(x){x$picking[[1]]$metric})

labels = c(expression(10^1), expression(10^2), expression(10^3),
           expression(10^4), expression(10^5), expression(10^6))

grfxtools::Quartz(file.path(SAVEPATH, "side-results",
                            "monte_carlo_picking_stability_N=3.pdf"))

plot(nmc, optimal.correlation, type = "b", log = "x", axes = FALSE,
     ylim = c(0.35, 0.5), xlab = "", ylab = "")
mtext("Number of Monte Carlo picks", 1, 3.5, cex = par()$cex.lab)
mtext("Max. correlation (N = 3)", 2, 3.5, cex = par()$cex.lab, las = 0)
axis(1, at = c(10, 100, 1e3, 1e4, 1e5, 1e6),
     labels = labels)
axis(2)

dev.off()

# ------------------------------------------------------------------------------
# Do the picking for EDML, Vostok, Dome F and WDC

N <- c(1, 3, 5)
nmc <- 1e5

edml.picks <- pickNCores(N = N, target = "edml",
                         target.field = model$t2m,
                         study.field = model$lnd.oxy.pw,
                         nmc = nmc, return.all = FALSE,
                         use = "pairwise")
vostok.picks <- pickNCores(N = N, target = "vostok",
                           target.field = model$t2m,
                           study.field = model$lnd.oxy.pw,
                           nmc = nmc, return.all = FALSE,
                           use = "pairwise")
domef.picks  <- pickNCores(N = N, target = "domef",
                           target.field = model$t2m,
                           study.field = model$lnd.oxy.pw,
                           nmc = nmc, return.all = FALSE,
                           use = "pairwise")
wdc.picks    <- pickNCores(N = N, target = "wdc",
                           target.field = model$t2m,
                           study.field = model$lnd.oxy.pw,
                           nmc = nmc, return.all = FALSE,
                           use = "pairwise")

if (SAVE) {

  saved <- list(
    edml.picks   = edml.picks,
    vostok.picks = vostok.picks
  )
  attr(saved, "version") <- Sys.Date()

  saveRDS(saved, file = "analysis/picking.rds")

}

# ------------------------------------------------------------------------------
# Plot EDML and Vostok results (run above code or load "./analysis/picking.rds")

col.scale <- grfxtools::ColorPal("OrRd")

ggplt <- list()
data <- edml.picks
for (i in 1 : length(N)) {

  guide = FALSE
  if (i == length(N)) guide = TRUE

  ggplt[[i]] <- plotPicking(data, N[i], min.lon = -60, max.lon = 105,
                            colour.scale = col.scale, guide = guide,
                            plotPickingCircle = TRUE)
}

j <- length(N)
data <- vostok.picks
for (i in 1 : length(N)) {

  guide = FALSE
 if (i == length(N)) guide = TRUE

  ggplt[[i + j]] <- plotPicking(data, N[i], min.lon = 30, max.lon = 180,
                                colour.scale = col.scale, name = "",
                                guide = guide, plotPickingCircle = TRUE)
}

labels <- c(expression("(" * bold("a") * ") " * "N = 1"),
            expression("(" * bold("b") * ") " * "N = 3"),
            expression("(" * bold("c") * ") " * "N = 5"),
            expression("(" * bold("d") * ") " * "N = 1"),
            expression("(" * bold("e") * ") " * "N = 3"),
            expression("(" * bold("f") * ") " * "N = 5"))

grfxtools::Quartz(height = 10.6364, width = 26)
egg::ggarrange(plots = ggplt, nrow = 2, ncol = 3, labels = labels,
               label.args = list(gp = grid::gpar(cex = 2)))

dev.copy2pdf(file = file.path(SAVEPATH, "main", "echam5_mpiom_wiso_picking.pdf"))

# ------------------------------------------------------------------------------
# Distance of optimal single core from the target; all Antarctic sites

continental.latitudes  <- GetLat(model$lnd.t2m)
continental.longitudes <- GetLon(model$lnd.t2m)

optimal.distances <- numeric(length = ncol(model$lnd.t2m))
n.cells <- list()
for (i in 1 : ncol(model$lnd.t2m)) {

  print(i)

  lat0 <- continental.latitudes[i]
  lon0 <- continental.longitudes[i]

  target.site <- setTarget(model$t2m, site = NULL, lat0 = lat0, lon0 = lon0)

  site <- pickNCores(N = 1, target = NULL, lat0 = lat0, lon0 = lon0,
                     target.field = model$t2m, study.field = model$lnd.oxy.pw,
                     return.all = FALSE, verbose = FALSE)

  optimal.distances[i] <- GetDistance(lat0, lon0,
                                      site$picking[[1]]$sample$lat,
                                      site$picking[[1]]$sample$lon)

  n.cells[[i]] <- getNumberOfCells(distance.field = target.site$dist,
                                   start = 0, end = 1900, binsize = 100)

  if (i == ncol(model$lnd.t2m)) n.cells <- simplify2array(n.cells)

}

weights <- 1 / rowMeans(n.cells)

phist <- hist(optimal.distances, breaks = seq(0, 2000, 100), plot = FALSE)

phist.weighted <- phist
phist.weighted$counts <- (phist.weighted$counts * weights) / sum(weights)


op <- grfxtools::Quartz(
  file.path(SAVEPATH, "side-results",
            "optimal_picking_distance_single_core_antarctica.pdf"),
  height = 6, width = 8.5, mar = c(5, 5, 0.5, 5))

plot(phist, main = "", xlab = "Distance from target (km)", ylab = "Counts",
     col = adjustcolor("black", 0.2), border = "dimgrey")

par(new = TRUE)

plot(phist.weighted, main = "", xlab = "", ylab = "", axes = FALSE,
     ylim = c(0, 2), col = adjustcolor("firebrick4", 0.2), border = "dimgrey")
axis(4, col = "firebrick4", col.axis = "firebrick4")
text(x = 2350, y = 1, labels = "Relative counts", srt = -90, xpd = NA,
     cex = par()$cex.lab * par()$cex, col = "firebrick4")

legend("topleft", c("Unweighted", "Weighted"),
       col = adjustcolor(c("black", "firebrick4"), 0.2),
       lty = 1, lwd = 10, inset = c(0, 0.05), bty = "n")

par(op)
dev.off()


sum(phist$counts[1 : 5]) / sum(phist$counts)
sum(phist$counts[6 : 15]) / sum(phist$counts)

sum(phist.weighted$counts[1 : 5]) / sum(phist.weighted$counts)
sum(phist.weighted$counts[6 : 15]) / sum(phist.weighted$counts)

