##
## aim:
## script to randomly pick cores (grid cells) from a climate model field to find
## the optimal sites to reconstruct a target time series.
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

source("init.R")

library(magrittr)
library(ggplot2)

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

Quartz(file.path(SAVEPATH, "side-results",
                 "monte_carlo_picking_stability_N=3.pdf"))
par(LoadGraphicsPar())

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
                         nmc = nmc, return.all = FALSE)
vostok.picks <- pickNCores(N = N, target = "vostok",
                           target.field = model$t2m,
                           study.field = model$lnd.oxy.pw,
                           nmc = nmc, return.all = FALSE)
domef.picks  <- pickNCores(N = N, target = "domef",
                           target.field = model$t2m,
                           study.field = model$lnd.oxy.pw,
                           nmc = nmc, return.all = FALSE)
wdc.picks    <- pickNCores(N = N, target = "wdc",
                           target.field = model$t2m,
                           study.field = model$lnd.oxy.pw,
                           nmc = nmc, return.all = FALSE)

if (SAVE) {

  saved <- list(
    edml.picks   = edml.picks,
    vostok.picks = vostok.picks
  )
  attr(saved, "version") <- Sys.Date()

  saveRDS(saved, file = "analysis/picking.rds")

}

# ------------------------------------------------------------------------------
# Plotting of EDML and Vostok results

colour.scale <- rev(RColorBrewer::brewer.pal(10, "RdYlBu"))

ggplt <- list()
data <- edml.picks
for (i in 1 : length(N)) {

  guide = FALSE
  if (i == length(N)) guide = TRUE

  ggplt[[i]] <- plotPicking(data, N[i], min.lon = -60, max.lon = 90,
                            colour.scale = colour.scale, guide = guide)
}

j <- length(N)
data <- vostok.picks
for (i in 1 : length(N)) {

  guide = FALSE
 if (i == length(N)) guide = TRUE

  ggplt[[i + j]] <- plotPicking(data, N[i], min.lon = 30, max.lon = 180,
                                colour.scale = colour.scale,
                                name = "", guide = guide)
}


Quartz(height = 9, width = 22)
egg::ggarrange(plots = ggplt, nrow = 2, ncol = 3)

dev.copy2pdf(file = file.path(SAVEPATH, "main", "echam5_mpiom_wiso_picking.pdf"))

# ------------------------------------------------------------------------------
# distance of optimal single core from the target; all Antarctic sites

continental.latitudes  <- GetLat(model$lnd.t2m)
continental.longitudes <- GetLon(model$lnd.t2m)

optimal.distances <- numeric(length = ncol(model$lnd.t2m))

for (i in 1 : ncol(model$lnd.t2m)) {

  print(i)

  lat0 <- continental.latitudes[i]
  lon0 <- continental.longitudes[i]

  site <- pickNCores(N = 1, target = NULL, lat0 = lat0, lon0 = lon0,
                     target.field = model$t2m, study.field = model$lnd.oxy.pw,
                     return.all = FALSE, verbose = FALSE)

  optimal.distances[i] <- GetDistance(lat0, lon0,
                                      site$picking[[1]]$sample$lat,
                                      site$picking[[1]]$sample$lon)

}

Quartz(file.path(SAVEPATH, "side-results",
                 "optimal_picking_distance_single_core_antarctica.pdf"))
op <- par(LoadGraphicsPar())

phist <- hist(optimal.distances, breaks = seq(0, 2000, 50),
              main = "", xlab = "Distance from target (km)", ylab = "Counts",
              col = "darkgrey", border = "dimgrey")

dev.off()
par(op)

sum(phist$counts[c(9 : 12, 17 : 20)]) / sum(phist$counts)
