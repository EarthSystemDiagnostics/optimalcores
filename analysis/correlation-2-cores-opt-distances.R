##
## aim:
## analyse pairwise ring combinations and find the optimal distance of the
## second core when the first core is fixed to the central ring; do this for all
## Antarctic sites.
##
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

source("init.R")

# ------------------------------------------------------------------------------
# Select data

model  <- selectData()
field  <- model$lnd.oxy.pw

continental.latitudes  <- GetLat(field)
continental.longitudes <- GetLon(field)

# ------------------------------------------------------------------------------
# Get optimal distance of second core for all Antarctic sites

# Antarctic grid cells to study:
# exclude three sites which exhibit empty rings causing the code to crash
sites <- (1 : ncol(field))[-c(1, 2, 32)]
n <- length(sites)

optimal.distances <- numeric(n)
for (i in 1 : n) {

  if ((i %% 50) == 0) message(sprintf("%s sites.", i))

  target <- setTarget(model$t2m, site = NULL,
                      lat0 = continental.latitudes[sites[i]],
                      lon0 = continental.longitudes[sites[i]])

  tmp <- sampleTwoFromRings(field = model$lnd.oxy.pw, target = target$dat,
                            distance.field = target$dist,
                            .fix.central = TRUE) %>%
    processCores(n.optim = 1)

  optimal.distances[i] <-
    tmp$input$ring.distances[tmp$optimal.rings$combinations][2]

}

# ------------------------------------------------------------------------------
# Plotting

# as histogram

Quartz(file.path(
  SAVEPATH, "side-results", "ring_sampling_N=2_fix_central_distance_hist.pdf"))
op <- par(LoadGraphicsPar())

phist <- hist(optimal.distances, breaks = seq(0, 2000, 250),
              main = "", xlab = "Inner ring distance of second core (km)",
              ylab = "Counts", col = "darkgrey", border = "dimgrey")

dev.off()
par(op)

# how many cores are placed at least into the second ring?
length(which(optimal.distances > 125)) / n

# as map

ring.dist <- tmp$input$ring.distances
col.scale <- RColorBrewer::brewer.pal(length(ring.dist), "Reds")
names(col.scale) <- ring.dist

df <- data.frame(lon = continental.longitudes[sites],
                 lat = continental.latitudes[sites],
                 dat = factor(as.character(optimal.distances),
                              ordered = TRUE, levels = ring.dist))

Quartz(file.path(
  SAVEPATH, "side-results", "ring_sampling_N=2_fix_central_distance_map.pdf"),
  height = 6, width = 6)

p <- ggplot()

p <- p +

    geom_tile(aes(x = lon, y = lat, fill = dat, colour = dat),
                 data = df, colour = "transparent") +

    scale_fill_manual(values = col.scale,
                      na.value = "transparent",
                      name = "d (km)") +

    theme(legend.key.height = unit(0.75, units = "inches"),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18),
          text = element_text(size = 15))

p <- ecustools::ggpolar(pole = "S", max.lat = -60, min.lat = -90,
                        n.lat.labels = 3,
                        longitude.spacing = 45,
                        land.fill.colour = "transparent",
                        size.outer = 0.5,
                        lat.ax.labs.pos = 180, ax.labs.size = 4.5,
                        data.layer = p)

print(p)

dev.off()

