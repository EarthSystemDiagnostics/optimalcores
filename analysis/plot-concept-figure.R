##
## aim:
## script to plot conceptual figure 1
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

source("init.R")

library(ggplot2)

# ------------------------------------------------------------------------------
# Settings

# Set some target site
lat0  <- -80
lon0  <- 45
radii <- seq(250, 1000, 250)

# Calculate ring coordinates
focus  <- data.frame(lat = lat0, lon = lon0)
df.lst <- list()
for (i in 1 : length(radii)) {

  tmp <- ecustools::CircleCoordinates(lat0, lon0, radii[i])
  tmp$id <- i
  df.lst[[i]] <- tmp
  
}
rings <- do.call(rbind, df.lst)

# ------------------------------------------------------------------------------
# Plotting

Quartz(file.path(SAVEPATH, "main", "fig_01.pdf"), height = 6, width = 6)
op <- par(LoadGraphicsPar())

p <- ggplot()

p <- p +

  geom_path(data = rings, aes(x = lon, y = lat, group = id),
                     col = "red", size = 0.75) +
  geom_point(data = focus, aes(x = lon, y = lat),
                      col = "black", size = 2.5, pch = 3, stroke = 1)

p <- ecustools::ggpolar(pole = "S", max.lat = -60, min.lat = -90,
                        n.lat.labels = 3,
                        longitude.spacing = 45,
                        land.fill.colour = "transparent",
                        size.outer = 0.5,
                        plt.lat.axes = FALSE,
                        ax.labs.size = 4.75,
                        data.layer = p)

print(p)

dev.off()
par(op)
