##
## aim:
## script to plot conceptual figure 1
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

source("init.R")

# ------------------------------------------------------------------------------
# Settings

# Set some target site
lat0  <- -82
lon0  <- 107.5
radii <- seq(250, 1000, 250)

# Calculate ring coordinates
target <- data.frame(lat = lat0, lon = lon0)
df.lst <- list()
for (i in 1 : length(radii)) {

  tmp <- geostools::CircleCoordinates(lat0, lon0, radii[i])
  tmp$id <- i
  df.lst[[i]] <- tmp
  
}
rings <- do.call(rbind, df.lst)

# Main study sites and regions
edml <- data.frame(lon = 0, lat = -75)
vost <- data.frame(lon = 106.83, lat = -78.47)

edml.region <- data.frame(
  lon = edml$lon + 17.5 * c(-1, 1, 1, -1),
  lat = edml$lat + 5 * c(1, 1, -1, -1)
)
vost.region <- data.frame(
  lon = vost$lon + 17.5 * c(-1, 1, 1, -1),
  lat = vost$lat + 5 * c(1, 1, -1, -1)
)

# ------------------------------------------------------------------------------
# Create panel (a): Antarctica figure

p1 <- ggplot()

p1 <- p1 +

  geom_point(data = edml, aes(x = lon, y = lat),
             col = "black", bg = "grey", size = 2, pch = 24, stroke = 1) +
  geom_point(data = vost, aes(x = lon, y = lat),
             col = "black", bg = "grey", size = 2, pch = 25, stroke = 1) +
  geom_polygon(data = edml.region, aes(x = lon, y = lat),
               col = "black", fill = "transparent", size = 0.75) +
  geom_polygon(data = vost.region, aes(x = lon, y = lat),
               col = "black", fill = "transparent", size = 0.75) +
  geom_path(data = rings, aes(x = lon, y = lat, group = id),
            col = "red", size = 0.75) +
  geom_point(data = target, aes(x = lon, y = lat),
             col = "black", size = 2.5, pch = 3, stroke = 1)

p1 <- grfxtools::ggpolar(pole = "S", max.lat = -60, min.lat = -90,
                        n.lat.labels = 3,
                        longitude.spacing = 45,
                        land.fill.colour = "transparent",
                        size.outer = 0.5,
                        plt.lat.axes = FALSE,
                        ax.labs.size = 4.75,
                        data.layer = p1)

# ------------------------------------------------------------------------------
# Create panel (b): Conceptual sketch of ring sampling approach

# wrapper function to create circle coordinates
drawCircle <- function(r0, r) {

  phi <- seq(0, 360, 0.1)

  data.frame(
    x = r0 + r * cos(phi),
    y = r0 + r * sin(phi)
  )
}
# wrapper function to create coordinates of selected grid cell combinations
specifyComb <- function(x, y, dx, pnts) {

  data.frame(
    x = mean(pnts) + x * dx,
    y = mean(pnts) + y * dx
  )
}

# create grid
dx <- 2
pnts <- seq(0, 24, dx)

grid <- data.frame(
  x = rep(pnts, times = length(pnts)),
  y = rep(pnts, each = length(pnts))
)

# set a target site
target <- data.frame(x = mean(pnts), y = mean(pnts))

# create selected grid cell combinations
c11 <- specifyComb(c(-1, 1), c(0, 0), dx, pnts)
c12 <- specifyComb(c(0, -1), c(-1, 1), dx, pnts)
c21 <- specifyComb(c(0, 2), c(1, 2), dx, pnts)
c22 <- specifyComb(c(-1, -2), c(-1, 2), dx, pnts)
c31 <- specifyComb(c(0, -3), c(-5, -3), dx, pnts)
c32 <- specifyComb(c(-4, 4), c(-4, 1), dx, pnts)

# create rings
r1 <- drawCircle(mean(pnts), r = 3)
r2 <- drawCircle(mean(pnts), r = 6)
r3 <- drawCircle(mean(pnts), r = 9)
r4 <- drawCircle(mean(pnts), r = 12)

# make ggplot
p2 <- ggplot()

p2 <- p2 +

  theme(line = element_blank(), title = element_blank(),
        text = element_blank(), panel.background = element_blank()) +

  geom_point(aes(x = x, y = y), grid, col = "darkgrey",
             pch = 19, size = 2.5, stroke = 1) +

  geom_point(aes(x = x, y = y), target, col = "black",
             pch = 3, size = 5, stroke = 1.25) +

  geom_path(aes(x = x, y = y), r1, col = "red") +
  geom_path(aes(x = x, y = y), r2, col = "red") +
  geom_path(aes(x = x, y = y), r3, col = "red") +
  geom_path(aes(x = x, y = y), r4, col = "red") +

  geom_point(aes(x = x, y = y), c11, col = "black",
             pch = 1, size = 3, stroke = 1.5) +
  geom_point(aes(x = x, y = y), c12, col = "black",
             pch = 19, size = 3, stroke = 1) +

  geom_point(aes(x = x, y = y), c21, col = "dodgerblue",
             pch = 1, size = 3, stroke = 1.5) +
  geom_point(aes(x = x, y = y), c22, col = "dodgerblue",
             pch = 19, size = 3, stroke = 1) +

  geom_point(aes(x = x, y = y), c31, col = "#d95f02",
             pch = 1, size = 3, stroke = 1.5) +
  geom_point(aes(x = x, y = y), c32, col = "#d95f02",
             pch = 19, size = 3, stroke = 1)

# ------------------------------------------------------------------------------
# Plot entire figure

grfxtools::Quartz(file.path(SAVEPATH, "main", "fig_01.pdf"),
                  height = 6.25, width = 12)

egg::ggarrange(plots = list(p1, p2), nrow = 1, ncol = 2, labels = c("a", "b"),
               label.args = list(gp = grid::gpar(cex = 1.5, font = 2)))

dev.off()
