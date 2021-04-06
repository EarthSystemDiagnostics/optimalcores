##
## aim:
## script to analyse local correlation between t2m and d18O.pw
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

source("init.R")

# ------------------------------------------------------------------------------
# Correlation analyses between t2m and precip-weighted d18O

model <- selectData()

correlation.field <- pfields::ApplyFields(model$t2m, model$oxy.pw,
                                          cor, use = "pairwise")
correlation.field.pw <- pfields::ApplyFields(model$t2m.pw, model$oxy.pw,
                                             cor, use = "pairwise")
correlation.field.df <- pField2df(correlation.field)
correlation.field.pw.df <- pField2df(correlation.field.pw)

difference <- correlation.field.pw.df
difference$dat <- correlation.field.pw.df$dat - correlation.field.df$dat

# ------------------------------------------------------------------------------
# Some numbers

# mean and range of correlation across Antarctica
range(correlation.field.df$dat, na.rm = TRUE)
mean(correlation.field.df$dat, na.rm = TRUE)

range(correlation.field.pw.df$dat, na.rm = TRUE)
mean(correlation.field.pw.df$dat, na.rm = TRUE)

# correlation percentage below 0.4
n <- length(na.omit(correlation.field.df$dat))

length(which(correlation.field.df$dat < 0.4)) / n

# see correlation at EDML and Vostok grid cells
pfields::SelPoint(correlation.field, lat = -75, lon = 0)
pfields::SelPoint(correlation.field, lat = -78.47, lon = 106.83)

# ------------------------------------------------------------------------------
# Prepare plotting

# set colour scale
col.scale <- RColorBrewer::brewer.pal(9, "OrRd")

# for t2m
p <- ggplot()

p <- p +

    geom_tile(aes(x = lon, y = lat, fill = dat, colour = dat),
                 data = correlation.field.df, colour = "transparent") +

    scale_fill_gradientn(colours = col.scale,
                         na.value = "transparent",
                         limits = c(0, 0.8001), name = "Correlation") +

    theme(legend.key.height = unit(0.75, units = "inches"),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18),
          text = element_text(size = 15))

p1 <- grfxtools::ggpolar(pole = "S", max.lat = -60, min.lat = -90,
                         n.lat.labels = 3,
                         longitude.spacing = 45,
                         land.fill.colour = "transparent",
                         size.outer = 0.5,
                         lat.ax.labs.pos = 180, ax.labs.size = 4.5,
                         data.layer = p)

# for t2m(pw)
correlation.field.pw.df$dat[which(correlation.field.pw.df$dat < 0)] <- 0
p <- ggplot()

p <- p +

    geom_tile(aes(x = lon, y = lat, fill = dat, colour = dat),
                 data = correlation.field.pw.df, colour = "transparent") +

    scale_fill_gradientn(colours = col.scale,
                         na.value = "transparent",
                         limits = c(0, 0.8001), name = "Correlation") +

    theme(legend.key.height = unit(0.75, units = "inches"),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18),
          text = element_text(size = 15))

p2 <- grfxtools::ggpolar(pole = "S", max.lat = -60, min.lat = -90,
                         n.lat.labels = 3,
                         longitude.spacing = 45,
                         land.fill.colour = "transparent",
                         size.outer = 0.5,
                         lat.ax.labs.pos = 180, ax.labs.size = 4.5,
                         data.layer = p)

# for the difference
col.scale <- rev(RColorBrewer::brewer.pal(10, "RdYlBu"))

p <- ggplot()

p <- p +

    geom_tile(aes(x = lon, y = lat, fill = dat, colour = dat),
                 data = difference, colour = "transparent") +

    scale_fill_gradientn(colours = col.scale,
                         na.value = "transparent",
                         limits = c(-0.2, 0.41), name = "Difference") +

    theme(legend.key.height = unit(0.75, units = "inches"),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18),
          text = element_text(size = 15))

p3 <- grfxtools::ggpolar(pole = "S", max.lat = -60, min.lat = -90,
                         n.lat.labels = 3,
                         longitude.spacing = 45,
                         land.fill.colour = "transparent",
                         size.outer = 0.5,
                         lat.ax.labs.pos = 180, ax.labs.size = 4.5,
                         data.layer = p)

# ------------------------------------------------------------------------------
# Combine and save plots

# OBS!
# Load prepared plot of the results from 'temperature-decorrelation.R' as 'p0'
# to include in this figure

Quartz(height = 11, width = 14)
egg::ggarrange(plots = list(p0, p1, p2, p3), nrow = 2, ncol = 2,
               labels = c("a", "b", "c", "d"),
               label.args = list(gp = grid::gpar(font = 2, cex = 1.6)))

dev.copy2pdf(file = file.path(SAVEPATH, "main", "fig_02.pdf"))
