##
## aim:
## script to analyse local correlation between t2m and d18O.pw
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

source("init.R")

library(ggplot2)

# ------------------------------------------------------------------------------
# Correlation analyses between t2m and precip-weighted d18O

model <- selectData()

correlation.field <- pfields::ApplyFields(model$t2m, model$oxy.pw,
                                          cor, use = "pairwise")
correlation.field.df <- pField2df(correlation.field)

# ------------------------------------------------------------------------------
# Some numbers

# correlation range all Antarctica
range(correlation.field.df$dat, na.rm = TRUE)

# see correlation at EDML and Vostok grid cells
pfields::SelPoint(correlation.field, lat = -75, lon = 0)
pfields::SelPoint(correlation.field, lat = -78.47, lon = 106.83)

# ------------------------------------------------------------------------------
# Plotting

# set colour scale
col.scale <- rev(RColorBrewer::brewer.pal(10, "RdYlBu"))

Quartz(file.path(SAVEPATH, "echam5_mpiom_wiso_fig_03.pdf"),
       height = 6, width = 6)
op <- par(LoadGraphicsPar())

p <- ggplot()

p <- p +
        
    geom_tile(aes(x = lon, y = lat, fill = dat, colour = dat),
                 data = correlation.field.df, colour = "transparent") +

    scale_fill_gradientn(colours = col.scale,
                         na.value = "transparent",
                         limits = c(0, 0.6001), name = "Correlation") +

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
par(op)

