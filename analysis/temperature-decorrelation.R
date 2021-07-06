##
## aim:
## script to analyse temperature decorrelation lengths
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

source("init.R")

# ------------------------------------------------------------------------------
# Estimate decorrelation lengths for t2m

model <- selectData()

t2m.decorrelation <- pfields::Decor.pField(model$t2m, print.progress = 100)
t2m.decorrelation.df <- pfields::pField2df(t2m.decorrelation)

# ------------------------------------------------------------------------------
# Plotting

# set colour scale
col.scale <- grfxtools::ColorPal("OrRd")

# set upper/lower quantiles to constant value for plotting purposes
quantile(t2m.decorrelation.df$dat, na.rm = TRUE, probs = c(0.05, 0.95))

t2m.decorrelation.df$dat[t2m.decorrelation.df$dat <= 800] <- 800
t2m.decorrelation.df$dat[t2m.decorrelation.df$dat >= 2500] <- 2500

p0 <- ggplot()

p0 <- p0 +

    geom_tile(aes(x = lon, y = lat, fill = dat, colour = dat),
                 data = t2m.decorrelation.df, colour = "transparent") +

    scale_fill_gradientn(colours = col.scale,
                         na.value = "transparent",
                         limits = c(500, 2500),
                         name = bquote(tau * " (km)")) +

    theme(legend.key.height = unit(0.75, units = "inches"),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18),
          text = element_text(size = 15))

p0 <- grfxtools::ggpolar(pole = "S", max.lat = -60, min.lat = -90,
                        n.lat.labels = 3,
                        longitude.spacing = 45,
                        land.fill.colour = "transparent",
                        size.outer = 0.5,
                        lat.ax.labs.pos = 180, ax.labs.size = 4.5,
                        data.layer = p0)

# -> results are plotted as hardcopy figure together with the
#    correlation maps in 'local-t2m-oxy-correlation.R'
