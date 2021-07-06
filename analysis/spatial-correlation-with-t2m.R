##
## aim:
## script to analyse spatial correlation between t2m and
## fields of t2m, t2m.pw, d18O and d18O.pw
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

source("init.R")

# ------------------------------------------------------------------------------
# Select data

model <- selectData()

target.field <- model$t2m

# ------------------------------------------------------------------------------
# Spatial correlation between t2m at EDML and the fields

edml <- data.frame(lat = -75, lon = 0)

target.site <- setTarget(target.field, site = NULL,
                         lat0 = edml$lat, lon0 = edml$lon)

correlation.map.edml <- list(

  t2m    = pfields::cor.pTs(point = target.site$dat,
                            field = model$lnd.t2m,
                            use = "pairwise"),
  t2m.pw = pfields::cor.pTs(point = target.site$dat,
                            field = model$lnd.t2m.pw,
                            use = "pairwise"),
  oxy    = pfields::cor.pTs(point = target.site$dat,
                            field = model$lnd.oxy,
                            use = "pairwise"),
  oxy.pw = pfields::cor.pTs(point = target.site$dat,
                            field = model$lnd.oxy.pw,
                            use = "pairwise")
)

# ------------------------------------------------------------------------------
# Spatial correlation between t2m at Vostok and the fields

vost <- data.frame(lat = -78.47, lon = 106.83)

target.site <- setTarget(target.field, site = NULL,
                         lat0 = vost$lat, lon0 = vost$lon)

correlation.map.vost <- list(

  t2m    = pfields::cor.pTs(point = target.site$dat,
                            field = model$lnd.t2m,
                            use = "pairwise"),
  t2m.pw = pfields::cor.pTs(point = target.site$dat,
                            field = model$lnd.t2m.pw,
                            use = "pairwise"),
  oxy    = pfields::cor.pTs(point = target.site$dat,
                            field = model$lnd.oxy,
                            use = "pairwise"),
  oxy.pw = pfields::cor.pTs(point = target.site$dat,
                            field = model$lnd.oxy.pw,
                            use = "pairwise")
)

# ------------------------------------------------------------------------------
# Build plot list

ggplt <- list()

ggplt[[1]] <- plotSpatialT2mCorrelation(
  map = pField2df(correlation.map.edml$t2m),
  target = edml, guide = FALSE)
ggplt[[2]] <- plotSpatialT2mCorrelation(
  map = pField2df(correlation.map.edml$t2m.pw),
  target = edml, guide = TRUE)
ggplt[[3]] <- plotSpatialT2mCorrelation(
  map = pField2df(correlation.map.edml$oxy),
  target = edml, guide = FALSE)
ggplt[[4]] <- plotSpatialT2mCorrelation(
  map = pField2df(correlation.map.edml$oxy.pw),
  target = edml, guide = TRUE, name = "")

ggplt[[5]] <- plotSpatialT2mCorrelation(
  map = pField2df(correlation.map.vost$t2m),
  target = vost, guide = FALSE)
ggplt[[6]] <- plotSpatialT2mCorrelation(
  map = pField2df(correlation.map.vost$t2m.pw),
  target = vost, guide = TRUE, name = "")
ggplt[[7]] <- plotSpatialT2mCorrelation(
  map = pField2df(correlation.map.vost$oxy),
  target = vost, guide = FALSE)
ggplt[[8]] <- plotSpatialT2mCorrelation(
  map = pField2df(correlation.map.vost$oxy.pw),
  target = vost, guide = TRUE, name = "")

# ------------------------------------------------------------------------------
# Save plot

labels <- c(expression("(" * bold("a") * ") " * italic("T")["2m"]),
            expression("(" * bold("b") * ") " * italic("T")["2m"]^{"(pw)"}),
            expression("(" * bold("c") * ") " * delta^{18} * "O"),
            expression("(" * bold("d") * ") " * delta^{18} * "O"^{"(pw)"}),
            expression("(" * bold("e") * ") " * italic("T")["2m"]),
            expression("(" * bold("f") * ") " * italic("T")["2m"]^{"(pw)"}),
            expression("(" * bold("g") * ") " * delta^{18} * "O"),
            expression("(" * bold("h") * ") " * delta^{18} * "O"^{"(pw)"}))

grfxtools::Quartz(height = 24, width = 14)
egg::ggarrange(plots = ggplt, nrow = 4, ncol = 2, labels = labels,
               label.args = list(gp = grid::gpar(cex = 2)))

dev.copy2pdf(file = file.path(SAVEPATH, "main", "fig_03.pdf"))
