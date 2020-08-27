##
## aim:
## script to plot conceptual model results.
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

source("init.R")

# ------------------------------------------------------------------------------
# Model parameters

# ring radii
r <- seq(0, 2000, 10)

# t2m decorrelation length (DML/Vostok region estimate)
tau <- 1900

# oxy decorrelation length (DML/Vostok region estimate)
tau.d <- 1100

# decorrelation length of precipitation intermittency effect (unconstrained)
tau.pw <- 500

# intermittency parameter (DML region estimate)
xi <- 0.73
# intermittency parameter (Vostok region estimate)
xi <- 0.82

# t2m-oxy correlation parameters DML (simple linear decay with distance)
c0 <- 0.4
c1 <- 0.
d0 <- 6000
# t2m-oxy correlation parameters Vostok
c0 <- 0.5
c1 <- 0.
d0 <- 2500

# ------------------------------------------------------------------------------
# Run model for the various fields

t2m    <- runConceptualModel(r, tau = tau, tau.d = tau.d, tau.pw = tau.pw,
                             xi = xi, c0 = c0, c1 = c1, d0 = d0,
                             field = "t2m")

t2m.pw <- runConceptualModel(r, tau = tau, tau.d = tau.d, tau.pw = tau.pw,
                             xi = xi, c0 = c0, c1 = c1, d0 = d0,
                             field = "t2m.pw")

oxy.pw <- runConceptualModel(r, tau = tau, tau.d = tau.d, tau.pw = tau.pw,
                             xi = xi, c0 = c0, c1 = c1, d0 = d0,
                             field = "oxy.pw")

# ------------------------------------------------------------------------------
# Plot model results

color.palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdYlBu")))
distances <- r
label <- c(expression("(" * bold("a") * ") " * italic("T")["2m"]),
           expression("(" * bold("b") * ") " * italic("T")["2m"]^{"(pw)"}),
           expression("(" * bold("c") * ") " * delta^{18} * "O"^{"(pw)"}))

filebase <- "conceptual_model_correlation_dml_"
filebase <- "conceptual_model_correlation_vostok_"

Quartz(dpi = 300, file = file.path(
  SAVEPATH, "main", paste0(filebase, "01_t2m.png")))
plotCorrelationContours(t2m, distances, color.palette, zlim = c(0.2, 1),
                        xlab.pos = seq(0, 2000, 500), label = label[1],
                        dx = seq(250, 1750, 500))
dev.off()

Quartz(dpi = 300, file = file.path(
  SAVEPATH, "main", paste0(filebase, "02_t2m.pw.png")))
plotCorrelationContours(t2m.pw, distances, color.palette, zlim = c(0.2, 0.6),
                        xlab.pos = seq(0, 2000, 500), label = label[2],
                        dx = seq(250, 1750, 500))
dev.off()

Quartz(dpi = 300, file = file.path(
  SAVEPATH, "main", paste0(filebase, "03_oxy.pw.png")))
plotCorrelationContours(oxy.pw, distances, color.palette, zlim = c(0.15, 0.25),
                        xlab.pos = seq(0, 2000, 500), label = label[3],
                        dx = seq(250, 1750, 500),
                        key.axes = axis(4, seq(0.15, 0.25, by = 0.02)))
dev.off()

