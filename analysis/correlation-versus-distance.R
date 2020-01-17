##
## aim:
## script to analyse distance dependence of the correlation with temperature
## for different model variables and averaged over rings and target sites.
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

source("init.R")

# ------------------------------------------------------------------------------
# Select data and define region of target sites

model  <- selectData()
dml    <- setTargetRegion(field = model$lnd.t2m)

# ------------------------------------------------------------------------------
# Example of sampling 1 : 3 cores for a single target site

target <- setTarget(model$t2m, site = "edml")

singles <- sampleOneFromRings(field = model$lnd.t2m, target = target$dat,
                              distance.field = target$dist) %>%
           processCores()

doubles <- sampleTwoFromRings(field = model$lnd.t2m, target = target$dat,
                              distance.field = target$dist) %>%
           processCores()

triples <- sampleNFromRings(N = 3, field = model$lnd.t2m, target = target$dat,
                            distance.field = target$dist) %>%
           processCores()

# ------------------------------------------------------------------------------
# t2m

study.field <- model$lnd.t2m

run <- analyseTargetRegion(region = dml, target.field = model$t2m,
                           study.field = study.field, N = 1)

t2m <- data.frame(
  
  ring.distances = run[[1]]$ring.distances.sampled$core1,

  cor = run %>%
    sapply(function(x) {x$correlation$mean}) %>%
    apply(1, mean)

)

# ------------------------------------------------------------------------------
# t2m.pw

study.field <- model$lnd.t2m.pw

run <- analyseTargetRegion(region = dml, target.field = model$t2m,
                           study.field = study.field, N = 1)

t2m.pw <- data.frame(
  
  ring.distances = run[[1]]$ring.distances.sampled$core1,

  cor = run %>%
    sapply(function(x) {x$correlation$mean}) %>%
    apply(1, mean)

)

# ------------------------------------------------------------------------------
# oxy

study.field <- model$lnd.oxy

run <- analyseTargetRegion(region = dml, target.field = model$t2m,
                           study.field = study.field, N = 1)

oxy <- data.frame(
  
  ring.distances = run[[1]]$ring.distances.sampled$core1,

  cor = run %>%
    sapply(function(x) {x$correlation$mean}) %>%
    apply(1, mean)

)

# ------------------------------------------------------------------------------
# oxy.pw

study.field <- model$lnd.oxy.pw

run <- analyseTargetRegion(region = dml, target.field = model$t2m,
                           study.field = study.field, N = 1)

oxy.pw <- data.frame(
  
  ring.distances = run[[1]]$ring.distances.sampled$core1,

  cor = run %>%
    sapply(function(x) {x$correlation$mean}) %>%
    apply(1, mean)

)

# ------------------------------------------------------------------------------
# Plotting

label <- c(expression(italic("T")["2m"]),
           expression(delta^{18} * "O"),
           expression(delta^{18} * "O"^{"(pw)"}))

Quartz(file.path(SAVEPATH, "echam5_mpiom_wiso_fig_correlation_distance.pdf"))
par(LoadGraphicsPar())

plot(t2m$ring.distances, t2m$cor, type = "n",
     xlab = "", ylab = "",
     xlim = c(0, 2125), ylim = c(0, 1))
mtext("Distance (km)", side = 1, line = 3.5, cex = par()$cex.lab)
mtext("Correlation", side = 2, line = 3.5, cex = par()$cex.lab, las = 0)

lines(t2m$ring.distances, t2m$cor, col = "black", type = "b", lwd = 2.5)
lines(oxy$ring.distances, oxy$cor, col = "green4", type = "b", lwd = 2.5)
lines(oxy.pw$ring.distances, oxy.pw$cor, col = "dodgerblue3", type = "b", lwd = 2.5)

legend("topright", label,
       col = c("black", "green4", "dodgerblue3"),
       lty = 1, lwd = 2.5, bty = "n")

dev.off()

