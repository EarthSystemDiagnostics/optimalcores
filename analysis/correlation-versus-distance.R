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
# wrapper function to assess spatial correlation structure

runSpatialCorrelation <- function(region, target.field, study.field) {

  run <- analyseTargetRegion(region = region, target.field = target.field,
                             study.field = study.field, N = 1)

  prc <- data.frame(
  
    ring.distances = run[[1]]$ring.distances.sampled$core1,

    cor = run %>%
      sapply(function(x) {x$correlation$mean}) %>%
      apply(1, mean)

  )

  return(prc)

}

# ------------------------------------------------------------------------------
# determine DML correlation structure for t2m, t2m.pw, oxy and oxy.pw

t2m    <- runSpatialCorrelation(region = dml, target.field = model$t2m,
                                study.field = model$lnd.t2m)

t2m.pw <- runSpatialCorrelation(region = dml, target.field = model$t2m,
                                study.field = model$lnd.t2m.pw)

oxy    <- runSpatialCorrelation(region = dml, target.field = model$t2m,
                                study.field = model$lnd.oxy)

oxy.pw <- runSpatialCorrelation(region = dml, target.field = model$t2m,
                                study.field = model$lnd.oxy.pw)

# ------------------------------------------------------------------------------
# Plotting

label <- c(expression(italic("T")["2m"]),
           expression(delta^{18} * "O"),
           expression(delta^{18} * "O"^{"(pw)"}))

Quartz(file.path(SAVEPATH, "main", "fig_03.pdf"))
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

