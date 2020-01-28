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
# Determine DML correlation structure for t2m, t2m.pw, oxy and oxy.pw

target.field <- model$t2m

t2m <- dml %>%
  analyseTargetRegion(target.field, study.field = model$lnd.t2m) %>%
  processRegionalMean()

t2m.pw <- dml %>%
  analyseTargetRegion(target.field, study.field = model$lnd.t2m.pw) %>%
  processRegionalMean()

oxy <- dml %>%
  analyseTargetRegion(target.field, study.field = model$lnd.oxy) %>%
  processRegionalMean()

oxy.pw <- dml %>%
  analyseTargetRegion(target.field, study.field = model$lnd.oxy.pw) %>%
  processRegionalMean()

# ------------------------------------------------------------------------------
# Plotting

label <- c(expression(italic("T")["2m"]),
           expression(delta^{18} * "O"),
           expression(delta^{18} * "O"^{"(pw)"}))

col1 <- "black"
col2 <- "green4"
col3 <- "dodgerblue3"

Quartz(file.path(SAVEPATH, "main", "fig_03.pdf"))
par(LoadGraphicsPar())

plot(t2m$bins, t2m$samples$cor, type = "n",
     xlab = "", ylab = "",
     xlim = c(0, 2125), ylim = c(0, 1))
mtext("Distance (km)", side = 1, line = 3.5, cex = par()$cex.lab)
mtext("Correlation", side = 2, line = 3.5, cex = par()$cex.lab, las = 0)

lines(t2m$bins, t2m$samples$cor, col = col1, type = "b", lwd = 2.5)
lines(oxy$bins, oxy$samples$cor, col = col2, type = "b", lwd = 2.5)
lines(oxy.pw$bins, oxy.pw$samples$cor, col = col3, type = "b", lwd = 2.5)

legend("topright", label, col = c(col1, col2, col3),
       lty = 1, lwd = 2.5, bty = "n")

dev.off()

