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

model <- selectData()

dml.region <- setTargetRegion(field = model$lnd.t2m)
vos.region <- setTargetRegion(field = model$lnd.t2m,
                              min.lat = -83.5, max.lat = -73.5,
                              min.lon = 89.5, max.lon = 124.5)

target.field <- model$t2m

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

dml <- list()

dml$t2m <- dml.region %>%
  analyseTargetRegion(target.field, study.field = model$lnd.t2m) %>%
  processRegionalMean()

dml$t2m.pw <- dml.region %>%
  analyseTargetRegion(target.field, study.field = model$lnd.t2m.pw) %>%
  processRegionalMean()

dml$oxy <- dml.region %>%
  analyseTargetRegion(target.field, study.field = model$lnd.oxy) %>%
  processRegionalMean()

dml$oxy.pw <- dml.region %>%
  analyseTargetRegion(target.field, study.field = model$lnd.oxy.pw) %>%
  processRegionalMean()

# ------------------------------------------------------------------------------
# Determine Vostok region correlation structure for t2m, t2m.pw, oxy and oxy.pw

vos <- list()

vos$t2m <- vos.region %>%
  analyseTargetRegion(target.field, study.field = model$lnd.t2m) %>%
  processRegionalMean()

vos$t2m.pw <- vos.region %>%
  analyseTargetRegion(target.field, study.field = model$lnd.t2m.pw) %>%
  processRegionalMean()

vos$oxy <- vos.region %>%
  analyseTargetRegion(target.field, study.field = model$lnd.oxy) %>%
  processRegionalMean()

vos$oxy.pw <- vos.region %>%
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

label1 <- expression(bold("(a) ") * "DML")
label2 <- expression(bold("(b) ") * "Vostok")

Quartz(file.path(SAVEPATH, "main", "fig_03.pdf"),
       height = 5.5, width = 12)
op <- par(LoadGraphicsPar(mfcol = c(1, 2),
                          mar = c(0, 0, 0, 0),
                          oma = c(5, 5, 2, 0.5)))

plot(dml$t2m$distances$mean, dml$t2m$correlation$mean, type = "n",
     xlab = "", ylab = "",
     xlim = c(0, 2125), ylim = c(0, 1))
mtext("Distance (km)", side = 1, line = 3.5, cex = par()$cex.lab)
mtext("Correlation", side = 2, line = 3.5, cex = par()$cex.lab, las = 0)
mtext(label1, side = 3, line = 1.1, cex = par()$cex.lab, adj = 0.02, padj = 1)

lines(dml$t2m$distances$mean, dml$t2m$correlation$mean,
      col = col1, type = "b", lwd = 2.5)
lines(dml$oxy$distances$mean, dml$oxy$correlation$mean,
      col = col2, type = "b", lwd = 2.5)
lines(dml$oxy.pw$distances$mean, dml$oxy.pw$correlation$mean,
      col = col3, type = "b", lwd = 2.5)

plot(vos$t2m$distances$mean, vos$t2m$correlation$mean, type = "n",
     xlab = "", ylab = "", yaxt = "n",
     xlim = c(0, 2125), ylim = c(0, 1))
mtext("Distance (km)", side = 1, line = 3.5, cex = par()$cex.lab)
mtext(label2, side = 3, line = 1.1, cex = par()$cex.lab, adj = 0.02, padj = 1)

lines(vos$t2m$distances$mean, vos$t2m$correlation$mean,
      col = col1, type = "b", lwd = 2.5)
lines(vos$oxy$distances$mean, vos$oxy$correlation$mean,
      col = col2, type = "b", lwd = 2.5)
lines(vos$oxy.pw$distances$mean, vos$oxy.pw$correlation$mean,
      col = col3, type = "b", lwd = 2.5)

legend("topright", label, col = c(col1, col2, col3),
       lty = 1, lwd = 2.5, bty = "n")

dev.off()
par(op)
