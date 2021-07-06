##
## aim:
## script to estimate parameters of the conceptual model.
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

# ------------------------------------------------------------------------------
# Decorrelation lengths

# t2m

target.field <- model$t2m

dml.t2m <- dml.region %>%
  analyseTargetRegion(target.field, study.field = model$lnd.t2m) %>%
  processRegionalMean()

vos.t2m <- vos.region %>%
  analyseTargetRegion(target.field, study.field = model$lnd.t2m) %>%
  processRegionalMean()

m1 <- lm(log(dml.t2m$correlation$mean) ~ dml.t2m$distances$mean + 0)
m2 <- lm(log(vos.t2m$correlation$mean) ~ vos.t2m$distances$mean + 0)

tau.dml.t2m <- -1 / coefficients(m1)
tau.vos.t2m <- -1 / coefficients(m2)

# t2m.pw

dml.t2m.pw <- dml.region %>%
  analyseTargetRegion(target.field, study.field = model$lnd.t2m.pw) %>%
  processRegionalMean()

m3 <- lm(log(dml.t2m.pw$correlation$mean) ~ dml.t2m.pw$distances$mean)# + 0)

tau.dml.t2m.pw <- -1 / coefficients(m3)[2]
xi.dml.t2m.pw   <- 1 - exp(2 * coefficients(m3)[1])

# oxy.pw

target.field <- model$oxy

dml.oxy <- dml.region %>%
  analyseTargetRegion(target.field, study.field = model$lnd.oxy) %>%
  processRegionalMean()

vos.oxy <- vos.region %>%
  analyseTargetRegion(target.field, study.field = model$lnd.oxy) %>%
  processRegionalMean()

m4 <- lm(log(dml.oxy$correlation$mean) ~ dml.oxy$distances$mean + 0)
m5 <- lm(log(vos.oxy$correlation$mean) ~ vos.oxy$distances$mean + 0)

tau.dml.oxy <- -1 / coefficients(m4)
tau.vos.oxy <- -1 / coefficients(m5)

# ------------------------------------------------------------------------------
# xi parameter from local t2m-t2m.pw correlation

t2m.vs.t2m.pw.cor <- pfields::ApplyFields(model$t2m, model$t2m.pw, cor)
avg.correlation   <- rowMeans(t2m.vs.t2m.pw.cor, na.rm = TRUE)

dml.region <- setTargetRegion(field = t2m.vs.t2m.pw.cor)
vos.region <- setTargetRegion(field = t2m.vs.t2m.pw.cor,
                              min.lat = -83.5, max.lat = -73.5,
                              min.lon = 89.5, max.lon = 124.5)

avg.correlation.dml <- rowMeans(t2m.vs.t2m.pw.cor[, dml.region$field.indices],
                                na.rm = TRUE)
avg.correlation.vos <- rowMeans(t2m.vs.t2m.pw.cor[, vos.region$field.indices])

# xi parameter averaged across Antarctica
xi <- 1 - (avg.correlation)^2

# xi parameter averaged across DML and Vostok region, respectively
xi.dml <- 1 - (avg.correlation.dml)^2
xi.vos <- 1 - (avg.correlation.vos)^2

# ------------------------------------------------------------------------------
# Plot comparison of DML temperature decorrelation lengths

x <- seq(min(dml.t2m$distances$mean), max(dml.t2m$distances$mean), 1)

col1 <- "black"
col2 <- "darkgrey"
col3 <- "green4"
col4 <- "dodgerblue3"

label <- c(expression(italic("T")["2m"]^{"(pw)"}),
           "independent fit",
           "parameter curve",
           expression(italic("T")["2m"] * " with fit"))

op <- grfxtools::Quartz(file.path(SAVEPATH, "side-results",
                                  "t2m_decorrelation_fits_dml.pdf"))

plot(dml.t2m$distances$mean, dml.t2m$correlation$mean, type = "n",
     xlab = "", ylab = "",
     xlim = c(0, 2125), ylim = c(0, 1))
mtext("Distance (km)", side = 1, line = 3.5, cex = par()$cex.lab)
mtext("Correlation", side = 2, line = 3.5, cex = par()$cex.lab, las = 0)

lines(x, exp(-x/tau.dml.t2m), lty = 2, lwd = 1.5, col = col2)

lines(x, sqrt(1 - xi.dml.t2m.pw) * exp(-x/tau.dml.t2m.pw),
      lty = 2, lwd = 1.5, col = col3)
lines(x, sqrt(1 - xi.dml) * exp(-x/tau.dml.t2m),
      lty = 2, lwd = 1.5, col = col4)

lines(dml.t2m$distances$mean, dml.t2m$correlation$mean,
      col = col2, type = "b", lwd = 2.5)
lines(dml.t2m.pw$distances$mean, dml.t2m.pw$correlation$mean,
      col = col1, type = "b", lwd = 2.5)

legend("topright", label, col = c(col1, col3, col4, col2),
       lty = c(1, 2, 2, 1), lwd = c(2.5, 1.5, 1.5, 2), bty = "n")

par(op)
dev.off()

# ------------------------------------------------------------------------------
# RMSD of t2m.pw fits

x <- dml.t2m.pw$distances$mean

independent.fit <- sqrt(1 - xi.dml.t2m.pw) * exp(-x/tau.dml.t2m.pw)
parameter.fit   <- sqrt(1 - xi.dml) * exp(-x/tau.dml.t2m)

stattools::rmsd(dml.t2m.pw$correlation$mean, independent.fit)
stattools::rmsd(dml.t2m.pw$correlation$mean, parameter.fit)
