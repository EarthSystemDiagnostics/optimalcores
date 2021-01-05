##
## aim:
## script to plot the increase in correlation with the number of cores averaged.
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

source("init.R")

# ------------------------------------------------------------------------------
# Produce or load correlation data

model <- selectData()
target1 <- setTarget(model$t2m, site = "edml")
target2 <- setTarget(model$t2m, site = "vostok")

N <- c(1, 2, 3, 5, 7, 10)

edml <- list(optim = numeric(length(N)), local = numeric(length(N)))
vost <- list(optim = numeric(length(N)), local = numeric(length(N)))

# N = 1

single.edml <- sampleOneFromRings(field = model$lnd.oxy.pw, target = target1$dat,
                                 distance.field = target1$dist) %>%
  processCores(n.optim = 1)

single.vost <- sampleOneFromRings(field = model$lnd.oxy.pw, target = target2$dat,
                                 distance.field = target2$dist) %>%
  processCores(n.optim = 1)

edml$optim[1] <- single.edml$optimal.rings$correlation$mean
edml$local[1] <- single.edml$correlation$mean[1]

vost$optim[1] <- single.vost$optimal.rings$correlation$mean
vost$local[1] <- single.vost$correlation$mean[1]

# N = 2

double.edml <- sampleTwoFromRings(field = model$lnd.oxy.pw, target = target1$dat,
                                 distance.field = target1$dist) %>%
  processCores(n.optim = 1)

double.vost <- sampleTwoFromRings(field = model$lnd.oxy.pw, target = target2$dat,
                                 distance.field = target2$dist) %>%
  processCores(n.optim = 1)

edml$optim[2] <- double.edml$optimal.rings$correlation$mean
edml$local[2] <- double.edml$correlation$mean[1]

vost$optim[2] <- double.vost$optimal.rings$correlation$mean
vost$local[2] <- double.vost$correlation$mean[1]

# N = 3, 5

dat <- readRDS("analysis/ring_occurrences_edml_vostok_N=3_N=5.rds")

edml$optim[3] <- dat$edml$N3$optimal.rings$correlation$mean[1]
edml$optim[4] <- dat$edml$N5$optimal.rings$correlation$mean[1]

edml$local[3] <- dat$edml$N3$correlation$mean[1]
edml$local[4] <- dat$edml$N5$correlation$mean[1]

vost$optim[3] <- dat$vost$N3$optimal.rings$correlation$mean[1]
vost$optim[4] <- dat$vost$N5$optimal.rings$correlation$mean[1]

vost$local[3] <- dat$vost$N3$correlation$mean[1]
vost$local[4] <- dat$vost$N5$correlation$mean[1]

# N = 7

dat <- readRDS("analysis/ring_occurrences_edml_vostok_N=7.rds")

edml$optim[5] <- dat$edml$optimal.rings$correlation$mean
edml$local[5] <- dat$edml$correlation$mean[1]

vost$optim[5] <- dat$vost$optimal.rings$correlation$mean
vost$local[5] <- dat$vost$correlation$mean[1]

# N = 10

dat <- readRDS("analysis/ring_occurrences_edml_vostok_N=10.rds")

edml$optim[6] <- dat$edml$optimal.rings$correlation$mean
edml$local[6] <- dat$edml$correlation$mean[1]

vost$optim[6] <- dat$vost$optimal.rings$correlation$mean
vost$local[6] <- dat$vost$correlation$mean[1]

# ------------------------------------------------------------------------------
# Distribution of single correlations for optimal ring combination (N = 3)

model  <- selectData()
field  <- model$lnd.oxy.pw

target  <- setTarget(model$t2m, site = "edml")
system.time(
n3.edml <- sampleNFromRings(N = 3, nmc = 10^5, field = field,
                            target = target$dat,
                            distance.field = target$dist,
                            default.ring.combination = rbind(c(1, 4, 5)),
                            .parallel = FALSE)
)

target  <- setTarget(model$t2m, site = "vostok")
system.time(
n3.vost <- sampleNFromRings(N = 3, nmc = 10^5, field = field,
                            target = target$dat,
                            distance.field = target$dist,
                            default.ring.combination = rbind(c(2, 3, 3)),
                            .parallel = FALSE)
)

dist1 <- n3.edml$correlations[1, ]
dist2 <- n3.vost$correlations[1, ]

# how many of these are above the local correlation?
length(which(dist1 > edml$local[1])) / length(dist1)
length(which(dist2 > vost$local[1])) / length(dist2)

# ------------------------------------------------------------------------------
# Plotting

col1 = "firebrick4"
col2 = "dodgerblue4"

label1 <- expression(bold("a"))
label2 <- expression(bold("b"))

Quartz(file.path(SAVEPATH, "main", "fig_08.pdf"),
       height = 5, width = 12)
op <- par(LoadGraphicsPar(mfcol = c(1, 2),
                          mar = c(5, 6, 2.5, 1), xaxs = "i", yaxs = "i"))

plot(N, edml$optim, type = "n", xlab = "", ylab = "", log = "x", axes = FALSE,
     xlim = range(N), ylim = c(0.2, 0.5))
axis(1, at = N)
axis(2)
mtext("Number of averaged locations", side = 1, line = 3.5, cex = par()$cex.lab)
mtext("Correlation", side = 2, line = 4.5, cex = par()$cex.lab, las = 0)
mtext(label1, side = 3, line = -1.5, cex = par()$cex.lab,
      adj = 0.02, padj = -0.35)

grfxtools::Polyplot(N, y1 = edml$local, y2 = edml$optim, col = col1)
grfxtools::Polyplot(N, y1 = vost$local, y2 = vost$optim, col = col2)

lines(N, edml$local, col = col1, lty = 2, lwd = 1.5)
lines(N, edml$optim, col = col1, lty = 1, lwd = 2.5)

lines(N, vost$local, col = col2, lty = 2, lwd = 1.5)
lines(N, vost$optim, col = col2, lty = 1, lwd = 2.5)

legend("bottomleft", c("EDML site", "Vostok site"), col = c(col1, col2),
       lty = 1, lwd = 2.5, bty = "n")
legend("bottomright", c("Local ring", "Optimal combination"), col = "darkgrey",
       lty = c(2, 1), lwd = 1.5, bty = "n")

hist(dist1, xlab = "", ylab = "", main = "", axes = FALSE,
     xlim = c(0, 0.6), ylim = c(0, 25000), breaks = 24,
     col = adjustcolor(col1, alpha = 0.2))
axis(1)
axis(2, at = seq(0, 25000, 5000), labels = seq(0, 25, 5))
mtext("Correlation", side = 1, line = 3.5, cex = par()$cex.lab)
mtext(expression("Counts ("%*%10^3 * ")"), side = 2, line = 3.5,
      cex = par()$cex.lab, las = 0)
mtext(label2, side = 3, line = -1.5, cex = par()$cex.lab,
      adj = 0.02, padj = -0.35)

lines(x = rep(edml$local[1], 2), y = c(0, 25000), lwd = 2.5, col = col1)

hist(dist2, add = TRUE, col = adjustcolor(col2, alpha = 0.2))
lines(x = rep(vost$local[1], 2), y = c(0, 25000), lwd = 2.5, col = col2)

par(op)
dev.off()
