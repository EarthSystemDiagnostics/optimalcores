##
## aim:
## script to test how many Monte Carlo iterations are needed for sampling two
## cores from rings in order to obtain similar results as with the exact
## solution. In addition, we test the convergence of the optimal correlation
## value for N = 3 depending on the number of Monte Carlo iterations.
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores;
## Muench, Werner and Laepple (2019), in prep.
##

source("init.R")

# ------------------------------------------------------------------------------
# Select data and target site

model  <- selectData()
target <- setTarget(model$t2m, site = "edml")

# ------------------------------------------------------------------------------
# exact solution for N = 2

doubles.exact <- sampleTwoFromRings(field = model$lnd.oxy.pw,
                                    target = target$dat,
                                    distance.field = target$dist) %>%
  processCores()

# ------------------------------------------------------------------------------
# Monte Carlo solutions for N = 2

nmc <- c(10, 100, 1000, 10000, 100000)

doubles.mc <- list()
for (i in 1 : length(nmc)) {

  print(nmc[i])

  doubles.mc[[i]] <- sampleNFromRings(N = 2, nmc = nmc[i],
                                      field = model$lnd.oxy.pw,
                                      target = target$dat,
                                      distance.field = target$dist,
                                      mc.cores = 2) %>%
    processCores()

}

mismatch <- sapply(doubles.mc, function(x, exact) {
  stattools::rmsd(x$correlation$mean, exact$correlation$mean)},
  exact = doubles.exact)

# ------------------------------------------------------------------------------
# Plotting

op <- grfxtools::Quartz(file.path(
  SAVEPATH, "side-results", "monte_carlo_mismatch_ring_sampling_N=2.pdf"),
  mar = c(5, 6, 0.5, 0.5))
          
plot(nmc, mismatch, type = "b", log = "xy", xlab = "", ylab = "",
     ylim = c(0.0005, 0.02), axes = FALSE)
axis(1, at = c(1e1, 1e2, 1e3, 1e4, 1e5),
     labels = c(expression(10^{1}), expression(10^{2}),
                expression(10^{3}), expression(10^{4}),
                expression(10^{5})))
axis(2, at = c(2e-2, 1e-2, 5e-3, 2e-3, 1e-3, 5e-4),
     labels = c(expression(2 %.% 10^{-2}), expression(1 %.% 10^{-2}),
                expression(5 %.% 10^{-3}), expression(2 %.% 10^{-3}),
                expression(1 %.% 10^{-3}), expression(5 %.% 10^{-4})))
mtext("Number of Monte Carlo iterations", 1, 3.5, cex = par()$cex.lab)
mtext("Correlation mismatch (RMSD)", 2, 4.5, cex = par()$cex.lab, las = 0)

par(op)
dev.off()

# ------------------------------------------------------------------------------
# Monte Carlo convergence for N = 3 (w/ optimal ring correlation value for EDML)

nmc <- c(10, 50, 100, 500, 1000, 5000, 10000, 50000, 1e5, 5e5)

triples.mc <- numeric(length = length(nmc))
for (i in 1 : length(nmc)) {

  print(nmc[i])

  triples.mc[i] <- sampleNFromRings(N = 3, nmc = nmc[i],
                                    field = model$lnd.oxy.pw,
                                    target = target$dat,
                                    distance.field = target$dist,
                                    default.ring.combination =
                                      rbind(c(1, 4, 6)),
                                    .parallel = FALSE) %>%
    processCores(n.optim = 1) %>%
    .$optimal.rings %>%
    .$correlation %>%
    .$mean

}

# ------------------------------------------------------------------------------
# Plotting

op <- grfxtools::Quartz(file.path(
  SAVEPATH, "side-results", "monte_carlo_mismatch_ring_sampling_N=3.pdf"),
  mar = c(5, 6, 0.5, 0.5))

plot(nmc, triples.mc - triples.mc[length(triples.mc)],
     type = "b", log = "x", xlab = "", ylab = "",
     xlim = c(1e1, 1e6), ylim = c(-0.02, 0.02), axes = FALSE)
axis(1, at = c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6),
     labels = c(expression(10^{1}), expression(10^{2}),
                expression(10^{3}), expression(10^{4}),
                expression(10^{5}), expression(10^{6})))
axis(2)
abline(h = 0, lty = 2, col = "darkgrey")
mtext("Number of Monte Carlo iterations", 1, 3.5, cex = par()$cex.lab)
mtext("Correlation difference", 2, 4.5, cex = par()$cex.lab, las = 0)

par(op)
dev.off()
