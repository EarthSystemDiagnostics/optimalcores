##
## aim:
## script to test how many Monte Carlo iterations are needed for sampling two
## cores from rings in order to obtain similar results as with the exact
## solution.
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
# exact solution

doubles.exact <- sampleTwoFromRings(field = model$lnd.oxy.pw,
                                    target = target$dat,
                                    distance.field = target$dist) %>%
  processCores()

# ------------------------------------------------------------------------------
# Monte Carlo solutions

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

Quartz(file.path(
  SAVEPATH, "side-results", "monte_carlo_mismatch_ring_sampling_N=2.pdf"))
op <- par(LoadGraphicsPar(mar = c(5, 6, 0.5, 0.5)))
          
plot(nmc, mismatch, type = "b", log = "xy", xlab = "", ylab = "",
     axes = FALSE)
axis(1, at = c(1e1, 1e2, 1e3, 1e4, 1e5),
     labels = c(expression(10^{1}), expression(10^{2}),
                expression(10^{3}), expression(10^{4}),
                expression(10^{5})))
axis(2)
mtext("Number of Monte Carlo iterations", 1, 3.5, cex = par()$cex.lab)
mtext("Correlation mismatch (RMSD)", 2, 4.5, cex = par()$cex.lab, las = 0)

par(op)
dev.off()
