##
## aim:
## conceptual model for the decorrelation structure of the temperature and
## precipitation-weighted temperature and isotope fields
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores
##

##' Conceptual temperature correlation
##'
##' Conceptual model for the correlation of a target temperature time series
##' with the average of the temperature time series from two sites at given
##' distances from the target assuming a temperature field that follows an
##' isotropic exponential decorrelation with a constant decorrelation length.
##'
##' Note that the function is vectorised on the input \code{r1}, \code{r2} and
##' \code{alpha}, so ascertain to provide meaningful combinations of these
##' vectors.
##'
##' @param r1 numeric; radial distance in km between target and first site.
##' @param r2 numeric; radial distance in km between target and second site.
##' @param alpha numeric; angle in [0, 2pi] between the \code{r1} and \code{r2}
##'   position vectors.
##' @param tau numeric; assumed decorrelation length of the temperature field in
##'   km.
##' @return correlation between the average temperature time series of sites
##'   \code{r1} and \code{r2} and the target site temperature.
##' @author Thomas Münch
modelT2mCorrelation <- function(r1, r2, alpha, tau = 2000) {

  d  <- sqrt(r1^2 + r2^2 - 2 * r1 * r2 * cos(alpha))
  
  c1 <- exp(-r1 / tau) + exp(-r2 / tau)
  c2 <- 2 * (1 + exp(-d / tau))

  c2 <- sqrt(c2)

  return(c1 / c2)
}

##' Conceptual precipitation-weighted temperature correlation
##'
##' Conceptual model for the correlation of a target temperature time series
##' with the average of the precipitation-weighted temperature time series from
##' two sites at given distances from the target assuming (i) a temperature
##' field that follows an isotropic exponential decorrelation with a constant
##' decorrelation length, (ii) that the effect of precipitation weighting
##' results in the creation of white noise, and (iii) that the effect of
##' precipitation weighting decays isotropically and exponentially with a given
##' length scale.
##'
##' Note that the function is vectorised on the input \code{r1}, \code{r2} and
##' \code{alpha}, so ascertain to provide meaningful combinations of these
##' vectors.
##'
##' @param r1 numeric; radial distance in km between target and first site.
##' @param r2 numeric; radial distance in km between target and second site.
##' @param alpha numeric; angle in [0, 2pi] between the \code{r1} and \code{r2}
##'   position vectors.
##' @param tau numeric; assumed decorrelation length of the temperature field in
##'   km.
##' @param tau.pw numeric; assumed decorrelation length of the effect of
##'   precipitation weighting in km.
##' @param xi numeric value in [0, 1] giving the fraction of the total
##'   temperature variance which is redistributed into white noise through the
##'   effect of precipitation weighting, where '0' means no creation of white
##'   noise and '1' means full conversion into white noise.
##' @return correlation between the average precipitation-weighted temperature
##'   time series of sites \code{r1} and \code{r2} and the target site
##'   temperature.
##' @author Thomas Münch
modelT2mPWCorrelation <- function(r1, r2, alpha, tau = 2000,
                                  tau.pw = 500, xi = 0.7) {

  d <- sqrt(r1^2 + r2^2 - 2 * r1 * r2 * cos(alpha))
  
  c1 <- sqrt(1 - xi) * (exp(-r1 / tau) + exp(-r2 / tau))
  c2 <- 2 * (1 + (1 - xi) * exp(-d / tau) + xi * exp(-d / tau.pw))

  c2 <- sqrt(c2)

  return(c1 / c2)

}

##' Conceptual precipitation-weighted isotope-temperature correlation
##'
##' Conceptual model for the correlation of a target temperature time series
##' with the average of the precipitation-weighted oxygen isotope time series from
##' two sites at given distances from the target assuming (i) an oxygen isotope
##' field that follows an isotropic exponential decorrelation with a constant
##' decorrelation length, (ii) a linear decrease with distance of the
##' correlation between temperature and oxygen isotopes until some critical
##' distance, (iii) that the effect of precipitation weighting results in the
##' creation of white noise, and (iv) that the effect of precipitation weighting
##' decays isotropically and exponentially with a given length scale.
##'
##' Note that the function is vectorised on the input \code{r1}, \code{r2} and
##' \code{alpha}, so ascertain to provide meaningful combinations of these
##' vectors.
##'
##' @param r1 numeric; radial distance in km between target and first site.
##' @param r2 numeric; radial distance in km between target and second site.
##' @param alpha numeric; angle in [0, 2pi] between the \code{r1} and \code{r2}
##'   position vectors.
##' @param tau.d numeric; assumed decorrelation length of the oxygen isotope
##'   field in km.
##' @param tau.pw numeric; assumed decorrelation length of the effect of
##'   precipitation weighting in km.
##' @param xi numeric value in [0, 1] giving the fraction of the total
##'   isotope variance which is redistributed into white noise through the
##'   effect of precipitation weighting, where '0' means no creation of white
##'   noise and '1' means full conversion into white noise.
##' @param c0 temperature-isotope correlation at zero distance.
##' @param c1 temperature-isotope correlation at distances larger than the
##'   critical distance.
##' @param d0 critical distance after which the temperature-isotope correlation
##'   stays constant.
##' @return correlation between the average precipitation-weighted oxygen
##'   isotope time series of sites \code{r1} and \code{r2} and the target site
##'   temperature.
##' @author Thomas Münch
modelOxyPWCorrelation <- function(r1, r2, alpha, tau.d = 1000,
                                  tau.pw = 500, xi = 0.7,
                                  c0 = 0.3, c1 = 0.175, d0 = 1500) {

  d <- sqrt(r1^2 + r2^2 - 2 * r1 * r2 * cos(alpha))

  c1 <- sqrt(1 - xi) *
    sum(t2m_oxy_correlation(c(r1, r2), c0 = c0, c1 = c1, d0 = d0))

  c2 <- 2 * (1 + (1 - xi) * exp(-d / tau.d) + xi * exp(-d / tau.pw))

  c2 <- sqrt(c2)

  return(c1 / c2)

}

##' Temperature-isotope correlation
##'
##' This function provides a simple model for the average correlation between
##' temperature and oxygen isotope composition as a function of distance. The
##' model assumes a linear decay in correlation till a critical distance after
##' which the correlation stays constant.
##'
##' @param d vector of distances between temperature and isotope records.
##' @param c0 temperature-isotope correlation at zero distance.
##' @param c1 temperature-isotope correlation at distances larger than the
##'   critical distance.
##' @param d0 critical distance after which the temperature-isotope correlation
##'   stays constant.
##' @param gamma the slope of the linear correlation decay; per default
##'   calculated from the values for \code{c0}, \code{c1} and \code{d0}.
##' @return the temperature-isotope correlation value at distance \code{d}.
##' @author Thomas Münch
t2m_oxy_correlation <- function(d, c0 = 0.3, c1 = 0.175, d0 = 1500,
                                gamma = (c0 - c1) / d0) {

  ifelse(d <= d0, c0 - gamma * d, c1)

}

##' Run conceptual model
##'
##' This is a wrapper function to run the conceptual model for a given climate
##' field to obtain the expected correlation between the target site temperature
##' and the field variable averaged across combinations of two sites being
##' distributed along concentric rings around the target.
##'
##' @param r numeric vector of ring radii to analyse; the two ring sites are
##'   either placed on the same ring or on two different rings.
##' @param alpha numeric vector of polar angle increments between the two sites;
##'   defaults to analysing the full ring circumference in increments of 10
##'   degree.
##' @param tau numeric; assumed decorrelation length of the temperature field in
##'   km (only needed for fields "t2m" and "t2m.pw").
##' @param tau.d numeric; assumed decorrelation length of the oxygen isotope
##'   field in km (only needed for field "oxy.pw").
##' @param tau.pw numeric; assumed decorrelation length of the effect of
##'   precipitation weighting in km (only needed for fields "t2m.pw" and
##'   "oxy.pw").
##' @param xi numeric value in [0, 1] giving the fraction of the total variance
##'   which is redistributed into white noise through the effect of
##'   precipitation weighting, where '0' means no creation of white noise and
##'   '1' means full conversion into white noise (only needed for fields
##'   "t2m.pw" and "oxy.pw").
##' @param c0 temperature-isotope correlation for single sites at zero distance
##'   (only needed for field "oxy.pw").
##' @param c1 temperature-isotope correlation for single sites at distances
##'   larger than the critical distance (only needed for field "oxy.pw").
##' @param d0 critical distance after which the temperature-isotope correlation
##'   for single sites stays constant (only needed for field "oxy.pw").
##' @param field character flag to signal which climate field to analyse;
##'   possible options are "t2m", "t2m.pw" and "oxy.pw".
##' @return a symmetric \code{n x n} matrix with the modelled correlation
##'   between the average across two climate field time series sampled from
##'   rings and the target site temperature, where 'n' is the number of ring
##'   radii \code{r}.
##' @author Thomas Münch
runConceptualModel <- function(r = seq(0, 2000, 10),
                               alpha = seq(0, 350, 10) * (pi / 180),
                               tau = 1900, tau.d = 1000, tau.pw = 500, xi = 0.7,
                               c0 = 0.3, c1 = 0.175, d0 = 1500,
                               field = "t2m") {

  if (field == "t2m") {
    FUN <- function(r1, r2, alpha, tau, ...) {
      modelT2mCorrelation(r1, r2, alpha, tau)}
  } else if (field == "t2m.pw") {
    FUN <- function(r1, r2, alpha, tau, ...) {
      modelT2mPWCorrelation(r1, r2, alpha, tau, tau.pw, xi)}
  } else if (field == "oxy.pw") {
    FUN <- function(r1, r2, alpha, tau, ...) {
      modelOxyPWCorrelation(r1, r2, alpha, ...)}
  } else {
    stop("Unknown field request.")
  }

  n <- length(r)
  m <- length(alpha)

  model <- array(dim = c(m, n, n))
  for (i in 1 : n) {
    for (j in 1 : n) {

      model[, i, j] <- FUN(r1 = r[i], r2 = r[j], alpha = alpha,
                           tau = tau, tau.pw = tau.pw, xi = xi,
                           tau.d = tau.d, c0 = c0, c1 = c1, d0 = d0)

    }
  }

  apply(model, c(2, 3), mean)

}
