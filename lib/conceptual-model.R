##
## aim:
## conceptual model for the decorrelation structure of the temperature and
## precipitation-weighted temperature fields
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
modelCorrelationT2m <- function(r1, r2, alpha, tau = 2000) {

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
modelCorrelationT2mPW <- function(r1, r2, alpha, tau = 2000,
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
##'   temperature variance which is redistributed into white noise through the
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
modelCorrelationOxyPW <- function(r1, r2, alpha, tau.d = 1000,
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
