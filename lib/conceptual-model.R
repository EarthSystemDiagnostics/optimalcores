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
