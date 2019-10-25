##
## aim:
## initialize the "optimalcores" R function bundle by loading all needed
## packages and library functions.
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores
##
## Thomas Muench, 10/2019
##

# ------------------------------------------------------------------------------
# Set path to your base directory of "optimalcores"

path <- "~/research/awi/comb-i/echam-wiso/optimalcores"


# ------------------------------------------------------------------------------
# Load all required packages

if (!require("ecustools", quietly = TRUE)) {
  stop("Package \"ecustools\" needed for \"optimalcores\".",
       call. = FALSE)
}

if (!require("arrangements", quietly = TRUE)) {
  stop("Package \"arrangements\" needed for \"optimalcores\".",
       call. = FALSE)
}

if (!require("abind", quietly = TRUE)) {
  stop("Package \"abind\" needed for \"optimalcores\".",
       call. = FALSE)
}


# ------------------------------------------------------------------------------
# Define function to source entire directory with .R files

sourceDir <- function (path, trace = TRUE, local = FALSE, ...) {
  cat("Sourcing files...\n")
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if (trace) cat(nm, ":")
    source(file.path(path, nm), local = local, ...)
    if (trace) cat("\n")
  }
}


# ------------------------------------------------------------------------------
# Source the "optimalcores" library directory

sourceDir(file.path(path, "lib"))

