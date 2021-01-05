##
## aim:
## initialize the "optimalcores" R function bundle by loading all needed
## packages and library functions.
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores
##

# ------------------------------------------------------------------------------
# Load required packages

required.packages <- c("abind", "arrangements", "egg", "geostools", "grfxtools",
                       "magrittr", "pfields", "RColorBrewer", "stattools")

for (pkg in required.packages) {

  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(sprintf("Package \"%s\" needed for \"optimalcores\". ", pkg),
         "Please install it using \"dependencies.R\".", call. = FALSE)
  }
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

sourceDir("lib")

