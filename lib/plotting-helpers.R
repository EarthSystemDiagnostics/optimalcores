##
## aim:
## helper functions for plotting
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores
##

#' Open quartz device
#'
#' Wrapper to open a quartz device with default dimensions for on-screen
#' plotting or saving to a file.
#' @param file path to a file for storing a hardcopy of the plot including a
#'   supported file extension to set the \code{type} of output (e.g. ".pdf" or
#'   ".png"); i.e. the function extracts the extension from \code{file} and
#'   uses it as the \code{type} argument for the call to
#'   \code{\link{quartz}}. Defaults to \code{NULL} for on-screen plotting.
#' @param type the type of output to use. Defaults to \code{"native"} for
#'   on-screen plotting. If \code{file} is not \code{NULL}, \code{type} is
#'   determined from its extension.
#' @param height the height of the plotting area in inches.  Default ‘6’.
#' @param width the width of the plotting area in inches.  Default ‘8’.
#' @param ... further arguments passed on to \code{\link{quartz}}.
#' @seealso \code{\link{quartz}}
#' @examples
#'
#' # Create an empty on-screen quartz device
#' Quartz()
#'
#' # Store empty plot in pdf format in local directory
#' \dontrun{
#' Quartz(file = file.path(getwd(), "test-quartz.pdf"))
#' dev.off()
#' }
#' @author Thomas Münch
Quartz <- function(file = NULL, type = "native",
                   height = 6, width = 8, ...) {

    # Determine file type from file extension
    if (!is.null(file)) {

        type <- tools::file_ext(file)
        if (nchar(type) == 0)
            stop("No file extension found for setting 'type'.")
        
    }

    # Open device
    quartz(height = height, width = width, file = file, type = type, ...)

}

#' Load graphical parameters
#'
#' This function returns a list of the graphical parameters specified as its
#' function arguments, which can then be set via a call to
#' \code{par}. Change default parameters by passing them in \code{tag = value}
#' form to the function; additional parameters can be specified via
#' \code{...}. This wrapper function provides a convenient way to set a bunch
#' of default together with new graphical parameters and save their old values
#' for later restoring at the same time; see the example. See \code{?par} for
#' information on the individual parameters.
#' @return A list of graphical parameters to be used with \code{par()};
#'   i.e. per default:
#'   \itemize{
#'     \item mar = c(5, 5, 0.5, 0.5)
#'     \item lwd = 1
#'     \item las = 1
#'     \item font.lab = 1
#'     \item font.axis = 1
#'     \item cex.main = 1.5
#'     \item cex.lab = 1.5
#'     \item cex.axis = 1.25
#' }
#' @author Thomas Münch
#' @seealso \code{\link{par}}
#' @examples
#' op <- par(LoadGraphicsPar())
#' plot(1 : 10, main = "Example plot",
#'      xlab = "X title", ylab = "Y title", type = "l")
#' par(op)
LoadGraphicsPar <- function(mar = c(5, 5, 0.5, 0.5), lwd = 1, las = 1,
                            font.lab = 1, font.axis = 1, cex.main = 1.5,
                            cex.lab = 1.5, cex.axis = 1.25, ...) {

    par <- c(as.list(environment()), list(...))
    return(par)

}
