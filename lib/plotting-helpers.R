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

##' Plot picking results
##'
##' Create a ggplot2 plot with the locations of N picked cores (grid cells)
##' plotted on top of the map of correlations between the target time series and
##' all grid cells of the model data field from which the cores were picked.
##'
##' @param data the output of \code{\link{pickNCores}}.
##' @param N integer; the number of picked cores. Determines which result to
##'   plot from \code{data}.
##' @param cor.min numeric; lower correlation value to set as threshold for the
##'   correlation map (for visual puropses, all lower values are set to this
##'   value in the plotted map).
##' @param cor.max numeric; upper correlation value to use in the colour scale
##'   of the map plot.
##' @param min.lon numeric; minimum longitude to plot.
##' @param max.lon numeric; maximum longitude to plot.
##' @param colour.scale vector of colours to use for the correlation map.
##' @param name name for the colour bar legend; defaults to "Correlation".
##' @param guide logical; should the colour bar legend be plotted?
##' @return a ggplot2 object.
##' @author Thomas Münch
plotPicking <- function(data, N, cor.min = 0, cor.max = 0.5,
                        min.lon, max.lon, colour.scale,
                        name = "Correlation", guide = TRUE) {

  Ncores <- sapply(data$picking, function(x) {x$N})
  i <- which(Ncores == N)

  # correlation field as data frame, and crop
  map <- pField2df(data$correlation.map,
                   lon.min = min.lon, lon.max = max.lon)

  # constrain correlation range for visual purposes
  map$dat[map$dat < cor.min] <- cor.min

  picking <- data$picking[[i]]

  p <- ggplot()

  p <- p +

    geom_tile(aes(x = lon, y = lat, fill = dat),
              data = map, colour = "transparent") +

    geom_point(data = data$target, aes(x = lon, y = lat),
               col = "black", size = 5, pch = 3, stroke = 1.5) +

    geom_point(data = picking$sample, aes(x = lon, y = lat),
               col = "black", size = 4, pch = 19)

  if (guide) {

    p <- p +

      scale_fill_gradientn(colours = colour.scale,
                           na.value = "transparent",
                           limits = c(cor.min, cor.max),
                           name = name)
  } else {

    p <- p +

      scale_fill_gradientn(colours = colour.scale,
                           na.value = "transparent",
                           limits = c(cor.min, cor.max),
                           name = name, guide = guide)
  }

  p <- p +

    theme(legend.key.height = unit(0.75, units = "inches"),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          text = element_text(size = 18))

  p <- ecustools::ggpolar(pole = "S", max.lat = -60, min.lat = -90,
                          n.lat.labels = 3,
                          min.lon = min.lon, max.lon = max.lon, rotate = TRUE,
                          longitude.spacing = 30,
                          land.fill.colour = "transparent",
                          size.outer = 0.5,
                          lat.ax.labs.pos = min.lon - c(5, 10),
                          ax.labs.size = 4.75,
                          data.layer = p)

  p

}

##' Plot ring correlation contours
##'
##' Produce a filled contour plot of the expected correlation with a target time
##' series of the average of two cores sampled from consecutive rings around a
##' target site.
##'
##' @param correlation a square matrix of correlations for sampling two cores
##'   from the (ring) bins given by \code{distances}.
##' @param distances numeric vector of the distances of the sampling bins around
##'   the target site. Its length must match the dimensions of
##'   \code{corelation}.
##' @param color.palette a colour palette function to be used to assign colors
##'   in the plot.
##' @param dx numeric vector of minor axis tick positions; the default
##'   \code{NULL} means to use internal tick positions.
##' @param xlab.pos numeric vector of major axis tick positions and labels.
##' @param zlim correlation limits for the plot. If \code{NULL} use default
##'   values from \code{filled.contour}.
##' @param label character string providing a label for the plot.
##' @param ... further parameters passed on to \code{filled.contour}.
##' @author Thomas Münch
plotCorrelationContours <- function(correlation, distances, color.palette,
                                    dx = NULL, xlab.pos = seq(500, 2000, 500),
                                    zlim = NULL, label = "", ...) {

  if (diff(dim(correlation)) != 0) {
    stop("Correlation input is not a square matrix.")
  }
  if (nrow(correlation) != length(distances)) {
    stop("Length of sampling bins must match matrix dimensions.")
  }

  if (!length(zlim)) zlim <- range(correlation, finite = TRUE)
  if (!length(dx)) dx <- (distances - distances[1])[seq(2, 8, 2)]

  op <- par(LoadGraphicsPar(mar = c(0, 0, 0, 0), oma = c(5, 6.75, 2, 6.25)))
  on.exit(par(op))

  filled.contour(x = distances, y = distances, z = correlation,
                 color.palette = color.palette, zlim = zlim,
                 plot.title = {
                   mtext("Distance of first core (km)", side = 1, line = 3.5,
                         cex = par()$cex.lab);
                   mtext("Distance of second core (km)", side = 2, line = 4,
                         cex = par()$cex.lab, las = 0)},
                 plot.axes = {
                   axis(1, at = xlab.pos);
                   axis(1, at = dx, labels = FALSE);
                   axis(2, at = xlab.pos);
                   axis(2, at = dx, labels = FALSE)},
                 ...)

  op.usr <- par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  text(1.16, 0.5, labels = "Correlation",
       srt = -90, xpd = NA, cex = par()$cex.lab)
  text(0.0105, 1.007, labels = label,
       adj = c(0, 0), xpd = NA, cex = par()$cex.lab)
  par(op.usr)

}

##' Plot ring bin sampling occurrence
##'
##' Produce a plot of the number of cores each ring around a a target site has
##' been sampled in the optimal cases.
##'
##' @param data the output of \code{processCores} called with data from
##'   \code{sampleNFromRings}.
##' @param xlab character string with the title for the x axis.
##' @param ylab character string with the title for the y axis.
##' @param xlim limits for the x axis.
##' @param ylim limits for the y axis; if \code{NULL} (the default) the y limits
##'   are calculated internally.
##' @param xaxt single character which specifies the x axis type; see \code{par}
##'   for details.
##' @param yaxt single character which specifies the y axis type; see \code{par}
##'   for details.
##' @param pch either an integer or a single character specifying the point
##'   symbol to use.
##' @param cex numerical value giving the amount by which the point symbols
##'   should be magnified relative to the default.
##' @param col1 character name of the first colour to use for shading the plot
##'   to differentiate between the sampling ranks.
##' @param col2 character name of the second colour to use for shading the plot
##'   to differentiate between the sampling ranks.
##' @param alpha opacity factor for \code{col1} and \code{col2} within [0,1].
##' @author Thomas Münch
plotRingOccurrences <- function(data,
                                xlab = "Ring distance (km)", ylab = "Rank",
                                xlim = c(0, 2250), ylim = NULL,
                                xaxt = "s", yaxt = "s",
                                pch = 16, cex = 1.5,
                                col1 = "grey", col2 = "lightgrey",
                                alpha = 0.5) {

  ring.distances   <- unique(data$ring.distances.sampled$core1)
  ring.occurrences <- arrangeRingOccurrences(data$optimal.rings$counts,
                                             ring.distances)

  nrank <- nrow(ring.occurrences)

  if (!length(ylim)) ylim <- c(0.5, nrank + 0.5)

  shading <- rep(c(col1, col2), nrank)

  plot(1 : 10, type = "n", axes = FALSE, yaxs = "i", xaxs = "i",
       xlab = "", ylab = "", xlim = xlim, ylim = ylim)
  if (xaxt != "n") axis(1)
  if (yaxt != "n") axis(2, at = 1 : nrank)
  mtext(xlab, 1, 3.5, cex = par()$cex.lab)
  mtext(ylab, 2, 3.5, cex = par()$cex.lab, las = 0)
  box()

  for (i in 1 : nrank) {
    ni <- length(ring.occurrences[i, ])

    ecustools::Polyplot(x = xlim, y1 = rep(i - 0.5, 2), y2 = rep(i + 0.5, 2),
                        col = shading[i], alpha = alpha)
    points(ring.occurrences[i, ], rep(i, ni), pch = pch, cex = cex)
  }

  abline(v = (ring.distances - ring.distances[1])[-1], col = "darkgrey")

}
