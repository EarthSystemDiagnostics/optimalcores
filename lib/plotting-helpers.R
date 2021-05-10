##
## aim:
## helper functions for plotting
## relation:
## https://github.com/EarthSystemDiagnostics/optimalcores
##

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
##'   correlation map (for visual purposes, all lower values are set to this
##'   value in the plotted map).
##' @param cor.max numeric; upper correlation value to use in the colour scale
##'   of the map plot.
##' @param min.lon numeric; minimum longitude to plot.
##' @param max.lon numeric; maximum longitude to plot.
##' @param colour.scale vector of colours to use for the correlation map.
##' @param name name for the colour bar legend; defaults to "Correlation".
##' @param guide logical; should the colour bar legend be plotted?
##' @param plotCircle logical; should the circle be plotted from within which
##'   the sites could be picked?
##' @return a ggplot2 object.
##' @author Thomas Münch
plotPicking <- function(data, N, cor.min = 0, cor.max = 0.5,
                        min.lon, max.lon, colour.scale,
                        name = "Correlation", guide = TRUE,
                        plotPickingCircle = FALSE) {

  getPickingCircle <- function(target.lat, target.lon, pick.radius,
                               min.lon, max.lon, max.lat = -60) {

    circle <- geostools::CircleCoordinates(lat0 = target.lat, lon0 = target.lon,
                                           radius.circle = pick.radius,
                                           return.pi.interval = TRUE)

    circle$id <- 1

    circle <- circle[which(circle$lon >= min.lon & circle$lon <= max.lon), ]

    if (any(i <- which(circle$lat > max.lat))) {

      circle[which(circle$lon >= mean(circle$lon[i])), "id"] <- 2
      circle <- circle[-i, ]
    }

    return(circle)
  }

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

  if (plotPickingCircle) {

    circle <- getPickingCircle(data$target$lat, data$target$lon, data$radius,
                               min.lon, max.lon)

    p <- p +
      geom_line(data = circle, aes(x = lon, y = lat, group = id),
                col = "black", size = 0.75)
  }

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

  p <- grfxtools::ggpolar(pole = "S", max.lat = -60, min.lat = -90,
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

  op <- grfxtools::Par(mar = c(0, 0, 0, 0), oma = c(5, 6.75, 2, 6.25))
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
##' Produce a plot of the number of cores each ring around a target site has
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

    grfxtools::Polyplot(x = xlim, y1 = rep(i - 0.5, 2), y2 = rep(i + 0.5, 2),
                        col = shading[i], alpha = alpha)
    points(ring.occurrences[i, ], rep(i, ni), pch = pch, cex = cex)
  }

  abline(v = (ring.distances - ring.distances[1])[-1], col = "darkgrey")

}

##' Plot spatial correlation with target site temperature
##'
##' Produce a ggplot2 map plot of Antarctica with a certain model variable's
##' spatial correlation to the temperature at a given target site, incl. the
##' target site location as well as correlation contour lines.
##'
##' @param map a data frame of the three columns \code{"lat"}, \code{"lon"} and
##'   \code{"dat"} with the spatial correlation field.
##' @param target a data frame of one row and the two columns \code{"lat"} and
##'   \code{"lon"} with the coordinates of the target site.
##' @param cor.min numeric; lower correlation value to set as threshold for the
##'   correlation map (for visual purposes, all lower values are set to this
##'   value in the plotted map).
##' @param cor.max numeric; upper correlation value to use in the colour scale
##'   of the map plot.
##' @param binwidth numeric; the bin width for the correlation contour lines.
##' @param colour.scale vector of colours to use for the correlation map.
##' @param guide logical; should the colour bar legend be plotted?
##' @param name name for the colour bar legend; defaults to "Correlation".
##' @return a ggplot2 object.
##' @author Thomas Münch
plotSpatialT2mCorrelation <- function(map, target, cor.min = -0.5, cor.max = 1,
                                      binwidth = 0.1, colour.scale = NULL,
                                      guide = TRUE, name = "Correlation") {

  if (is.null(colour.scale)) {
    colour.scale <- grfxtools::ColorPal("RdYlBu", n.in = 10, rev = TRUE)
  }

  # constrain correlation range for visual purposes
  map$dat[map$dat < cor.min] <- cor.min

  p <- ggplot()

  p <- p +

  geom_tile(aes(x = lon, y = lat, fill = dat),
            data = map, colour = "transparent") +

  geom_contour(aes(x = lon, y = lat, z = dat),
               data = map, colour = "black", binwidth = binwidth) +

  geom_point(data = target, aes(x = lon, y = lat),
             col = "black", size = 4, pch = 3, stroke = 1.5)

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

  p <- grfxtools::ggpolar(pole = "S", max.lat = -60, min.lat = -90,
                          n.lat.labels = 3,
                          longitude.spacing = 45,
                          land.fill.colour = "transparent",
                          size.outer = 0.5,
                          lat.ax.labs.pos = 180, ax.labs.size = 4.5,
                          data.layer = p)

  p

}
