## Common style of the enetLTS plots.
## Colors: regular observations are blue, flagged observations red
## (always together with a different shape, so the distinction never relies on
## color alone); classes use blue, red and green, and are always shown by
## shape as well, so that they remain distinguishable under color vision deficiency.

.enetLTS_pal <- list(
  regular = "#2a78d6",
  outlier = "#e34948",
  line    = "#2a78d6",
  classes = c("#2a78d6", "#e34948", "#1baf7a", "#eda100",
              "#e87ba4", "#008300", "#6250d6", "#eb6834")
)

.enetLTS_theme <- function(legend = TRUE) {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.05)),
                   plot.subtitle = ggplot2::element_text(colour = "gray30"),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = if (legend) "right" else "none")
}

## color scale for regular ("0") versus flagged ("1") observations
.diagColour <- function(lab.val)
  ggplot2::scale_colour_manual(values = c("0" = .enetLTS_pal$regular, "1" = .enetLTS_pal$outlier),
                               name = "Diagnostics", breaks = c("0", "1"), labels = lab.val)

## color scale for the classes of the multinomial family
.classColour <- function(K, name = "Groups") {
  cols <- if (K <= length(.enetLTS_pal$classes)) .enetLTS_pal$classes[seq_len(K)] else grDevices::hcl.colors(K, "Dark 3")
  ggplot2::scale_colour_manual(values = unname(cols), name = name)
}

## Points of the multinomial plots: the class is shown by color AND shape,
## the outlier status by the fill (open = regular, filled = flagged), so that
## no information relies on color alone.
.classShapes <- c(21, 22, 24, 23, 25)
## Zoom: when the flagged observations are far away from the regular ones, the
## plot gets a second panel restricted to the range of the regular observations,
## so that the model and the cut-offs remain visible. zx / zy: zoom in x / y.
.zoomLevels <- c("All observations", "Zoom on the regular observations")
.zoomData <- function(d, zx = TRUE, zy = TRUE, zoom = TRUE, factor = 5) {
  d$panel <- factor(.zoomLevels[1], levels = .zoomLevels)
  r <- d[!d$flag, , drop = FALSE]
  if (!zoom || nrow(r) < 2) return(list(d = d, zoom = FALSE))
  rx <- range(r$X); ry <- range(r$Y)
  far <- (zx && diff(range(d$X)) > factor * diff(rx)) || (zy && diff(range(d$Y)) > factor * diff(ry))
  if (!far) return(list(d = d, zoom = FALSE))
  ## keep everything within one range-width of the regular observations: flagged
  ## observations close to the regular ones stay visible, remote ones are left out
  px <- diff(rx); py <- diff(ry)
  keep <- (!zx | (d$X >= rx[1] - px & d$X <= rx[2] + px)) & (!zy | (d$Y >= ry[1] - py & d$Y <= ry[2] + py))
  z <- d[keep, , drop = FALSE]; z$panel <- factor(.zoomLevels[2], levels = .zoomLevels)
  list(d = rbind(d, z), zoom = TRUE)
}
.zoomFacet <- function(zd) if (zd$zoom) ggplot2::facet_wrap(~ panel, nrow = 1, scales = "free") else NULL

## gaussian: regular (blue x) versus flagged (red triangle) observations
.diagPlot <- function(d, lab.val, under = NULL, zx = TRUE, zy = TRUE, zoom = TRUE) {
  zd <- .zoomData(d, zx, zy, zoom)
  dd <- zd$d; dd$ind <- factor(as.integer(dd$flag), levels = c(0, 1))
  ggplot2::ggplot(dd, ggplot2::aes(x = X, y = Y, colour = ind, shape = ind)) +
    under +
    ggplot2::geom_point(size = 2) +
    .diagColour(lab.val) +
    ggplot2::scale_shape_manual(values = c("0" = 4, "1" = 17), name = "Diagnostics",
                                breaks = c("0", "1"), labels = lab.val) +
    .zoomFacet(zd) +
    .enetLTS_theme()
}

.classPlot <- function(d, K, flagLabel = "flagged as outliers", under = NULL,
                       zx = TRUE, zy = TRUE, zoom = TRUE, curveFun = NULL) {
  ## under: reference layers (lines) drawn below the points
  ## curveFun: a curve y = curveFun(x) drawn below the points in every panel,
  ##           over the x-range of that panel
  zd   <- .zoomData(d, zx, zy, zoom)
  d    <- zd$d
  if (!is.null(curveFun)) {
    cv <- do.call(rbind, lapply(split(d, d$panel, drop = TRUE), function(p) {
      xs <- seq(min(p$X), max(p$X), length.out = 500)
      data.frame(X = xs, P = curveFun(xs), panel = p$panel[1])
    }))
    under <- c(under, list(ggplot2::geom_line(data = cv, ggplot2::aes(x = X, y = P), inherit.aes = FALSE,
                                              linetype = "dotted", colour = "gray30")))
  }
  lev  <- levels(d$cls)
  cols <- if (K <= length(.enetLTS_pal$classes)) .enetLTS_pal$classes[seq_len(K)] else grDevices::hcl.colors(K, "Dark 3")
  cols <- stats::setNames(cols, lev)                 # colors tied to the class labels
  shp  <- stats::setNames(rep_len(.classShapes, K), lev)
  ## second legend for the outlier status: open = regular, filled = flagged
  slab <- c("regular observation", if (grepl("best subset", flagLabel)) "not in best subset" else "outlier")
  d$status <- factor(ifelse(d$flag, slab[2], slab[1]), levels = slab)
  ggplot2::ggplot(d, ggplot2::aes(x = X, y = Y, colour = cls, shape = cls)) +
    under +
    ggplot2::geom_point(data = d[!d$flag, , drop = FALSE], ggplot2::aes(alpha = status),
                        fill = "white", size = 2, stroke = 0.7) +
    ggplot2::geom_point(data = d[d$flag, , drop = FALSE], ggplot2::aes(fill = cls, alpha = status),
                        size = 2.6, stroke = 0.7) +
    ggplot2::scale_colour_manual(values = cols, name = "Class", limits = lev, drop = FALSE) +
    ggplot2::scale_shape_manual(values = shp, name = "Class", limits = lev, drop = FALSE) +
    ggplot2::scale_fill_manual(values = cols, guide = "none", limits = lev, drop = FALSE) +
    ggplot2::scale_alpha_manual(values = c(1, 1), limits = slab, drop = FALSE, name = "Diagnostics",
      guide = ggplot2::guide_legend(override.aes = list(shape = 21, colour = "gray30",
                                                        fill = c("white", "gray30"), size = c(2, 2.6)))) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1, override.aes = list(fill = "white", alpha = 1)),
                    shape  = ggplot2::guide_legend(order = 1)) +
    .zoomFacet(zd) +
    .enetLTS_theme()
}
