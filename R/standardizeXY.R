## -----------------------------------------------------------------------------
## Centering and scaling
##
## F.S. Kurnaz
## -----------------------------------------------------------------------------

## Center the response; the center is returned as attribute "center".
## `index` restricts the computation of the center to a subset of observations.
center.y <- function(y, index = NULL, centerFun = mean) {
  center <- centerFun(if (is.null(index)) y else y[index])
  y <- y - center
  attr(y, "center") <- center
  y
}

## A column without spread (e.g. a 0/1 dummy with more than half zeros has
## MAD = 0) cannot be scaled: use the standard deviation instead and, if that is
## zero as well, leave the column unscaled.
fixScale <- function(scale, xc) {
  bad <- !is.finite(scale) | scale <= sqrt(.Machine$double.eps)
  if (any(bad)) {
    alt <- sqrt(colSums(xc[, bad, drop = FALSE]^2) / max(1, nrow(xc) - 1))
    alt[!is.finite(alt) | alt <= sqrt(.Machine$double.eps)] <- 1
    scale[bad] <- alt
  }
  scale
}

## Center and scale the columns of `x`; center and scale are returned as
## attributes "center" and "scale". `index` restricts their computation to a
## subset of observations, but all rows of `x` are transformed.
##
## The two cases used by the package are handled directly:
##   centerFun = mean,   scaleFun = sd    classical standardization
##   centerFun = median, scaleFun = mad   robust standardization
standardize.x <- function(x, index = NULL, centerFun = mean, scaleFun = sd) {
  rows   <- if (is.null(index)) seq_len(nrow(x)) else index
  center <- if (identical(centerFun, mean)) colMeans(x[rows, , drop = FALSE])
            else apply(x[rows, , drop = FALSE], 2, centerFun)
  x <- x - rep(center, each = nrow(x))

  xc <- if (is.null(index)) x else x[rows, , drop = FALSE]
  scale <- if (identical(centerFun, mean) && identical(scaleFun, sd)) {
    sqrt(colSums(xc^2) / max(1, nrow(xc) - 1))          # = sd(), as xc is centered
  } else if (identical(centerFun, median) && identical(scaleFun, mad)) {
    apply(xc, 2, mad, center = 0)                        # median already swept out
  } else {
    apply(xc, 2, scaleFun)
  }
  scale <- fixScale(scale, xc)
  x <- x / rep(scale, each = nrow(x))

  attr(x, "center") <- center
  attr(x, "scale")  <- scale
  x
}
