## Class-free screening for leverage points
##
## The class-wise robust distances of the reweighting step break down when one
## class consists mostly of leverage points (e.g. contamination that falls
## into a single class): the MCD then describes the contamination, not the
## class. Leverage points are however outlying in the predictors whatever
## their class, so they can be found without the labels: robust
## standardization of the columns, spatial-sign PCA (Locantore et al., 1999),
## and robust (MCD) distances of the first principal component scores.
##
## Returns the logical vector of flagged observations. At most a fraction
## `maxfrac` of the observations is flagged (those with the largest distances).
leverageScreen <- function(x, kmax = 10, varp = 0.9, level = 0.975, maxfrac = 0.25) {
  n <- nrow(x); p <- ncol(x)
  if (n < 20 || p < 2) return(rep(FALSE, n))
  xr <- standardize.x(x, centerFun = median, scaleFun = mad)
  nr <- sqrt(rowSums(xr^2))
  if (any(!is.finite(nr)) || all(nr == 0)) return(rep(FALSE, n))
  U  <- xr / ifelse(nr > 0, nr, 1)
  ev <- svd(U, nu = 0)
  d2 <- ev$d^2
  k  <- min(kmax, which(cumsum(d2) / sum(d2) >= varp)[1], n %/% 3, p)
  if (!is.finite(k) || k < 1) return(rep(FALSE, n))
  V  <- ev$v[, seq_len(k), drop = FALSE]
  tt <- xr %*% V
  sd <- tryCatch({
    mcd <- robustbase::covMcd(tt, alpha = 0.75, nsamp = 100)
    sqrt(mahalanobis(tt, mcd$center, mcd$cov))
  }, error = function(e) NULL)
  if (is.null(sd) || any(!is.finite(sd)) || median(sd) <= 0) return(rep(FALSE, n))
  sd <- sd * sqrt(qchisq(0.5, k)) / median(sd)
  flag <- sd > sqrt(qchisq(level, k))
  maxn <- floor(maxfrac * n)
  if (sum(flag) > maxn) flag <- rank(-sd, ties.method = "first") <= maxn
  flag
}
