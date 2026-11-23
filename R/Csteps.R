## -----------------------------------------------------------------------------
## C-steps
##
## A C-step ("concentration step") fits the elastic net on the current subset
## `indx`, ranks all observations by their fit to this model and returns the `h`
## best observations as the new subset together with the value of the criterion
## on that subset. Repeating C-steps never makes the criterion worse (for the
## gaussian and binomial family), which is what the algorithm relies on.
##
## All three functions share the signature
##     CStep.<family>(x, y, indx, h, hsize, alpha, lambda, scal)
## and return list(object, index, ...). `object` is the criterion: it is
## minimized for gaussian and multinomial and maximized for binomial (see
## `familySpec()`). A non-finite `object` flags a degenerate fit.
## -----------------------------------------------------------------------------


## Elastic-net fit without intercept on already centered / standardized data.
## glmnet can fail on a subset (e.g. no convergence for the only lambda); the
## result is then NULL and the C-step reports a degenerate fit instead of
## stopping the whole run.
enetFit <- function(x, y, alpha, lambda, family = "gaussian") {
  tryCatch(
    glmnet(x, y, family = family, alpha = alpha, lambda = lambda,
           standardize = FALSE, intercept = FALSE),
    error = function(e) NULL)
}


## Data used inside a C-step: standardized with the mean and sd of the current
## subset if `scal` is TRUE, otherwise unchanged. For the gaussian family the
## response is centered as well.
cstepData <- function(x, y, indx, scal, center.response = FALSE) {
  if (!scal) return(list(x = x, y = y))
  list(x = standardize.x(x, index = indx, centerFun = mean, scaleFun = sd),
       y = if (center.response) center.y(y, index = indx, centerFun = mean) else y)
}


## ---- gaussian ---------------------------------------------------------------
CStep.gaus <- function(x, y, indx, h, hsize, alpha, lambda, scal) {
  d <- cstepData(x, y, indx, scal, center.response = TRUE)

  fit   <- enetFit(d$x[indx, , drop = FALSE], d$y[indx], alpha, lambda)
  if (is.null(fit)) return(list(object = Inf, index = indx))
  beta  <- as.matrix(fit$beta)
  resid <- d$y - d$x %*% beta

  indxnew <- order(abs(resid))[seq_len(h)]
  obj     <- objGaus(d$x, d$y, beta, indxnew, alpha, lambda)
  list(object = obj, index = indxnew, residu = resid, beta = beta)
}


## ---- binomial ---------------------------------------------------------------
CStep.binom <- function(x, y, indx, h, hsize, alpha, lambda, scal) {
  d <- cstepData(x, y, indx, scal)

  fit   <- enetFit(d$x[indx, , drop = FALSE], d$y[indx], alpha, lambda, "binomial")
  if (is.null(fit)) return(list(object = -Inf, index = indx))
  beta  <- as.matrix(fit$beta)
  eta   <- d$x %*% beta
  resid <- -(d$y * eta) + log1pexp(eta)      # deviance contributions
  if (all(beta == 0)) {
    return(list(object = -Inf, index = indx, residu = resid, beta = beta))
  }

  # keep the class proportions: the h0 best observations of class 0 and the
  # h - h0 best observations of class 1
  ord     <- order(resid)
  h0      <- floor((sum(y == 0) + 1) * hsize)
  h1      <- h - h0
  index0  <- ord[y[ord] == 0][seq_len(h0)]
  index1  <- ord[y[ord] == 1][seq_len(h1)]
  indxnew <- c(index0, index1)

  obj <- objBinom(d$x, d$y, beta, indxnew, alpha, lambda)
  list(object = obj, index = indxnew, residu = resid, beta = beta)
}


## ---- multinomial ------------------------------------------------------------
## Observations are ranked by a robust distance in the space of the linear
## predictors: within each class the Mahalanobis distances based on the MCD are
## rescaled by their median, and the best `floor(n_l * hsize + 1)` observations
## of every class are kept.
CStep.multinom <- function(x, y, indx, h, hsize, alpha, lambda, scal) {
  k  <- ncol(y)
  hk <- floor(colSums(y) * hsize + 1)
  d  <- cstepData(x, y, indx, scal)

  fit  <- enetFit(d$x[indx, , drop = FALSE], d$y[indx, , drop = FALSE], alpha, lambda, "multinomial")
  if (is.null(fit)) return(list(object = Inf, index = indx))
  beta <- do.call(cbind, lapply(fit$beta, as.matrix))     # p x k
  eta  <- d$x %*% beta                                    # linear predictors, n x k
  md   <- rep(NA_real_, nrow(d$x))

  zj.svd  <- svd(scale(eta, TRUE, FALSE))
  svdrank <- sum(zj.svd$d > 1e-4)
  if (svdrank < k - 1) {
    return(list(object = Inf, index = indx, md = md, beta = beta))
  }

  # The MCD uses random subsamples; seeding it from the subset makes the C-step
  # a deterministic function of its input, whichever worker evaluates it. The
  # score space has at most k dimensions, so 50 subsamples are used instead of
  # the default 500 to save computing time.
  zj <- zj.svd$u[, seq_len(svdrank), drop = FALSE]
  localSeed(subsetSeed(indx), {
    for (j in seq_len(k)) {
      in.j   <- d$y[, j] == 1
      mcd    <- covMcd(zj[in.j, , drop = FALSE], nsamp = 50)
      md.j   <- sqrt(mahalanobis(zj[in.j, , drop = FALSE], mcd$center, mcd$cov))
      md[in.j] <- md.j * sqrt(qchisq(0.5, svdrank)) / median(md.j)
    }
  })

  ord     <- order(md)
  indxnew <- unlist(lapply(seq_len(k), function(l) ord[y[ord, l] == 1][seq_len(hk[l])]))

  obj <- objMultinomEta(eta, d$y, beta, indxnew, alpha, lambda)
  list(object = obj, index = indxnew, md = md, beta = beta)
}
