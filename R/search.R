## -----------------------------------------------------------------------------
## Search for outlier-free subsets
##
## The search runs in three phases (Kurnaz, Hoffmann & Filzmoser 2018):
##
##   1. `initialSubsets()`  draws `nsamp` elemental subsets and applies two
##                          C-steps to each of them (first on the elemental
##                          subset, then on h observations);
##   2. `startCsteps()`     keeps the `s1` most promising ones, iterates
##                          C-steps on them until convergence and returns the
##                          best subset for the first (alpha, lambda) pair;
##   3. `warmCsteps.mod()`  walks over the (alpha, lambda) grid, warm-starting
##                          every C-step sequence from the neighboring subset.
##
## Everything that differs between the three model families is collected in
## `familySpec()`.
## -----------------------------------------------------------------------------


## ---- family-specific building blocks ---------------------------------------

## Random elemental subsets: one column per subset.
##   gaussian     3 observations
##   binomial     2 observations from each class
##   multinomial  2 observations from each class
subsets.gaus <- function(x, y, nsamp) {
  replicate(nsamp, sample.int(nrow(x), 3))
}

subsets.binom <- function(x, y, nsamp) {
  replicate(nsamp, c(sample(which(y == 1), 2), sample(which(y == 0), 2)))
}

subsets.multinom <- function(x, y, nsamp) {
  replicate(nsamp, unlist(lapply(seq_len(ncol(y)), function(l) sample(which(y[, l] == 1), 2))))
}


## The penalty parameter is scaled differently in the different phases. These
## conventions are part of the method and are collected here in one place:
##
##                    first C-step   second C-step   iterated C-steps
##                    (elemental)    (h obs.)        (start subsets)
##   gaussian         lambda / 3     lambda / h      lambda / h
##   binomial         lambda / 4     lambda / h      lambda / h
##   multinomial      lambda         lambda          lambda
##
## The C-steps along the (alpha, lambda) grid and the final fit use `lambda`.
familySpec <- function(family) {
  switch(family,
    gaussian = list(
      cstep = CStep.gaus, subsets = subsets.gaus, maximize = FALSE,
      lambda.first = function(lambda, h) lambda / 3,
      lambda.second = function(lambda, h) lambda / h,
      lambda.iter = function(lambda, h) lambda / h,
      lambda.start = function(lambdas) lambdas[1]),
    binomial = list(
      cstep = CStep.binom, subsets = subsets.binom, maximize = TRUE,
      lambda.first = function(lambda, h) lambda / 4,
      lambda.second = function(lambda, h) lambda / h,
      lambda.iter = function(lambda, h) lambda / h,
      lambda.start = function(lambdas) lambdas[1]),
    multinomial = list(
      cstep = CStep.multinom, subsets = subsets.multinom, maximize = FALSE,
      lambda.first = function(lambda, h) lambda,
      lambda.second = function(lambda, h) lambda,
      lambda.iter = function(lambda, h) lambda,
      lambda.start = function(lambdas) lambdas[ceiling(length(lambdas) / 2)]),
    stop("unknown family: ", family))
}


## ---- C-step iteration -------------------------------------------------------

## Iterate C-steps from `index` until the criterion stops improving (at most
## `nCsteps` steps) and return the best subset seen together with its criterion
## value. `maximize` says in which direction the criterion improves.
iterateCsteps <- function(CStep.fun, x, y, index, h, hsize, alpha, lambda, scal,
                          nCsteps, maximize = FALSE) {
  better <- if (maximize) function(new, old) new > old else function(new, old) new < old

  cs        <- CStep.fun(x, y, index, h, hsize, alpha, lambda, scal)
  best.obj  <- cs$object
  best.idx  <- cs$index
  countloop <- 1
  while (is.finite(best.obj) && countloop < nCsteps) {
    cs <- CStep.fun(x, y, best.idx, h, hsize, alpha, lambda, scal)
    if (!is.finite(cs$object) || !better(cs$object, best.obj)) break
    best.obj  <- cs$object
    best.idx  <- cs$index
    countloop <- countloop + 1
  }
  list(index = best.idx, object = best.obj)
}


## ---- phase 1: elemental subsets --------------------------------------------

## glmnet warns that a class has fewer than 8 observations; this is expected for
## the tiny elemental subsets (2 + 2 observations) and is not shown.
quietElemental <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("fewer than 8", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
}

initialSubsets <- function(spec, x, y, h, hsize, alpha, lambda, nsamp, scal, ncores, seed) {
  if (!is.null(seed)) set.seed(seed)
  elemental <- spec$subsets(x, y, nsamp)       # random numbers are drawn here,
                                               # in the main process only
  twoCsteps <- function(c) {
    cs1 <- quietElemental(
      spec$cstep(x, y, elemental[, c], h, hsize, alpha,
                 spec$lambda.first(lambda, h), scal = FALSE))
    cs2 <- spec$cstep(x, y, cs1$index, h, hsize, alpha,
                      spec$lambda.second(lambda, h), scal)
    list(obj = cs2$object, indx = cs2$index)
  }
  enetApply(seq_len(nsamp), twoCsteps, ncores = ncores)
}

## The `s1` best subsets of phase 1 (those with a finite criterion).
bestSubsets <- function(subsets, s1, maximize, alpha, lambda) {
  obj <- vapply(subsets, function(s) s$obj, numeric(1))
  ord <- order(obj, decreasing = maximize)
  ord <- ord[is.finite(obj[ord])]
  ord <- ord[seq_len(min(s1, length(ord)))]
  if (length(ord) == 0) {
    stop("Model is not suitable for alpha ", alpha, " lambda ", lambda,
         " for this data set. Choose another lambda.", call. = FALSE)
  }
  lapply(ord, function(i) subsets[[i]]$indx)
}


## ---- phase 2: C-steps on the best subsets ----------------------------------
startCsteps <- function(spec, x, y, h, hsize, alpha, lambda, nsamp, s1, scal,
                        ncores, nCsteps, seed) {
  subsets <- initialSubsets(spec, x, y, h, hsize, alpha, lambda, nsamp, scal, ncores, seed)
  best    <- bestSubsets(subsets, s1, spec$maximize, alpha, lambda)

  converge <- function(index) {
    iterateCsteps(spec$cstep, x, y, index, h, hsize, alpha,
                  spec$lambda.iter(lambda, h), scal, nCsteps, spec$maximize)
  }
  res <- enetApply(best, converge, ncores = ncores)

  obj <- vapply(res, function(r) r$object, numeric(1))
  res[[if (spec$maximize) which.max(obj) else which.min(obj)]]$index
}


## ---- phase 3: the (alpha, lambda) grid -------------------------------------

## Follow the grid along `values` (a sequence of alphas or lambdas), starting
## from `index`; every C-step sequence starts from the previous result.
## `param(v)` returns the (alpha, lambda) pair for the value `v`.
walkGrid <- function(spec, x, y, h, hsize, values, param, index, scal, nCsteps) {
  out <- matrix(NA_integer_, nrow = h, ncol = length(values))
  for (i in seq_along(values)) {
    p   <- param(values[i])
    res <- iterateCsteps(spec$cstep, x, y, index, h, hsize, p$alpha, p$lambda,
                         scal, nCsteps, spec$maximize)
    index     <- res$index
    out[, i]  <- index[seq_len(h)]
  }
  out
}

## Warm-started C-steps for one fixed alpha over the whole lambda sequence
## (columns 2, 3, ...); every lambda starts from the same subset `start`.
walkLambdas <- function(cc, spec, x, y, h, hsize, alphas, lambdas, start, scal, nCsteps) {
  out <- matrix(NA_integer_, nrow = h, ncol = length(lambdas) - 1)
  for (ii in seq_along(lambdas)[-1]) {
    res <- iterateCsteps(spec$cstep, x, y, start[, cc], h, hsize, alphas[cc],
                         lambdas[ii], scal, nCsteps, spec$maximize)
    out[, ii - 1] <- res$index[seq_len(h)]
  }
  cbind(start[, cc], out)
}

## Entry point used by the family functions: returns the subsets for the whole
## grid, a list with one h x length(lambdas) matrix per alpha (or a single
## subset if there is only one alpha and one lambda).
warmCsteps.mod <- function(x, y, h, hsize, alphas, lambdas, nsamp, s1, scal,
                           ncores, nCsteps, tol, seed, family) {
  spec   <- familySpec(family)
  lambda <- spec$lambda.start(lambdas)
  first  <- startCsteps(spec, x, y, h, hsize, alphas[1], lambda, nsamp, s1, scal,
                        ncores, nCsteps, seed)

  nA <- length(alphas)
  nL <- length(lambdas)
  if (nA == 1 && nL == 1) return(first)

  # Subsets in the grid have exactly h observations. (For the multinomial
  # family the class-wise subsets can be slightly longer, so truncate.)
  first.h <- first[seq_len(h)]

  # alpha-direction at the first lambda: one start subset for every alpha
  start <- cbind(first.h,
                 walkGrid(spec, x, y, h, hsize, alphas[-1],
                          function(a) list(alpha = a, lambda = lambdas[1]),
                          first, scal, nCsteps))

  # Result: a list with one h x length(lambdas) matrix per alpha.
  if (nL == 1) {
    return(lapply(seq_len(nA), function(i) start[, i, drop = FALSE]))
  }
  if (nA == 1) {
    grid <- cbind(first.h,
                  walkGrid(spec, x, y, h, hsize, lambdas[-1],
                           function(l) list(alpha = alphas, lambda = l),
                           first, scal, nCsteps))
    return(list(grid))
  }
  enetApply(seq_len(nA), function(cc) {
    walkLambdas(cc, spec, x, y, h, hsize, alphas, lambdas, start, scal, nCsteps)
  }, ncores = ncores)
}


## ---- default penalty grid for the multinomial family -----------------------

## Data-driven grid: the smallest lambda that sets all coefficients to zero for
## alpha = 1 (lambda_max of glmnet, on the robustly standardized data used by
## the C-steps) down to 1% of it, equally spaced on the log scale. A fixed grid
## cannot match the scale of the penalty, which depends on n, p and the data.
multinomLambdas <- function(x, y, nlambda = 20, ratio = 0.01) {
  lmax <- glmnet::glmnet(x, y, family = "multinomial", alpha = 1,
                         intercept = FALSE, standardize = FALSE)$lambda[1]
  lmax * exp(seq(0, log(ratio), length.out = nlambda))
}


## ---- reweighting rule for the multinomial family ---------------------------

## Per-observation multinomial deviance (negative log-likelihood of the observed
## class) from the linear predictors `eta` (n x k) and the indicator matrix `y`.
multinomDeviance <- function(eta, y) {
  m <- apply(eta, 1, max)
  m + log(rowSums(exp(eta - m))) - rowSums(y * eta)
}

## 0/1 weights of the reweighting step. An observation is removed if it is an
## outlier in the score space (Mahalanobis distance `md`, leverage points) or if
## the fitted model gives its observed class a probability below 1/6 (label
## outliers; approximately the cut-off of the Pearson residual rule of the binomial
## family with del = 0.0125).
## Leverage points are fitted well by a penalized model, so the deviance alone
## cannot find them; the label outliers do not show up in the distances.
multinomWeights <- function(md, eta, y, cut = log(6)) {
  as.integer(weight.multinomial(md) != 0 & multinomDeviance(eta, y) <= cut)
}


## Robust distances of the observations of every class in the score space `z`:
## location and scatter from the MCD of the class, distances rescaled by their
## median (so that the median equals the chi-square median).
multinomMD <- function(z, y) {
  md <- rep(NA_real_, nrow(z))
  for (j in seq_len(ncol(y))) {
    in.j  <- y[, j] == 1
    mcd   <- covMcd(z[in.j, , drop = FALSE])
    m     <- sqrt(mahalanobis(z[in.j, , drop = FALSE], mcd$center, mcd$cov))
    md[in.j] <- m * sqrt(qchisq(0.5, ncol(z))) / median(m)
  }
  md
}

## Weights of the reweighting step from the linear predictors `eta` (n x k) of a
## fit: first the principal components of the (centered) predictors, then the
## robust distances within the classes. If the predictors have less than k - 1
## dimensions (e.g. all coefficients are zero) there is no spread to judge the
## observations by, and only the deviance rule applies.
multinomReweight <- function(eta, y) {
  sv   <- svd(scale(eta, TRUE, FALSE))
  rank <- sum(sv$d > 1e-4)
  md   <- if (rank < ncol(y) - 1) rep(0, nrow(eta))
          else multinomMD(sv$u[, seq_len(rank), drop = FALSE], y)
  multinomWeights(md, eta, y)
}
