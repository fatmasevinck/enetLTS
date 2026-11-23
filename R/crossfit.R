## -----------------------------------------------------------------------------
## Cross-fitted trimming (gaussian family, optional: enetLTS(..., crossfit = TRUE))
##
## The ordinary C-step keeps the h observations that the *current* fit describes
## best. With many predictors and a weak penalty the fit adapts to the noise, so
## the retained observations are partly chosen *because* noise variables explain
## them; every later step (reweighting, choice of alpha and lambda) then inherits
## a subset that favours false positives. Here the residuals that decide who is
## kept are cross-fitted: an observation is never judged by a model that has seen
## it.
## -----------------------------------------------------------------------------

## Alphas used inside the refinement: positive values of the user's grid, thinned.
crossfitAlphas <- function(alphas, k = 6) {
  a <- sort(unique(alphas[alphas > 0]))
  if (length(a) == 0) return(max(alphas))
  if (length(a) > k) a <- a[unique(round(seq(1, length(a), length.out = k)))]
  a
}

## alpha and lambda by cross-validation on the observations `H`:
## best alpha by the smallest CV error, lambda by the one-standard-error rule.
crossfitPick <- function(x, y, H, alphas, nfold) {
  fid <- sample(rep(seq_len(nfold), length.out = length(H)))
  cvs <- lapply(alphas, function(a)
    cv.glmnet(x[H, , drop = FALSE], y[H], alpha = a, foldid = fid, standardize = FALSE))
  i <- which.min(vapply(cvs, function(z) min(z$cvm), numeric(1)))
  list(alpha = alphas[i], lambda = cvs[[i]]$lambda.1se)
}

## Cross-fitted residuals of all n observations from a model trained on `H`:
## observation i is predicted by a model fitted on H without i's fold.
crossfitResiduals <- function(x, y, H, alpha, lambda, nfold) {
  n <- length(y); fold <- sample(rep(seq_len(nfold), length.out = n)); r <- numeric(n)
  for (k in seq_len(nfold)) {
    te <- fold == k; tr <- setdiff(H, which(te))
    g <- glmnet(x[tr, , drop = FALSE], y[tr], alpha = alpha, lambda = lambda, standardize = FALSE)
    r[te] <- drop(y[te] - predict(g, x[te, , drop = FALSE]))
  }
  r
}

## Iterates: choose (alpha, lambda) on H, cross-fit residuals, keep the h smallest.
## `x` must be robustly standardized. Returns the final h-subset, the (alpha,
## lambda) used for it and the reweighting set W (cross-fitted |residual| within
## qnorm(1 - del) robust standard deviations).
crossfitTrim <- function(x, y, H0, h, alphas, nfold, del, iters = 5) {
  ag <- crossfitAlphas(alphas); H <- H0
  for (it in seq_len(iters)) {
    p <- crossfitPick(x, y, H, ag, nfold)
    r <- crossfitResiduals(x, y, H, p$alpha, p$lambda, nfold)
    Hn <- order(abs(r))[seq_len(h)]
    if (setequal(Hn, H)) break
    H <- Hn
  }
  s <- 1.4826 * median(abs(r))
  list(H = H, alpha = p$alpha, lambda = p$lambda, W = which(abs(r) <= qnorm(1 - del) * s))
}

## Final choice of alpha and lambda on the (already standardized) reweighted data:
## all alphas share the same folds, the CV error is averaged over `nrep` fold
## assignments, alpha minimizes it and lambda follows the one-standard-error rule.
crossfitFinal <- function(xw, yw, nfold, nrep = 5, alphas = seq(0, 1, length.out = 21)) {
  n <- nrow(xw)
  res <- lapply(alphas, function(a) {
    L <- glmnet(xw, yw, alpha = a, standardize = FALSE, intercept = FALSE)$lambda
    M <- vapply(seq_len(nrep), function(r) {
      fd <- sample(rep(seq_len(nfold), length.out = n))
      cv.glmnet(xw, yw, alpha = a, lambda = L, foldid = fd, standardize = FALSE,
                intercept = FALSE)$cvm }, numeric(length(L)))
    M <- matrix(M, nrow = length(L))
    list(L = L, cvm = rowMeans(M), se = apply(M, 1, sd) / sqrt(nrep))
  })
  ia <- which.min(vapply(res, function(z) min(z$cvm), numeric(1)))
  z <- res[[ia]]; il <- which.min(z$cvm)
  bound <- z$cvm[il] + max(z$se[il], sqrt(.Machine$double.eps))
  list(alpha = alphas[ia], lambda = z$L[min(which(z$cvm <= bound))])
}


## ---- multinomial ------------------------------------------------------------
## Same idea for the multinomial family; the loss is -log p(y_i) and the h-subset
## is class-balanced (hk_l = floor(n_l * hsize + 1) observations per class).

multinomCrossfitAlphas <- function(alphas) {
  a <- c(0.05, 0.25, 0.5, 0.75, 1)
  a <- a[a >= min(alphas) & a <= max(alphas)]
  if (length(a) == 0) max(alphas) else a
}

## Cross-fitted -log p(y_i) of every observation in `pool`; folds stratified by class.
multinomCrossfitLoss <- function(x, cls, K, H, pool, alpha, lambda, nfold) {
  fold <- integer(length(cls))
  for (l in seq_len(K)) { ix <- pool[cls[pool] == l]; fold[ix] <- sample(rep(seq_len(nfold), length.out = length(ix))) }
  loss <- rep(Inf, length(cls)); P <- matrix(NA_real_, length(cls), K)
  for (k in seq_len(nfold)) {
    te <- pool[fold[pool] == k]; tr <- setdiff(H, te)
    g <- glmnet(x[tr, , drop = FALSE], factor(cls[tr], levels = seq_len(K)), family = "multinomial",
                alpha = alpha, lambda = lambda, standardize = FALSE)
    eta <- matrix(predict(g, x[te, , drop = FALSE], type = "link"), nrow = length(te))
    eta <- eta - apply(eta, 1, max)
    lp  <- eta - log(rowSums(exp(eta)))
    loss[te] <- -lp[cbind(seq_along(te), cls[te])]
    P[te, ] <- exp(lp)
  }
  attr(loss, "P") <- P
  loss
}

## Model-calibrated cut: the largest c such that, if the class probabilities in P were
## the true ones, the expected share of observations with p(y_i) < c is at most `target`.
multinomProbCut <- function(P, target) {
  P <- P[stats::complete.cases(P), , drop = FALSE]
  f <- function(c) mean(rowSums(P * (P < c)))
  lo <- 0
  for (c in sort(unique(as.vector(P)))) { if (f(c) <= target) lo <- c else break }
  lo
}

multinomCrossfitPick <- function(x, cls, K, H, alphas, nfold) {
  fid <- sample(rep(seq_len(nfold), length.out = length(H)))
  cvs <- lapply(alphas, function(a)
    quietPathEnd(cv.glmnet(x[H, , drop = FALSE], factor(cls[H], levels = seq_len(K)), family = "multinomial",
                           alpha = a, foldid = fid, standardize = FALSE, type.measure = "deviance")))
  i <- which.min(vapply(cvs, function(z) min(z$cvm), numeric(1)))
  list(alpha = alphas[i], lambda = cvs[[i]]$lambda.1se)
}

## `x` robustly standardized, `cls` the class (1..K) of every row. `H0` is the
## starting subset (rows of x). Returns the final subset H, the (alpha, lambda)
## used for it and the reweighting set W: within every class, the observations
## whose cross-fitted sqrt(loss) lies within `wk` robust standard deviations
## (median / MAD) of the class [option enetLTS.wrule = "mad"; see below for the
## default rule]. wk = 3.5 (not qnorm(1 - del) = 2.24): in the simulations
## the narrower cut dropped well-labeled boundary observations and made the fit too
## sparse. `del` is used by the default rule (cut calibrated to 2 * del).
multinomCrossfitTrim <- function(x, cls, K, H0, hsize, alphas, nfold, del, iters = 4, wk = 3.5) {
  pool <- seq_len(nrow(x)); ag <- multinomCrossfitAlphas(alphas); H <- H0
  hk <- floor(tabulate(cls, K) * hsize + 1)
  for (it in seq_len(iters)) {
    p    <- multinomCrossfitPick(x, cls, K, H, ag, nfold)
    loss <- multinomCrossfitLoss(x, cls, K, H, pool, p$alpha, p$lambda, nfold)
    Hn <- unlist(lapply(seq_len(K), function(l) {
      ix <- pool[cls == l]; ix[order(loss[ix])][seq_len(min(hk[l], length(ix)))] }))
    if (setequal(Hn, H)) break
    H <- Hn
  }
  if (identical(getOption("enetLTS.wrule", "full"), "mad")) {
    ## option enetLTS.wrule = "mad": within every class, cross-fitted sqrt(loss)
    ## within wk robust standard deviations of the class median
    W <- unlist(lapply(seq_len(K), function(l) {
      ix <- pool[cls == l]; s <- sqrt(loss[ix])
      ix[s <= median(s) + wk * 1.4826 * mad(s)] }))
  } else {
    ## default: probabilities from a cross-fitted model on all (screened) observations.
    ## The model fitted on the trimmed subset H separates the classes almost perfectly
    ## and gives over-confident probabilities, so clean boundary observations looked
    ## outlying. The cut is calibrated: if the probabilities were correct, a share of
    ## 2 * del of the observations would be flagged.
    pf <- multinomCrossfitPick(x, cls, K, pool, ag, nfold)
    lf <- multinomCrossfitLoss(x, cls, K, pool, pool, pf$alpha, pf$lambda, nfold)
    cc <- multinomProbCut(attr(lf, "P"), 2 * del)
    W  <- pool[exp(-lf[pool]) >= cc]
  }
  list(H = H, alpha = p$alpha, lambda = p$lambda, W = sort(W))
}

## alpha and lambda on the reweighted (standardized) data: common folds, alpha
## with the smallest CV deviance, lambda at its minimum.
## The alpha grid of the final (reweighted) fit is restricted to the range of
## the user-supplied `alphas`; if no grid value falls in that range, the
## user-supplied values themselves are used.
finalAlphaGrid <- function(alphas, grid) {
  a <- grid[grid >= min(alphas) - 1e-12 & grid <= max(alphas) + 1e-12]
  if (length(a) == 0) sort(unique(alphas)) else a
}

multinomCrossfitFinal <- function(xw, yw, nfold, typ, alphas = c(0, 0.1, 0.3, 0.5, 0.7, 1)) {
  fid <- sample(rep(seq_len(nfold), length.out = nrow(xw)))
  cvs <- lapply(alphas, function(a)
    try(quietPathEnd(cv.glmnet(xw, yw, family = "multinomial", alpha = a, foldid = fid, standardize = FALSE,
                               intercept = TRUE, type.multinomial = typ, type.measure = "deviance")),
        silent = TRUE))
  ok  <- !vapply(cvs, inherits, logical(1), "try-error")
  cvs <- cvs[ok]; alphas <- alphas[ok]
  i   <- which.min(vapply(cvs, function(z) min(z$cvm), numeric(1)))
  list(alpha = alphas[i], lambda = cvs[[i]]$lambda.min)
}
