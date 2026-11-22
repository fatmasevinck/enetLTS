## -----------------------------------------------------------------------------
## Cross-validation over the (alpha, lambda) grid
##
## For every (alpha, lambda) pair the optimal subset found by the C-steps is
## split into `nfold` folds (`repl` times), the elastic net is fitted on the
## training folds and the loss is evaluated on the test fold. The pair with the
## smallest criterion wins. The three families only differ in
##
##   * how the folds are formed (random for gaussian, stratified by class for
##     binomial and multinomial),
##   * the penalty used on a training set,
##   * the loss and how it is summarized;
##
## see `cvSpec()`. The three public functions at the end of the file are thin
## wrappers around `cv.core()`.
## -----------------------------------------------------------------------------


## ---- helpers ----------------------------------------------------------------

## Fold membership of every observation: random within each group of
## `groups`, folds of (almost) equal size. This is `cvTools::cvFolds(type =
## "random")` applied group-wise and consumes random numbers in the same order,
## so a given seed gives the same folds as before.
cvFoldIds <- function(groups, nfold) {
  ids <- integer(length(groups))
  for (g in sort(unique(groups))) {
    idx <- which(groups == g)
    if (nfold < 2 || nfold > length(idx)) {
      stop("'nfold' outside allowable range for a group of ", length(idx),
           " observations", call. = FALSE)
    }
    ids[idx[sample.int(length(idx))]] <- rep_len(seq_len(nfold), length(idx))
  }
  ids
}

## Trimmed summary of a vector of losses: the mean over the smallest
## `1 - trim` fraction of the values.
uptrim <- function(x, trim = 0.1) {
  sort(x)[seq_len(length(x) * (1 - trim))]
}

## Numerically stable log(1 + exp(eta)).
log1pexp <- function(eta) pmax(eta, 0) + log1p(exp(-abs(eta)))


## ---- family-specific parts --------------------------------------------------

## groups:    grouping of the observations used to stratify the folds
## lambda:    penalty to use on a training set of `ntrain` observations
## loss:      per-observation loss of a model fitted on a training set
## criterion: summary of the losses of one repetition; `trimmed` is TRUE when
##            the subset was not found by the C-steps (`index = NULL`) and the
##            losses therefore still contain outliers
cvSpec <- function(family) {
  switch(family,
    gaussian = list(
      glm.family = "gaussian", label = c("RTMSPE", "RMSPE"),
      groups = function(y) rep(1L, length(y)),
      lambda = function(lambda, ntrain) lambda / ntrain,
      loss = function(fit, x, y) drop(y - x %*% as.matrix(fit$beta)),
      criterion = function(loss, groups, trimmed) {
        sqrt(if (trimmed) uptrimCrit(loss^2) else mean(loss^2, na.rm = TRUE))
      }),
    binomial = list(
      glm.family = "binomial", label = c("TMNLL", "MNLL"),
      groups = function(y) y,
      lambda = function(lambda, ntrain) lambda / ntrain,
      loss = function(fit, x, y) {
        eta <- drop(x %*% as.matrix(fit$beta))
        log1pexp(eta) - y * eta
      },
      criterion = function(loss, groups, trimmed) {
        if (trimmed) uptrimCrit(loss) else mean(loss, na.rm = TRUE)
      }),
    multinomial = list(
      glm.family = "multinomial", label = c("TMNLL", "MNLL"),
      groups = function(y) max.col(y, ties.method = "first"),
      lambda = function(lambda, ntrain) lambda,
      loss = function(fit, x, y) {
        eta <- matrix(predict(fit, newx = x, type = "link"), nrow = nrow(x))
        multinomDeviance(eta, y)
      },
      # the losses are trimmed within each class
      criterion = function(loss, groups, trimmed) {
        mean(unlist(lapply(split(loss, groups), uptrim)), na.rm = TRUE)
      }),
    stop("unknown family: ", family))
}

## Rows of `y` (vector or indicator matrix).
rowsOf <- function(y, i) if (is.matrix(y)) y[i, , drop = FALSE] else y[i]


## ---- one grid point ---------------------------------------------------------
cvEvalCombo <- function(task, spec, index, xx, yy, alphas, lambdas, nfold, trimmed) {
  alpha  <- alphas[task$i]
  lambda <- lambdas[task$j]

  if (is.null(index)) {
    x <- xx; y <- yy
  } else {
    rows <- index[[task$i]][, task$j]
    x <- xx[rows, , drop = FALSE]
    y <- rowsOf(yy, rows)
  }
  groups <- spec$groups(y)

  crit <- vapply(task$folds, function(fold.id) {
    loss <- rep(NA_real_, length(fold.id))
    for (f in seq_len(nfold)) {
      train <- fold.id != f
      test  <- !train
      fit <- tryCatch(
        glmnet(x[train, , drop = FALSE], rowsOf(y, train), family = spec$glm.family,
               alpha = alpha, lambda = spec$lambda(lambda, sum(train)),
               standardize = FALSE, intercept = FALSE),
        error = function(e) NULL)
      if (is.null(fit)) {
        message("CV broke off for alpha=", alpha, " and lambda=", lambda)
      } else {
        loss[test] <- spec$loss(fit, x[test, , drop = FALSE], rowsOf(y, test))
      }
    }
    spec$criterion(loss, groups, trimmed)
  }, numeric(1))

  mean(crit, na.rm = TRUE)
}


## ---- the cross-validation ---------------------------------------------------
cv.core <- function(family, index, xx, yy, alphas, lambdas, nfold, repl, ncores,
                    plot = TRUE) {
  if (any(alphas < 0 | alphas > 1)) stop("alphas can take the values only between 0 and 1")
  if (missing(lambdas) || length(lambdas) == 0) stop("provide a lambdas sequence")

  spec    <- cvSpec(family)
  trimmed <- is.null(index)
  nA <- length(alphas)
  nL <- length(lambdas)

  # Grid points, alpha varying fastest. The random fold assignments are drawn
  # here, in the main process and in grid order, so that the result does not
  # depend on the number of cores.
  combos <- expand.grid(i = seq_len(nA), j = seq_len(nL))
  tasks <- lapply(seq_len(nrow(combos)), function(r) {
    i <- combos$i[r]; j <- combos$j[r]
    y <- if (is.null(index)) yy else rowsOf(yy, index[[i]][, j])
    groups <- spec$groups(y)
    list(i = i, j = j, folds = replicate(repl, cvFoldIds(groups, nfold), simplify = FALSE))
  })

  evalCrit <- unlist(enetApply(tasks, cvEvalCombo, spec = spec, index = index,
                               xx = xx, yy = yy, alphas = alphas, lambdas = lambdas,
                               nfold = nfold, trimmed = trimmed, ncores = ncores))
  evalCrit <- matrix(evalCrit, nrow = nA)          # rows: alphas, columns: lambdas
  dimnames(evalCrit) <- list(paste("alpha", alphas), paste("lambda", lambdas))

  if (all(is.na(evalCrit))) stop("cross-validation failed for every (alpha, lambda) pair")
  optind <- which(evalCrit == min(evalCrit, na.rm = TRUE), arr.ind = TRUE)[1, ]
  a.opt  <- optind[[1]]
  l.opt  <- optind[[2]]

  out <- list(indexbest   = if (is.null(index)) NULL else index[[a.opt]][, l.opt],
              evalCrit    = evalCrit,
              minevalCrit = evalCrit[a.opt, l.opt],
              lambdaopt   = round(lambdas, 4)[l.opt],
              alphaopt    = round(alphas, 4)[a.opt])

  if (isTRUE(plot)) cvPlot(evalCrit, alphas, lambdas, family, trimmed, out$alphaopt, out$lambdaopt, spec)
  out
}


## Heat map of the criterion over the grid.
cvPlot <- function(evalCrit, alphas, lambdas, family, trimmed, alpha, lambda, spec) {
  RTMSPE <- RMSPE <- TMNLL <- MNLL <- NULL            # silence R CMD check
  message("optimal model: lambda = ", lambda, " alpha = ", alpha)

  colors <- colorRampPalette(c("black", "blue2", "purple", "orange", "yellow"))(
    length(alphas) * length(lambdas))
  label  <- spec$label[if (trimmed) 1 else 2]

  ## long format: one row per (alpha, lambda) pair; evalCrit has rows = alphas
  d <- data.frame(alpha  = rep(round(alphas, 4), times = length(lambdas)),
                  lambda = rep(round(lambdas, 4), each = length(alphas)),
                  crit   = as.vector(evalCrit))
  names(d)[3] <- label

  p <- ggplot(d, aes(x = as.factor(lambda), y = as.factor(alpha), fill = .data[[label]])) +
    geom_tile() +
    scale_fill_gradientn(colours = colors) +
    theme(axis.text.x = element_text(angle = -90)) +
    ggtitle(paste0(label, " (", if (trimmed) "minimum" else "optimal", " at lambda=",
                   lambda, ", alpha=", alpha, ", ", family, ")")) +
    xlab("lambda") + ylab("alpha")

  print(p)
}


## ---- wrappers used by the family functions and by cv.enetLTS() -------------
cv.gaussian.enetLTS <- function(index = NULL, xx, yy, alphas, lambdas, nfold, repl,
                                ncores, plot = TRUE) {
  cv.core("gaussian", index, xx, yy, alphas, lambdas, nfold, repl, ncores, plot)
}

cv.binomial.enetLTS <- function(index = NULL, xx, yy, alphas, lambdas, nfold, repl,
                                ncores, plot = TRUE) {
  cv.core("binomial", index, xx, yy, alphas, lambdas, nfold, repl, ncores, plot)
}

cv.multinom.enetLTS <- function(index = NULL, xx, yy, alphas, lambdas, nfold, repl,
                                ncores, plot = TRUE) {
  cv.core("multinomial", index, xx, yy, alphas, lambdas, nfold, repl, ncores, plot)
}
