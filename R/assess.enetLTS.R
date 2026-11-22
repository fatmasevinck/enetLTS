## Assessment of a fitted enetLTS model on (new) data -------------------------

.enetLTS_pred <- function(object, newX, vers) {
  family <- object$inputs$family
  cn     <- object$classnames
  ## predict() names the element "response" (reweighted) or "raw.response" (raw)
  P <- predict.enetLTS(object, newX, vers = vers, type = "response")[[1]]
  if (family == "gaussian") {
    list(family = family, pred = drop(P))
  } else if (family == "binomial") {
    list(family = family, cn = cn, pred = drop(P))
  } else {
    list(family = family, cn = cn, pred = unname(as.matrix(P)))
  }
}

.enetLTS_y01 <- function(newy, cn) {
  idx <- match(as.character(newy), cn)
  if (anyNA(idx)) stop("'newy' contains labels that are not in the classes of the fit: ",
                       paste(cn, collapse = ", "))
  idx
}

.auc <- function(score, y01) {
  n1 <- sum(y01 == 1); n0 <- sum(y01 == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[y01 == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

assess.enetLTS <- function(object, newX, newy, vers = c("both", "reweighted", "raw")) {
  ## Performance measures on test data, in the layout of glmnet::assess.glmnet:
  ## a list with one element per measure; each element has one value per fit
  ## (reweighted and/or raw) and the attribute "measure".
  vers  <- match.arg(vers)
  newX  <- as.matrix(newX)
  if (nrow(newX) != length(newy)) stop("'newX' and 'newy' must have the same number of observations")
  versions <- if (vers == "both") c("reweighted", "raw") else vers
  one <- function(v) {
    pr <- .enetLTS_pred(object, newX, v)
    if (pr$family == "gaussian") {
      r <- as.numeric(newy) - pr$pred
      c(mse = mean(r^2), rmse = sqrt(mean(r^2)), mae = mean(abs(r)))
    } else if (pr$family == "binomial") {
      y <- .enetLTS_y01(newy, pr$cn) - 1          # 0/1, second class = 1
      p <- pmin(pmax(pr$pred, 1e-15), 1 - 1e-15)
      c(deviance = -2 * mean(y * log(p) + (1 - y) * log(1 - p)),
        class = mean((p > 0.5) != (y == 1)), auc = .auc(p, y),
        mse = mean((y - p)^2), mae = mean(abs(y - p)))
    } else {
      idx <- .enetLTS_y01(newy, pr$cn)
      P   <- pmin(pmax(pr$pred, 1e-15), 1)
      Y   <- diag(ncol(P))[idx, , drop = FALSE]
      c(deviance = -2 * mean(log(P[cbind(seq_along(idx), idx)])),
        class = mean(max.col(pr$pred, ties.method = "first") != idx),
        mse = mean(rowSums((Y - pr$pred)^2)), mae = mean(rowSums(abs(Y - pr$pred))))
    }
  }
  M <- sapply(versions, one)                      # measures x versions
  if (!is.matrix(M)) M <- matrix(M, ncol = length(versions), dimnames = list(names(M), versions))
  labels <- c(mse = "Mean-Squared Error", rmse = "Root Mean-Squared Error",
              mae = "Mean Absolute Error", class = "Misclassification Error", auc = "AUC",
              deviance = if (object$inputs$family == "binomial") "Binomial Deviance" else "Multinomial Deviance")
  out <- lapply(rownames(M), function(m) structure(setNames(M[m, ], colnames(M)), measure = unname(labels[m])))
  names(out) <- rownames(M)
  out
}

confusion.enetLTS <- function(object, newX, newy, vers = c("reweighted", "raw")) {
  vers <- match.arg(vers)
  if (object$inputs$family == "gaussian") stop("'confusion' is available for binomial and multinomial fits only")
  pr  <- .enetLTS_pred(object, as.matrix(newX), vers)
  cn  <- pr$cn
  idx <- .enetLTS_y01(newy, cn)
  hat <- if (pr$family == "binomial") as.integer(pr$pred > 0.5) + 1L else max.col(pr$pred, ties.method = "first")
  tab <- table(Predicted = factor(cn[hat], levels = cn), True = factor(cn[idx], levels = cn))
  structure(list(table = tab, correct = sum(diag(tab)) / sum(tab)), class = "confusion.enetLTS")
}

print.confusion.enetLTS <- function(x, ...) {
  print(addmargins(x$table))
  cat("\n Percent correct:", round(x$correct, 4), "\n")
  invisible(x)
}

roc.enetLTS <- function(object, newX, newy, vers = c("reweighted", "raw"), plot = FALSE, ...) {
  vers <- match.arg(vers)
  if (object$inputs$family != "binomial") stop("'roc' is available for binomial fits only")
  pr <- .enetLTS_pred(object, as.matrix(newX), vers)
  y  <- .enetLTS_y01(newy, pr$cn) - 1
  o  <- order(pr$pred, decreasing = TRUE)
  tp <- cumsum(y[o] == 1); fp <- cumsum(y[o] == 0)
  keep <- c(!duplicated(pr$pred[o], fromLast = TRUE))
  res <- data.frame(FPR = c(0, fp[keep] / sum(y == 0)), TPR = c(0, tp[keep] / sum(y == 1)))
  attr(res, "auc") <- .auc(pr$pred, y)
  if (isTRUE(plot)) {
    fit <- if (vers == "reweighted") "reweighted fit" else "raw fit"
    g <- ggplot2::ggplot(res, ggplot2::aes(x = FPR, y = TPR)) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "gray60") +
      ggplot2::geom_step(colour = .enetLTS_pal$line, linewidth = 0.8) +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(title = paste0("ROC curve (", fit, ")"),
                    subtitle = sprintf("AUC = %.3f", attr(res, "auc")),
                    x = "False positive rate", y = "True positive rate") +
      .enetLTS_theme(legend = FALSE)
    print(g)
  }
  invisible(res)
}

## Outlier listing --------------------------------------------------------------

outliers.enetLTS <- function(object, vers = c("reweighted", "raw")) {
  vers   <- match.arg(vers)
  family <- object$inputs$family
  X <- object$inputs$xx
  if (isTRUE(object$inputs$intercept)) X <- X[, -1, drop = FALSE]
  n <- nrow(X)
  w <- if (vers == "reweighted") object$wt else object$raw.wt
  flagged <- as.integer(w) == 0L
  screened <- seq_len(n) %in% object$screened
  pr <- .enetLTS_pred(object, X, vers)
  if (family == "gaussian") {
    y   <- as.numeric(object$inputs$yy)
    hat <- pr$pred
    loss <- abs(y - hat)
    type <- ifelse(flagged, "outlier", "regular")
    out <- data.frame(index = seq_len(n), flagged = flagged, loss = loss,
                      observed = y, predicted = hat, type = type,
                      stringsAsFactors = FALSE)
  } else {
    cn  <- pr$cn
    yf  <- if (family == "multinomial") as.character(object$inputs$y_factor) else cn[as.integer(object$inputs$yy) + 1L]
    idx <- .enetLTS_y01(yf, cn)
    P   <- if (family == "binomial") cbind(1 - pr$pred, pr$pred) else pr$pred
    P   <- pmin(pmax(P, 1e-15), 1)
    ## loss = -log p(observed class), computed from the linear predictor so that
    ## it is not capped for very confident wrong predictions
    eta <- predict.enetLTS(object, X, vers = vers, type = "link")[[1]]
    if (family == "binomial") {
      eta  <- drop(eta)
      sgn  <- ifelse(idx == 2, -1, 1)                      # observed class 2 -> -eta
      z    <- sgn * eta
      loss <- ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))  # log(1 + exp(z))
    } else {
      eta  <- unname(as.matrix(eta))
      mx   <- apply(eta, 1, max)
      lse  <- mx + log(rowSums(exp(eta - mx)))
      loss <- lse - eta[cbind(seq_len(n), idx)]
    }
    hat  <- cn[max.col(P, ties.method = "first")]
    mis  <- hat != yf
    type <- rep("regular", n)
    type[flagged & !mis] <- "leverage point"
    type[flagged &  mis] <- "possible mislabel / outlier"
    type[flagged & screened & !mis] <- "leverage point (screened)"
    type[flagged & screened &  mis] <- "bad leverage point (screened)"
    out <- data.frame(index = seq_len(n), flagged = flagged, loss = loss,
                      observed = yf, predicted = hat, type = type,
                      stringsAsFactors = FALSE)
  }
  class(out) <- c("outliers.enetLTS", "data.frame")
  attr(out, "vers") <- vers
  attr(out, "family") <- family
  out
}

print.outliers.enetLTS <- function(x, ...) {
  d <- as.data.frame(x, stringsAsFactors = FALSE)
  class(d) <- "data.frame"
  fl <- d[d$flagged, , drop = FALSE]
  cat("enetLTS outliers (", attr(x, "vers"), " fit, ", attr(x, "family"), "): ",
      nrow(fl), " of ", nrow(d), " observations flagged\n", sep = "")
  if (nrow(fl) > 0) {
    fl$loss <- round(fl$loss, 3)
    if (is.numeric(fl$observed)) fl$observed <- round(fl$observed, 3)
    if (is.numeric(fl$predicted)) fl$predicted <- round(fl$predicted, 3)
    fl$flagged <- NULL
    print(fl, row.names = FALSE)
  }
  invisible(x)
}
