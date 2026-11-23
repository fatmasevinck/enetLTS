coef.enetLTS <-
  function(object, vers = c("reweighted", "raw"), zeros = TRUE, ...)
  {
    vers <- match.arg(vers)
    b    <- if (vers == "reweighted") object$coefficients else object$raw.coefficients
    labs <- colnames(object$inputs$xx)       # "(Intercept)", then the predictor names
    if (is.matrix(b)) {                      # multinomial: one column per class
      if (!is.null(labs) && length(labs) == nrow(b)) rownames(b) <- labs
      if (!is.null(object$classnames) && length(object$classnames) == ncol(b))
        colnames(b) <- object$classnames
      if (!isTRUE(zeros)) b <- b[rowSums(b != 0) > 0, , drop = FALSE]
    } else {
      b <- as.numeric(b)
      names(b) <- if (!is.null(labs) && length(labs) == length(b)) labs else seq_along(b)
      if (!isTRUE(zeros)) b <- b[b != 0]
    }
    b
  }
