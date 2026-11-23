print.enetLTS <-
   function(x, vers = c("reweighted", "raw"), digits = max(3, getOption("digits") - 3),
            maxprint = 20, ...){

      ## Compact summary: only the nonzero coefficients are shown, with the
      ## variable names (x1, x2, ... if the predictors have no column names).
      vers   <- match.arg(vers)
      family <- x$inputs$family
      cat("enetLTS estimator (family = \"", family, "\", ", vers, " fit)\n", sep = "")
      cat("\nCall: ", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")

      xx  <- x$inputs$xx
      int <- isTRUE(x$inputs$intercept)
      vn  <- colnames(xx)
      if (!is.null(vn)) vn <- vn[vn != "(Intercept)"]
      p   <- ncol(xx) - as.integer(int && !is.null(colnames(xx)) && "(Intercept)" %in% colnames(xx))
      if (is.null(vn) || length(vn) != p || all(vn == "")) vn <- paste0("x", seq_len(p))
      labs <- if (int) c("(Intercept)", vn) else vn

      b <- coef.enetLTS(x, vers = vers)
      if (is.list(b)) b <- b[[1]]

      if (is.matrix(b)) {
         rownames(b) <- labs[seq_len(nrow(b))]
         keep  <- which(rowSums(b != 0) > 0)
         slope <- if (int) setdiff(keep, 1) else keep
         nnz   <- length(slope)
         cat("\nNonzero coefficients (", nnz, " of ", p, " predictors in at least one class):\n", sep = "")
         show <- keep[seq_len(min(length(keep), maxprint))]
         print(round(b[show, , drop = FALSE], digits))
         if (length(keep) > maxprint)
            cat("... ", length(keep) - maxprint, " more rows; see coef(x)\n", sep = "")
      } else {
         b <- as.numeric(b); names(b) <- labs[seq_along(b)]
         keep  <- which(b != 0)
         slope <- if (int) setdiff(keep, 1) else keep
         nnz   <- length(slope)
         cat("\nNonzero coefficients (", nnz, " of ", p, " predictors):\n", sep = "")
         show <- keep[seq_len(min(length(keep), maxprint))]
         print(round(b[show], digits))
         if (length(keep) > maxprint)
            cat("... ", length(keep) - maxprint, " more; see coef(x)\n", sep = "")
      }

      w <- if (vers == "reweighted") x$wt else x$raw.wt
      cat("\nalpha = ", format(x$alpha, digits = digits),
          ", lambda = ", format(x$lambda, digits = digits),
          ", lambdaw = ", format(x$lambdaw, digits = digits), "\n", sep = "")
      cat("Observations flagged as outliers: ", sum(w == 0), " of ", length(w),
          " (see outliers.enetLTS())\n", sep = "")
      invisible(x)
}
