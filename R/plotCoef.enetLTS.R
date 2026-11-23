plotCoef.enetLTS <- function(object, vers = c("reweighted", "raw"), colors = NULL, ...){

  ## Needle plot of the slope coefficients against the predictor index: one vertical
  ## segment from zero per coefficient, a point marks the nonzero ones. The intercept
  ## is given in the subtitle. For the multinomial family there is one panel per class.
  ## `colors` is kept for compatibility and is not used.

  vers   <- match.arg(vers)
  family <- object$inputs$family
  B      <- if (vers == "reweighted") object$coefficients else object$raw.coefficients
  int    <- isTRUE(object$inputs$intercept)
  fit    <- if (vers == "reweighted") "reweighted fit" else "raw fit"
  fam    <- switch(family, gaussian = "linear regression", binomial = "logistic regression",
                   multinomial = "multinomial regression")

  if (family == "multinomial") {
    B <- as.matrix(B)
    a0 <- if (int) B[1, ] else rep(0, ncol(B))
    S  <- if (int) B[-1, , drop = FALSE] else B
    K  <- ncol(S)
    cl <- object$classnames
    if (is.null(cl) || length(cl) != K) cl <- as.character(seq_len(K))
    d  <- data.frame(index = rep(seq_len(nrow(S)), K),
                     class = factor(rep(paste("class", cl), each = nrow(S)), levels = paste("class", cl)),
                     coef  = as.vector(S))
    sub <- if (int) paste0("Intercepts: ", paste(paste0(cl, " = ", signif(a0, 3)), collapse = ", ")) else NULL
    nz  <- d[d$coef != 0, ]
    plot <- ggplot(d, aes(x = index, y = coef, colour = class)) +
      geom_hline(yintercept = 0, colour = "gray60") +
      geom_segment(aes(xend = index, yend = 0), linewidth = 0.5) +
      geom_point(data = nz, size = 1.6) +
      facet_wrap(~ class, ncol = 1) +
      .classColour(K, name = "Class") +
      scale_x_continuous("Predictor index", breaks = coefBreaks(nrow(S))) +
      labs(title = paste0("enetLTS coefficients, ", fam, " (", fit, ")"), subtitle = sub, y = "Coefficient") +
      .enetLTS_theme(legend = FALSE)
  } else {
    B  <- as.numeric(B)
    a0 <- if (int) B[1] else 0
    S  <- if (int) B[-1] else B
    d  <- data.frame(index = seq_along(S), coef = S)
    sub <- paste0(if (int) paste0("Intercept = ", signif(a0, 3), ";  ") else "",
                  sum(S != 0), " of ", length(S), " slopes nonzero")
    plot <- ggplot(d, aes(x = index, y = coef)) +
      geom_hline(yintercept = 0, colour = "gray60") +
      geom_segment(aes(xend = index, yend = 0), colour = .enetLTS_pal$line, linewidth = 0.5) +
      geom_point(data = d[d$coef != 0, ], colour = .enetLTS_pal$line, size = 1.6) +
      scale_x_continuous("Predictor index", breaks = coefBreaks(length(S))) +
      labs(title = paste0("enetLTS coefficients, ", fam, " (", fit, ")"), subtitle = sub, y = "Coefficient") +
      .enetLTS_theme(legend = FALSE)
  }
  print(plot)
  invisible(plot)
}

## x-axis breaks: every index for a few coefficients, otherwise about 10 rounded
## positions so that the labels stay readable
coefBreaks <- function(K) if (K <= 30) seq_len(K) else unique(pmax(1, pretty(c(1, K), n = 10)))
