
plotResid.enetLTS <- function(object,vers=c("reweighted","raw"),zoom=TRUE,...){

  vers   <- match.arg(vers)
  family <- object$inputs$family
  del    <- if (is.null(object$inputs$del)) 0.0125 else object$inputs$del
  cutoff <- qnorm(1 - del)
  plots  <- list()          # the drawn plots, returned invisibly


  if (family=="multinomial"){
    y.factor <- object$inputs$y_factor
  }

  x <- object$inputs$xx
  y <- object$inputs$yy

  coefficients     <- object$coefficients
  raw.coefficients <- object$raw.coefficients

  if (family=="multinomial"){
    if (vers=="reweighted"){
      mainfit          <- paste("Negative log-likelihood vs index, multinomial regression")
      xlab             <- "Index"
      ylab             <- "Negative log-likelihood"
      lab.val          <- c("regular observation","outliers")
      classification   <- factor(as.numeric(!object$wt==0))
      ind              <- factor(as.numeric(object$wt==0))
      yhat             <- x %*% coefficients
      probs            <- multinomProbs(yhat)
      residuals        <- (-apply(y*log(probs),1,sum))
      Index            <- as.numeric(c(1:length(residuals)))
      names(residuals) <- 1:length(residuals)
      residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                     llim=residuals,ulim=residuals)
    } else if (vers=="raw"){
      mainfit          <- paste("Negative log-likelihood vs index, multinomial regression (raw fit)")
      xlab             <- "Index"
      ylab             <- "Negative log-likelihood (raw fit)"
      lab.val          <- c("regular observation","outliers")
      yhat             <- x %*% raw.coefficients
      probs            <- multinomProbs(yhat)
      residuals        <- (-apply(y*log(probs),1,sum))
      classification   <- factor(as.numeric(c(1:length(residuals)) %in% object$best))
      ind              <- factor(as.numeric(!c(1:length(residuals)) %in% object$best))
      Index            <- as.numeric(c(1:length(residuals)))
      names(residuals) <- 1:length(residuals)
      residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                     llim=residuals,ulim=residuals)
    }
    d    <- data.frame(X = Index, Y = residuals$residuals, cls = factor(y.factor), flag = ind == "1")
    plot <- .classPlot(d, nlevels(factor(y.factor)),
                       if (vers == "raw") "not in the best subset" else "flagged as outliers",
                       under = list(geom_hline(yintercept=0, linetype="dashed", color="gray")),
                       zx = FALSE, zoom = zoom) +
      labs(title = mainfit, x = xlab, y = ylab)
    print(plot); plots <- c(plots, list(plot))
  } else if (family=="binomial"){
    if (vers=="reweighted"){
      mainfit          <- paste("Negative log-likelihood vs index, logistic regression")
      xlab             <- "Index"
      ylab             <- "Negative log-likelihood"
      lab.val          <- c("regular observation","outliers")
      classification   <- factor(as.numeric(!object$wt==0))
      ind              <- factor(as.numeric(object$wt==0))
      residuals        <- -(y * x %*% coefficients) + log1pexp(x %*% coefficients)
      Index            <- as.numeric(c(1:length(residuals)))
      names(residuals) <- 1:length(residuals)
      residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                     llim=residuals,ulim=residuals)
    } else if (vers=="raw"){
      mainfit          <- paste("Negative log-likelihood vs index, logistic regression (raw fit)")
      xlab             <- "Index"
      ylab             <- "Negative log-likelihood (raw fit)"
      lab.val          <- c("best subset","outliers")
      residuals        <- -(y * x %*% raw.coefficients) + log1pexp(x %*% raw.coefficients)
      classification   <- factor(as.numeric(c(1:length(residuals)) %in% object$best))
      ind              <- factor(as.numeric(!c(1:length(residuals)) %in% object$best))
      Index            <- as.numeric(c(1:length(residuals)))
      names(residuals) <- 1:length(residuals)
      residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                     llim=residuals,ulim=residuals)
    }
    d    <- data.frame(X = Index, Y = residuals$residuals, cls = factor(y, levels = c(0, 1)), flag = ind == "1")
    plot <- .classPlot(d, 2, if (vers == "raw") "not in the best subset" else "flagged as outliers", under = list(geom_hline(yintercept=0, linetype="dashed", color="gray")), zx = FALSE, zoom = zoom) +
      labs(title = mainfit, x = xlab, y = ylab)
    print(plot); plots <- c(plots, list(plot))
  } else if (family=="gaussian"){
    if (vers=="reweighted"){
      mainfit          <- paste("Standardized residuals vs index, linear regression")
      xlab             <- "Index"
      ylab             <- "Standardized residuals"
      lab.val          <- c("regular observation","outliers")
      classification   <- factor(as.numeric(!object$wt==0))
      ind              <- factor(as.numeric(object$wt==0))
      residuals        <- y - x %*% coefficients
      ## standardized as in the reweighting step, so that the cut-off lines
      ## +/- qnorm(1 - del) separate exactly the flagged observations
      wg               <- weight.gaussian(residuals, object$raw.wt == 1, del)
      residuals        <- (residuals - wg$mu) / wg$s
      Index            <- as.numeric(c(1:length(residuals)))
      names(residuals) <- 1:length(residuals)
      residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                     llim=residuals,ulim=residuals)
    } else if (vers=="raw"){
      mainfit          <- paste("Standardized residuals vs index, linear regression (raw fit)")
      xlab             <- "Index"
      ylab             <- "Standardized raw residuals"
      lab.val          <- c("best subset","outliers")
      residuals        <- y - x %*% raw.coefficients
      wg               <- weight.gaussian(residuals, object$best, del)
      residuals        <- (residuals - wg$mu) / wg$s
      classification   <- factor(as.numeric(c(1:length(residuals)) %in% object$best))
      ind              <- factor(as.numeric(!c(1:length(residuals)) %in% object$best))
      Index            <- as.numeric(c(1:length(residuals)))
      names(residuals) <- 1:length(residuals)
      residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                     llim=residuals,ulim=residuals)
    }
    d    <- data.frame(X = Index, Y = residuals$residuals, flag = ind == "1")
    plot <- .diagPlot(d, lab.val, under = list(geom_hline(yintercept=0, linetype="dashed", color="gray"),
                                               geom_hline(yintercept=c(-cutoff, cutoff), color="gray60")),
                      zx = FALSE, zoom = zoom) +
      labs(title = mainfit, x = xlab, y = ylab)
    print(plot); plots <- c(plots, list(plot))
  }


  if (family=="binomial" | family=="gaussian"){
    if (family=="binomial"){
      if (vers=="reweighted"){
        mainfit          <- expression(bold(paste("Negative log-likelihood vs ", X*hat(beta), ", logistic regression")))
        xlab             <- expression(X*hat(beta))
        ylab             <- "Negative log-likelihood"
        lab.val          <- c("regular observation","outliers")
        classification   <- factor(as.numeric(!object$wt==0))
        ind              <- factor(as.numeric(object$wt==0))
        residuals        <- -(y * x %*% coefficients) + log1pexp(x %*% coefficients)
        fitted.values    <- x %*% coefficients
        names(residuals) <- 1:length(residuals)
        residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                       llim=residuals,ulim=residuals)
      } else if (vers=="raw"){
        mainfit          <- expression(bold(paste("Negative log-likelihood vs ", X*hat(beta)[raw], ", logistic regression")))
        xlab             <- expression(X*hat(beta)[raw])
        ylab             <- "Negative log-likelihood (raw fit)"
        lab.val          <- c("best subset","outliers")
        residuals        <- -(y * x %*% raw.coefficients) + log1pexp(x %*% raw.coefficients)
        classification   <- factor(as.numeric(c(1:length(residuals)) %in% object$best))
        ind              <- factor(as.numeric(!c(1:length(residuals)) %in% object$best))
        fitted.values    <- x %*% raw.coefficients
        names(residuals) <- 1:length(residuals)
        residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                       llim=residuals,ulim=residuals)
      }
      d    <- data.frame(X = drop(fitted.values), Y = residuals$residuals,
                         cls = factor(y, levels = c(0, 1)), flag = ind == "1")
      plot <- .classPlot(d, 2, if (vers == "raw") "not in the best subset" else "flagged as outliers", under = list(geom_hline(yintercept=0, linetype="dashed", color="gray")), zoom = zoom) +
        labs(title = mainfit, x = xlab, y = ylab)
      print(plot); plots <- c(plots, list(plot))

    } else if (family=="gaussian"){
      if (vers=="reweighted"){
        mainfit          <- expression(bold(paste("Standardized residuals vs ", X*hat(beta), ", linear regression")))
        xlab             <- expression(X*hat(beta))
        ylab             <- "Standardized residuals"
        lab.val          <- c("regular observation","outliers")
        classification   <- factor(as.numeric(!object$wt==0))
        ind              <- factor(as.numeric(object$wt==0))
        residuals        <- y - x %*% coefficients
        wg               <- weight.gaussian(residuals, object$raw.wt == 1, del)
        residuals        <- (residuals - wg$mu) / wg$s
        fitted.values    <- x %*% coefficients
        names(residuals) <- 1:length(residuals)
        residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                       llim=residuals,ulim=residuals)
      } else if (vers=="raw"){
        mainfit          <- expression(bold(paste("Standardized residuals vs ", X*hat(beta)[raw], ", linear regression")))
        xlab             <- expression(X*hat(beta)[raw])
        ylab             <- "Standardized raw residuals"
        lab.val          <- c("best subset","outliers")
        residuals        <- y - x %*% raw.coefficients
        wg               <- weight.gaussian(residuals, object$best, del)
        residuals        <- (residuals - wg$mu) / wg$s
        classification   <- factor(as.numeric(c(1:length(residuals)) %in% object$best))
        ind              <- factor(as.numeric(!c(1:length(residuals)) %in% object$best))
        fitted.values    <- x %*% raw.coefficients
        names(residuals) <- 1:length(residuals)
        residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                       llim=residuals,ulim=residuals)
      }
      d    <- data.frame(X = drop(fitted.values), Y = residuals$residuals, flag = ind == "1")
      plot <- .diagPlot(d, lab.val, under = list(geom_hline(yintercept=0, linetype="dashed", color="gray"),
                                                 geom_hline(yintercept=c(-cutoff, cutoff), color="gray60")),
                        zoom = zoom) +
        labs(title = mainfit, x = xlab, y = ylab)
      print(plot); plots <- c(plots, list(plot))
    }

  }
  invisible(plots)
}
