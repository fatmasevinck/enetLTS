
plotDiagnostic.enetLTS <- function(object,vers=c("reweighted","raw"),zoom=TRUE,...){

  vers <- match.arg(vers)
  family <- object$inputs$family

  if (family=="multinomial"){
    y <- object$inputs$y_factor
  } else {
    y <- object$inputs$yy
  }
  x <- object$inputs$xx

  coefficients     <- object$coefficients
  raw.coefficients <- object$raw.coefficients

  if (family=="multinomial"){
    if (vers=="reweighted"){
      mainfit           <- paste("Principal components of the linear predictors (reweighted fit)")
      xlab              <- "Component 1"
      ylab              <- "Component 2"
      lab.val           <- c("regular observation","outliers")
      classification    <- factor(as.numeric(object$wt==0))
      ind               <- factor(as.numeric(object$wt==0))
      fitted.values     <- zj <- x %*% coefficients
      pc.scores         <- princomp(zj)$sco
      pc.scores         <- as.data.frame(pc.scores[,1:2])
    } else if (vers=="raw"){
      mainfit           <- paste("Principal components of the linear predictors (raw fit)")
      xlab              <- "Component 1"
      ylab              <- "Component 2"
      lab.val           <- c("best subset","outliers")
      classification    <- factor(as.numeric(!c(1:nrow(x)) %in% object$best))
      ind               <- factor(as.numeric(!c(1:nrow(x)) %in% object$best))
      fitted.values     <- zj <- x %*% raw.coefficients
      pc.scores         <- princomp(zj)$sco
      pc.scores         <- as.data.frame(pc.scores[,1:2])
    }

    d    <- data.frame(X = pc.scores[, 1], Y = pc.scores[, 2], cls = factor(y), flag = ind == "1")
    plot <- .classPlot(d, nlevels(factor(y)),
                       if (vers == "raw") "not in the best subset" else "flagged as outliers", zoom = zoom) +
      labs(title = mainfit, x = xlab, y = ylab)
    print(plot)
  }

  if (family=="binomial"){
    if (vers=="reweighted"){
      mainfit              <- expression(bold(paste("y vs ", X*hat(beta), ", logistic regression")))
      xlab                 <- expression(X*hat(beta))
      ylab                 <- "y"
      lab.val              <- c("regular observation","outliers")
      classification       <- factor(as.numeric(!object$wt==0))
      ind                  <- factor(as.numeric(object$wt==0))
      fitted.values        <- x %*% coefficients
      probs                <- drop(1/(1+exp(-fitted.values))); probs[ind==1]=NA
      names(fitted.values) <- 1:length(fitted.values)
      plotfitted.values    <- data.frame(fitted.values=fitted.values,
                                         nam=names(fitted.values),
                                         llim=fitted.values,
                                         ulim=fitted.values)
    } else if (vers=="raw"){
      mainfit              <- expression(bold(paste("y vs ", X*hat(beta)[raw], ", logistic regression")))
      xlab                 <- expression(X*hat(beta)[raw])
      ylab                 <- "y"
      lab.val              <- c("best subset","outliers")
      fitted.values        <- x %*% raw.coefficients
      classification       <- factor(as.numeric(c(1:length(fitted.values)) %in% object$best))
      ind                  <- factor(as.numeric(!c(1:length(fitted.values)) %in% object$best))
      probs                <- drop(1/(1+exp(-fitted.values))); probs[ind==1]=NA
      names(fitted.values) <- 1:length(fitted.values)
      plotfitted.values    <- data.frame(fitted.values=fitted.values,
                                         nam=names(fitted.values),
                                         llim=fitted.values,
                                         ulim=fitted.values)
    }
    d     <- data.frame(X = drop(fitted.values), Y = y, cls = factor(y, levels = c(0, 1)), flag = ind == "1")
    plot  <- .classPlot(d, 2, if (vers == "raw") "not in the best subset" else "flagged as outliers",
                        under = list(geom_hline(yintercept = c(0, 1), linetype = "dashed", color = "gray")),
                        zy = FALSE, zoom = zoom, curveFun = function(t) 1 / (1 + exp(-t))) +
      labs(title = mainfit, x = xlab, y = ylab)
    print(plot)

  } else if (family=="gaussian"){
    if (vers=="reweighted"){

      mainfit              <- expression(bold(paste("y vs ", X*hat(beta), ", linear regression")))
      xlab                 <- expression(X*hat(beta))
      ylab                 <- "y"
      lab.val              <- c("regular observation","outliers")
      classification       <- factor(as.numeric(!object$wt==0))
      ind                  <- factor(as.numeric(object$wt==0))
      fitted.values        <- x %*% coefficients
      names(fitted.values) <- 1:length(fitted.values)
      plotfitted.values    <- data.frame(fitted.values=fitted.values,
                                         nam=names(fitted.values),
                                         llim=fitted.values,
                                         ulim=fitted.values)
    } else if (vers=="raw"){
      mainfit              <- expression(bold(paste("y vs ", X*hat(beta)[raw], ", linear regression")))
      xlab                 <- expression(X*hat(beta)[raw])
      ylab                 <- "y"
      lab.val              <- c("best subset","outliers")
      fitted.values        <- x %*% raw.coefficients
      classification       <- factor(as.numeric(c(1:length(fitted.values)) %in% object$best))
      ind                  <- factor(as.numeric(!c(1:length(fitted.values)) %in% object$best))
      names(fitted.values) <- 1:length(fitted.values)
      plotfitted.values    <- data.frame(fitted.values=fitted.values,
                                         nam=names(fitted.values),
                                         llim=fitted.values,
                                         ulim=fitted.values)
    }
    d    <- data.frame(X = drop(fitted.values), Y = y, flag = ind == "1")
    plot <- .diagPlot(d, lab.val, under = list(geom_abline(intercept = 0, slope = 1, color = "gray")),
                      zoom = zoom) +
      labs(title = mainfit, x = xlab, y = ylab)
    print(plot)
  }
  invisible(plot)
}
