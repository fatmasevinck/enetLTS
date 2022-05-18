
plotResid.enetLTS <- function(object,vers=c("reweighted","raw"),...){

  vers   <- match.arg(vers)
  family <- object$inputs$family


  if (family=="multinomial"){
    y.factor <- object$inputs$y_factor
  }

  x <- object$inputs$xx
  y <- object$inputs$yy

  coefficients     <- object$coefficients
  raw.coefficients <- object$raw.coefficients

  if (family=="multinomial"){
    if (vers=="reweighted"){
      mainfit          <- paste("Deviances vs Index for Multinomial Logistic Regression")
      xlab             <- "Index"
      ylab             <- "Deviances"
      lab.val          <- c("regular observation","outliers")
      classification   <- factor(as.numeric(!object$raw.wt==0))
      ind              <- factor(as.numeric(object$raw.wt==0))
      yhat             <- x %*% coefficients
      probs            <- exp(yhat)/apply(exp(yhat),1,sum)
      residuals        <- (-apply(y*log(probs),1,sum))
      Index            <- as.numeric(c(1:length(residuals)))
      names(residuals) <- 1:length(residuals)
      residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                     llim=residuals,ulim=residuals)
    } else if (vers=="raw"){
      mainfit          <- paste("Raw Deviances vs Index for Multinomial Logistic Regression")
      xlab             <- "Index"
      ylab             <- "Raw deviances"
      lab.val          <- c("regular observation","outliers")
      yhat             <- x %*% raw.coefficients
      probs            <- exp(yhat)/apply(exp(yhat),1,sum)
      residuals        <- (-apply(y*log(probs),1,sum))
      classification   <- factor(as.numeric(c(1:length(residuals)) %in% object$best))
      ind              <- factor(as.numeric(!c(1:length(residuals)) %in% object$best))
      Index            <- as.numeric(c(1:length(residuals)))
      names(residuals) <- 1:length(residuals)
      residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                     llim=residuals,ulim=residuals)
    }
    plot <- ggplot(residuals, aes(x=Index,
                                  y=residuals,
                                  color=factor(y.factor),
                                  shape=ind)) +
      geom_point(size=2) +
      geom_hline(yintercept=0, linetype="dashed",color="gray") +
      labs(title = mainfit, x = xlab, y = ylab) +
      theme(plot.title   = element_text(size=rel(1.3),face="bold"),
            axis.title   = element_text(size = 13),
            axis.text.y  = element_text(size = 13),
            axis.text.x  = element_text(size = 13),
            legend.text  = element_text(size = 13),
            legend.title = element_text(size = 13)) +
      scale_colour_discrete(name= "Groups") +
      # scale_shape_discrete(name = "Diagnostics", breaks=c("0", "1"), labels = lab.val)
      scale_shape_manual(values=c(4, 17),name="Diagnostics",breaks=c("0", "1"), labels=lab.val)
    print(plot)
  } else if (family=="binomial"){
    if (vers=="reweighted"){
      mainfit          <- paste("Deviances vs Index for Binary Logistic Regression")
      xlab             <- "Index"
      ylab             <- "Deviances"
      lab.val          <- c("regular observation","outliers")
      classification   <- factor(as.numeric(!object$raw.wt==0))
      ind              <- factor(as.numeric(object$raw.wt==0))
      residuals        <- -(y * x %*% coefficients) + log(1+exp(x %*% coefficients))
      Index            <- as.numeric(c(1:length(residuals)))
      names(residuals) <- 1:length(residuals)
      residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                     llim=residuals,ulim=residuals)
    } else if (vers=="raw"){
      mainfit          <- paste("Raw Deviances vs Index for Binary Logistic Regression")
      xlab             <- "Index"
      ylab             <- "Raw deviances"
      lab.val          <- c("best subset","outliers")
      residuals        <- -(y * x %*% raw.coefficients) + log(1+exp(x %*% raw.coefficients))
      classification   <- factor(as.numeric(c(1:length(residuals)) %in% object$best))
      ind              <- factor(as.numeric(!c(1:length(residuals)) %in% object$best))
      Index            <- as.numeric(c(1:length(residuals)))
      names(residuals) <- 1:length(residuals)
      residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                     llim=residuals,ulim=residuals)
    }
    plot <- ggplot(residuals, aes(x=Index,
                                  y=residuals,
                                  color=classification,
                                  shape=ind)) +
      geom_point(size=2) +
      geom_hline(yintercept=0, linetype="dashed",color="gray") +
      labs(title = mainfit, x = xlab, y = ylab) +
      theme(plot.title   = element_text(size=rel(1.3),face="bold"),
            axis.title   = element_text(size = 13),
            axis.text.y  = element_text(size = 13),
            axis.text.x  = element_text(size = 13),
            legend.text  = element_text(size = 13),
            legend.title = element_text(size = 13)) +
      scale_colour_discrete(name="Diagnostics",breaks=c("1", "0"),labels=lab.val) +
      # scale_shape_discrete(name="Diagnostics",breaks=c("0", "1"),labels=lab.val)
      scale_shape_manual(values=c(4, 17),name="Diagnostics",breaks=c("0", "1"), labels=lab.val)
    print(plot)
  } else if (family=="gaussian"){
    if (vers=="reweighted"){
      mainfit          <- paste("Standardized Residuals vs Index for Regression")
      xlab             <- "Index"
      ylab             <- "Residuals"
      lab.val          <- c("regular observation","outliers")
      classification   <- factor(as.numeric(!object$raw.wt==0))
      ind              <- factor(as.numeric(object$raw.wt==0))
      residuals        <- y - x %*% coefficients
      residuals        <- residuals/mad(residuals)
      Index            <- as.numeric(c(1:length(residuals)))
      names(residuals) <- 1:length(residuals)
      residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                     llim=residuals,ulim=residuals)
    } else if (vers=="raw"){
      mainfit          <- paste("Standardized Raw Residuals vs Index for Regression")
      xlab             <- "Index"
      ylab             <- "Raw residuals"
      lab.val          <- c("best subset","outliers")
      residuals        <- y - x %*% raw.coefficients
      residuals        <- residuals/mad(residuals)
      classification   <- factor(as.numeric(c(1:length(residuals)) %in% object$best))
      ind              <- factor(as.numeric(!c(1:length(residuals)) %in% object$best))
      Index            <- as.numeric(c(1:length(residuals)))
      names(residuals) <- 1:length(residuals)
      residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                     llim=residuals,ulim=residuals)
    }
    plot <- ggplot(residuals, aes(x=Index,
                                  y=residuals,
                                  color=classification,
                                  shape=ind)) +
      geom_point(size=2) +
      geom_hline(yintercept=0, linetype="dashed",color="gray") +
      geom_hline(yintercept=-2.5,linetype="solid",color="gray") +
      geom_hline(yintercept=2.5,linetype="solid",color="gray") +
      labs(title = mainfit, x = xlab, y = ylab) +
      theme(plot.title   = element_text(size=rel(1.3),face="bold"),
            axis.title   = element_text(size = 13),
            axis.text.y  = element_text(size = 13),
            axis.text.x  = element_text(size = 13),
            legend.text  = element_text(size = 13),
            legend.title = element_text(size = 13)) +
      scale_colour_discrete(name="Diagnostics",breaks=c("1", "0"), labels=lab.val) +
      # scale_shape_discrete(name="Diagnostics",breaks=c("0", "1"), labels=lab.val)
      scale_shape_manual(values=c(4, 17),name="Diagnostics",breaks=c("0", "1"), labels=lab.val)
    print(plot)
  }


  if (family=="binomial" | family=="gaussian"){
    if (family=="binomial"){
      if (vers=="reweighted"){
        mainfit          <- expression(paste("Deviances vs ", X*hat(beta), " for Binary Logistic Regression"))
        xlab             <- expression(X*hat(beta))
        ylab             <- "Deviances"
        lab.val          <- c("regular observation","outliers")
        classification   <- factor(as.numeric(!object$raw.wt==0))
        ind              <- factor(as.numeric(object$raw.wt==0))
        residuals        <- -(y * x %*% coefficients) + log(1+exp(x %*% coefficients))
        fitted.values    <- x %*% coefficients
        names(residuals) <- 1:length(residuals)
        residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                       llim=residuals,ulim=residuals)
      } else if (vers=="raw"){
        mainfit          <- expression(paste("Raw Deviances vs ", X*hat(beta)[raw], " for Binary Logistic Regression"))
        xlab             <- expression(X*hat(beta)[raw])
        ylab             <- "Raw deviances"
        lab.val          <- c("best subset","outliers")
        residuals        <- -(y * x %*% raw.coefficients) + log(1+exp(x %*% raw.coefficients))
        classification   <- factor(as.numeric(c(1:length(residuals)) %in% object$best))
        ind              <- factor(as.numeric(!c(1:length(residuals)) %in% object$best))
        fitted.values    <- x %*% raw.coefficients
        names(residuals) <- 1:length(residuals)
        residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                       llim=residuals,ulim=residuals)
      }
      plot <- ggplot(residuals, aes(x=fitted.values,
                                    y=residuals,
                                    color=classification,
                                    shape=ind)) +
        geom_point(size=2) +
        geom_hline(yintercept=0, linetype="dashed",color="gray") +
        labs(title = mainfit, x = xlab, y = ylab) +
        theme(plot.title   = element_text(size=rel(1.3),face="bold"),
              axis.title   = element_text(size = 13),
              axis.text.y  = element_text(size = 13),
              axis.text.x  = element_text(size = 13),
              legend.text  = element_text(size = 13),
              legend.title = element_text(size = 13)) +
        scale_colour_discrete(name="Diagnostics",breaks=c("1", "0"),labels=lab.val) +
        scale_shape_manual(values=c(4, 17),name="Diagnostics",breaks=c("0", "1"), labels=lab.val)
      print(plot)

    } else if (family=="gaussian"){
      if (vers=="reweighted"){
        mainfit          <- expression(paste("Standardized Residuals vs ", X*hat(beta), " for Linear Regression"))
        xlab             <- expression(X*hat(beta))
        ylab             <- "Standarddized residuals"
        lab.val          <- c("regular observation","outliers")
        classification   <- factor(as.numeric(!object$raw.wt==0))
        ind              <- factor(as.numeric(object$raw.wt==0))
        residuals        <- y - x %*% coefficients
        residuals        <- residuals/mad(residuals)
        fitted.values    <- x %*% coefficients
        names(residuals) <- 1:length(residuals)
        residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                       llim=residuals,ulim=residuals)
      } else if (vers=="raw"){
        mainfit          <- expression(paste("Standardized Raw Residuals vs ", X*hat(beta)[raw], " for Linear Regression"))
        xlab             <- expression(X*hat(beta)[raw])
        ylab             <- "Stamdardized raw residuals"
        lab.val          <- c("best subset","outliers")
        residuals        <- y - x %*% raw.coefficients
        residuals        <- residuals/mad(residuals)
        classification   <- factor(as.numeric(c(1:length(residuals)) %in% object$best))
        ind              <- factor(as.numeric(!c(1:length(residuals)) %in% object$best))
        fitted.values    <- x %*% raw.coefficients
        names(residuals) <- 1:length(residuals)
        residuals        <- data.frame(residuals=residuals,nam=names(residuals),
                                       llim=residuals,ulim=residuals)
      }
      plot <- ggplot(residuals, aes(x=fitted.values,
                                    y=residuals,
                                    color=classification,
                                    shape=ind)) +
        geom_point(size=2) +
        geom_hline(yintercept=0, linetype="dashed",color="gray") +
        geom_hline(yintercept=-2.5,linetype="solid",color="gray") +
        geom_hline(yintercept=2.5,linetype="solid",color="gray") +
        labs(title = mainfit, x = xlab, y = ylab) +
        theme(plot.title   = element_text(size=rel(1.3),face="bold"),
              axis.title   = element_text(size = 13),
              axis.text.y  = element_text(size = 13),
              axis.text.x  = element_text(size = 13),
              legend.text  = element_text(size = 13),
              legend.title = element_text(size = 13)) +
        scale_colour_discrete(name="Diagnostics",breaks=c("1", "0"),labels=lab.val) +
        # scale_shape_discrete(name="Diagnostics",breaks=c("0", "1"),labels=lab.val)
        scale_shape_manual(values=c(4, 17),name="Diagnostics",breaks=c("0", "1"), labels=lab.val)
      print(plot)
    }

  }
}

