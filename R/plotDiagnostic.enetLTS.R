
plotDiagnostic.enetLTS <- function(object,vers=c("reweighted","raw"),...){
  
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
      mainfit           <- paste("First two components of estimated scores for multinomial logistic regression (reweighted)")
      xlab              <- "Component 1"
      ylab              <- "Component 2"
      lab.val           <- c("regular observation","outliers")
      classification    <- factor(as.numeric(object$raw.wt==0))
      ind               <- factor(as.numeric(object$raw.wt==0))
      fitted.values     <- zj <- x %*% coefficients
      pc.scores         <- princomp(zj)$sco
      pc.scores         <- as.data.frame(pc.scores[,1:2])
    } else if (vers=="raw"){
      mainfit           <- paste("First two components of estimated scores for multinomial logistic regression (raw)")
      xlab              <- "Component 1"
      ylab              <- "Component 2"
      lab.val           <- c("best subset","outliers")
      classification    <- factor(as.numeric(!c(1:nrow(x)) %in% object$best))
      ind               <- factor(as.numeric(!c(1:length(residuals)) %in% object$best))
      fitted.values     <- zj <- x %*% raw.coefficients
      pc.scores         <- princomp(zj)$sco
      pc.scores         <- as.data.frame(pc.scores[,1:2])
    }
    
    plot <- ggplot() +
      geom_point(data = pc.scores, aes(x=Comp.1, y=Comp.2, colour = factor(y),
                                       shape=classification), size=2) +
      labs(title = mainfit, x = xlab, y = ylab) +
      theme(plot.title=element_text(size=rel(1.3)),
            axis.title   = element_text(size = 13),
            axis.text.y  = element_text(size = 13),
            axis.text.x  = element_text(size = 13),
            legend.text  = element_text(size = 13),
            legend.title = element_text(size = 13),
            legend.position = "bottom") +
      guides(color=guide_legend(ncol=2, byrow=TRUE)) +
      scale_colour_discrete(name= "Groups") +
      scale_shape_manual(values=c(4, 17),name="Diagnostics",breaks=c("0", "1"), labels=lab.val)
    print(plot)
  }
  
  if (family=="binomial"){
    if (vers=="reweighted"){
      mainfit              <- expression(paste("y vs ", X*hat(beta), " for Binary Logistic Regression"))
      xlab                 <- expression(X*hat(beta))
      ylab                 <- "y"
      lab.val              <- c("regular observation","outliers")
      classification       <- factor(as.numeric(!object$raw.wt==0))
      ind                  <- factor(as.numeric(object$raw.wt==0))
      fitted.values        <- x %*% coefficients
      probs                <- drop(1/(1+exp(-fitted.values))); probs[ind==1]=NA
      names(fitted.values) <- 1:length(fitted.values)
      plotfitted.values    <- data.frame(fitted.values=fitted.values,
                                         nam=names(fitted.values),
                                         llim=fitted.values,
                                         ulim=fitted.values)
    } else if (vers=="raw"){
      mainfit              <- expression(paste("y vs ", X*hat(beta)[raw], " for Binary Logistic Regression"))
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
    plot <- ggplot(plotfitted.values, aes(x=fitted.values,
                                          y=y,
                                          color=factor(y),
                                          shape=ind)) +
      geom_point(size=2) +
      geom_line(aes(y=probs,x=fitted.values), linetype="dotted", color="black") +
      geom_hline(yintercept=1, linetype="dashed", color="gray") +
      geom_hline(yintercept=0, linetype="dashed", color="gray") +
      labs(title = mainfit, x = xlab, y = ylab) +
      theme(plot.title   = element_text(size=rel(1.3)),
            axis.title   = element_text(size = 13),
            axis.text.y  = element_text(size = 13),
            axis.text.x  = element_text(size = 13),
            legend.text  = element_text(size = 13),
            legend.title = element_text(size = 13),
            legend.position = "bottom") +
      guides(color=guide_legend(nrow=2, byrow=TRUE)) +
      scale_colour_discrete(name= "Groups") +
      # scale_colour_discrete(name = "Diagnostics", breaks = c("1", "0"), labels = lab.val) +
      scale_shape_manual(values=c(4, 17),name="Diagnostics",breaks=c("0", "1"), labels=lab.val); options(warn = - 1)
    print(plot)
    
  } else if (family=="gaussian"){
    if (vers=="reweighted"){
      
      mainfit              <- expression(paste("y vs ", X*hat(beta), " for Linear Regression"))
      xlab                 <- expression(X*hat(beta))
      ylab                 <- "y"
      lab.val              <- c("regular observation","outliers")
      classification       <- factor(as.numeric(!object$raw.wt==0))
      ind                  <- factor(as.numeric(object$raw.wt==0))
      fitted.values        <- x %*% coefficients
      names(fitted.values) <- 1:length(fitted.values)
      plotfitted.values    <- data.frame(fitted.values=fitted.values,
                                         nam=names(fitted.values),
                                         llim=fitted.values,
                                         ulim=fitted.values)
    } else if (vers=="raw"){
      mainfit              <- expression(paste("y vs ", X*hat(beta)[raw], " for Linear Regression"))
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
    plot <- ggplot(plotfitted.values, aes(x=fitted.values,
                                          y=y,
                                          color=classification,
                                          shape=ind)) +
      geom_point(size=2) +
      # geom_hline(yintercept=0, linetype="dashed",color="gray") +
      geom_abline(intercept = 0, slope = 1, color="gray", size = 0.5) +
      labs(title = mainfit, x = xlab, y = ylab) +
      theme(plot.title   = element_text(size=rel(1.3)),
            axis.title   = element_text(size = 13),
            axis.text.y  = element_text(size = 13),
            axis.text.x  = element_text(size = 13),
            legend.text  = element_text(size = 13),
            legend.title = element_text(size = 13),
            legend.position= "bottom") +
      scale_colour_discrete(name="Diagnostics",breaks=c("1", "0"), labels=lab.val) +
      scale_shape_manual(values=c(4, 17),name="Diagnostics",breaks=c("0", "1"), labels=lab.val)
    print(plot)
  }
  
}
