
residuals.enetLTS <-
  function(object,vers=c("reweighted","raw","both"),...){
    
    vers <- match.arg(vers)
    
    if (vers=="reweighted"){
      if (object$inputs$family=="multinomial"){
        yhat      <- object$inputs$x %*% object$coefficients
        probs     <- exp(yhat)/apply(exp(yhat),1,sum)
        deviances <- (-apply(object$inputs$y*log(probs),1,sum))
        deviances[is.nan(deviances)] <- 0
        residuals <- deviances
      } else if (object$inputs$family=="binomial"){
        residuals <-  -(object$inputs$y * object$inputs$x %*% object$coefficients) +
          log(1+exp(object$inputs$x %*% c(object$a0,object$coefficients)))
      } else if (object$inputs$family=="gaussian"){
        residuals <- object$inputs$y - object$inputs$x %*% object$coefficients
      }
      nfit <- list(residuals=residuals)
    } else if (vers=="raw"){
      if (object$inputs$family=="multinomial"){
        raw.yhat      <- object$inputs$x %*% object$raw.coefficients
        raw.probs     <- exp(raw.yhat)/apply(exp(raw.yhat),1,sum)
        raw.deviances <- (-apply(object$inputs$y*log(raw.probs),1,sum))
        raw.deviances[is.nan(raw.deviances)] <- 0
        raw.residuals <- raw.deviances
      } else if (object$inputs$family=="binomial"){
        raw.residuals <- -(object$inputs$y * object$inputs$x %*% object$raw.coefficients) +
          log(1+exp(object$inputs$x %*% object$raw.coefficients))
      } else if (object$inputs$family=="gaussian"){
        raw.residuals <- object$inputs$y - object$inputs$x %*% object$raw.coefficients
      }
      nfit <- list(raw.residuals=raw.residuals)
    } else if (vers=="both"){
      if (object$inputs$family=="multinomial"){
        yhat      <- object$inputs$x %*% object$coefficients
        probs     <- exp(yhat)/apply(exp(yhat),1,sum)
        deviances <- (-apply(object$inputs$y*log(probs),1,sum))
        deviances[is.nan(deviances)] <- 0
        residuals <- deviances
        
        raw.yhat      <- object$inputs$x %*% object$raw.coefficients
        raw.probs     <- exp(raw.yhat)/apply(exp(raw.yhat),1,sum)
        raw.deviances <- (-apply(object$inputs$y*log(raw.probs),1,sum))
        raw.deviances[is.nan(raw.deviances)] <- 0
        raw.residuals <- raw.deviances
        
      } else if (object$inputs$family=="binomial"){
        residuals     <-  -(object$inputs$y * object$inputs$x %*% c(object$a0,object$coefficients)) +
          log(1+exp(object$inputs$x %*% c(object$a0,object$coefficients)))
        raw.residuals <- -(object$inputs$y * object$inputs$x %*% c(object$a0,object$raw.coefficients)) +
          log(1+exp(object$inputs$x %*% c(object$a0,object$raw.coefficients)))
      } else if (object$inputs$family=="gaussian"){
        residuals     <- object$inputs$y - object$inputs$x %*% c(object$a0,object$coefficients)
        raw.residuals <- object$inputs$y - object$inputs$x %*% c(object$a0,object$raw.coefficients)
      }
      nfit <- list(residuals=residuals,raw.residuals=raw.residuals)
    }
    nfit
  }
