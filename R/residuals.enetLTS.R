
residuals.enetLTS <-
  function(object,vers=c("reweighted","raw"),...){
    
    vers <- match.arg(vers)
    
    if (vers=="reweighted"){
      if (object$inputs$family=="multinomial"){
        yhat      <- object$inputs$xx %*% object$coefficients
        probs     <- exp(yhat)/apply(exp(yhat),1,sum)
        deviances <- (-apply(object$inputs$yy*log(probs),1,sum))
        deviances[is.nan(deviances)] <- 0
        residuals <- deviances
      } else if (object$inputs$family=="binomial"){
        residuals <-  -(object$inputs$yy * object$inputs$xx %*% object$coefficients) +
          log(1+exp(object$inputs$xx %*% c(object$a0,object$coefficients)))
      } else if (object$inputs$family=="gaussian"){
        residuals <- object$inputs$yy - object$inputs$xx %*% object$coefficients
      }
      return(drop(residuals))
    }  
    if (vers=="raw"){
      if (object$inputs$family=="multinomial"){
        raw.yhat      <- object$inputs$xx %*% object$raw.coefficients
        raw.probs     <- exp(raw.yhat)/apply(exp(raw.yhat),1,sum)
        raw.deviances <- (-apply(object$inputs$yy*log(raw.probs),1,sum))
        raw.deviances[is.nan(raw.deviances)] <- 0
        raw.residuals <- raw.deviances
      } else if (object$inputs$family=="binomial"){
        raw.residuals <- -(object$inputs$yy * object$inputs$xx %*% object$raw.coefficients) +
          log(1+exp(object$inputs$xx %*% object$raw.coefficients))
      } else if (object$inputs$family=="gaussian"){
        raw.residuals <- object$inputs$yy - object$inputs$xx %*% object$raw.coefficients
      }
      return(drop(raw.residuals))
    } 
  }

