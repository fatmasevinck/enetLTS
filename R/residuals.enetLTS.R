
residuals.enetLTS <-
  function(object,vers=c("reweighted","raw","both"),...){

    vers <- match.arg(vers)
    coefficients <- object$coefficients
    raw.coefficients <- object$raw.coefficients
    xx <- object$inputs$xx
    yy <- object$inputs$yy


    if (vers=="reweighted"){
      if (object$inputs$family=="multinomial"){
        yhat      <- xx %*% coefficients
        probs     <- exp(yhat)/apply(exp(yhat),1,sum)
        deviances <- (-apply(yy*log(probs),1,sum))
        deviances[is.nan(deviances)] <- 0
        residuals <- deviances
      } else if (object$inputs$family=="binomial"){
        residuals <-  -(yy * xx %*% coefficients) +
          log(1+exp(xx %*% coefficients))
      } else if (object$inputs$family=="gaussian"){
        residuals <- yy - xx %*% coefficients
      }
      nfit <- list(residuals=drop(residuals))
    } else if (vers=="raw"){
      if (object$inputs$family=="multinomial"){
        raw.yhat      <- xx %*% raw.coefficients
        raw.probs     <- exp(raw.yhat)/apply(exp(raw.yhat),1,sum)
        raw.deviances <- (-apply(yy*log(raw.probs),1,sum))
        raw.deviances[is.nan(raw.deviances)] <- 0
        raw.residuals <- raw.deviances
      } else if (object$inputs$family=="binomial"){
        raw.residuals <- -(yy * xx %*% raw.coefficients) +
          log(1+exp(xx %*% raw.coefficients))
      } else if (object$inputs$family=="gaussian"){
        raw.residuals <- yy - xx %*% raw.coefficients
      }
      nfit <- list(raw.residuals=drop(raw.residuals))
    } else if (vers=="both"){
      if (object$inputs$family=="multinomial"){
        yhat      <- xx %*% coefficients
        probs     <- exp(yhat)/apply(exp(yhat),1,sum)
        deviances <- (-apply(yy*log(probs),1,sum))
        deviances[is.nan(deviances)] <- 0
        residuals <- deviances

        raw.yhat      <- xx %*% raw.coefficients
        raw.probs     <- exp(raw.yhat)/apply(exp(raw.yhat),1,sum)
        raw.deviances <- (-apply(yy*log(raw.probs),1,sum))
        raw.deviances[is.nan(raw.deviances)] <- 0
        raw.residuals <- raw.deviances

      } else if (object$inputs$family=="binomial"){
        residuals     <-  -(yy * xx %*% coefficients) +
          log(1+exp(object$inputs$x %*% coefficients))
        raw.residuals <- -(yy * xx %*% raw.coefficients) +
          log(1+exp(xx %*% raw.coefficients))
      } else if (object$inputs$family=="gaussian"){
        residuals     <- yy - xx %*% coefficients
        raw.residuals <- yy - xx %*% raw.coefficients
      }
      nfit <- list(residuals=drop(residuals),raw.residuals=drop(raw.residuals))
    }
    nfit
  }
