
predict.enetLTS <-
  function(object,newX,vers=c("reweighted","raw"),
           type=c("link","response","coefficients","nonzero","class"),...)
  {
    
    type <- match.arg(type)
    vers <- match.arg(vers)
    if (type=="class" & object$inputs$family=="gaussian") stop("'class' is not  available for 'gaussian' family")
    
    if (missing(newX)){
      if(!match(type,c("coefficients","nonzero"),FALSE)) stop("You need to supply a value for 'newX'")
    }
    
    if (vers=="reweighted"){
      coefficients     <- object$coefficients
    } else if (vers=="raw"){
      raw.coefficients <- object$raw.coefficients
    } 
    
    if (object$inputs$intercept=="TRUE"){
      newX <- cbind(1,newX)
    }
    
    # link
    if (type=="link"){
      if (object$inputs$family=="multinomial"){
        if (vers=="reweighted"){
          yhat               <- as.matrix(newX %*% coefficients)
          fit.link           <- drop(yhat)
          colnames(fit.link) <- paste0("class", 1:(length(object$classize))) 
          rownames(fit.link) <- 1:nrow(newX)
          link <- fit.link
        } else if (vers=="raw"){
          raw.yhat               <- as.matrix(newX %*% raw.coefficients)
          raw.fit.link           <- drop(raw.yhat)
          colnames(raw.fit.link) <- paste0("class", 1:(length(object$classize))) 
          rownames(raw.fit.link) <- 1:nrow(newX)
          link <- raw.fit.link
        } 
      } else if (object$inputs$family=="gaussian" | object$inputs$family=="binomial"){
        if (vers=="reweighted"){
          yhat         <- as.matrix(newX %*% coefficients)
          fit.link     <- drop(yhat)
          link         <- fit.link
        } else if (vers=="raw"){
          raw.yhat     <- as.matrix(newX %*% raw.coefficients)
          raw.fit.link <- drop(raw.yhat)
          link         <- raw.fit.link
        } 
      }
      return(link)
    }    
    
    # response
    if (type=="response"){
      if (object$inputs$family=="multinomial"){
        if (vers=="reweighted"){
          yhat            <- newX %*% coefficients
          probs           <- exp(yhat)/apply(exp(yhat),1,sum)
          colnames(probs) <- paste0("class", 1:(length(object$classize))) 
          rownames(probs) <- 1:nrow(newX)
          response    <- probs
        } else if (vers=="raw"){
          raw.yhat            <- newX %*% raw.coefficients
          raw.probs           <- exp(raw.yhat)/apply(exp(raw.yhat),1,sum)
          colnames(raw.probs) <- paste0("class", 1:(length(object$classize)))  
          rownames(raw.probs) <- 1:nrow(newX)
          response <- raw.probs
        } 
      } else if (object$inputs$family=="binomial"){
        if (vers=="reweighted"){
          probs        <- drop(1/(1+exp(-newX%*%coefficients)))
          response <- probs
        } else if (vers=="raw"){
          raw.probs    <- drop(1/(1+exp(-newX%*%raw.coefficients)))
          response     <- raw.probs
        } 
      } else if (object$inputs$family=="gaussian"){
        if (vers=="reweighted"){
          res          <- drop(as.matrix(newX%*%coefficients))
          response   <- res
        } else if (vers=="raw"){
          res          <- drop(as.matrix(newX%*%raw.coefficients))
          response    <- res
        } 
      }
      return(response)
    }    
    
    
    # coefff
    if (type=="coefficients"){
      if (object$inputs$family=="multinomial"){
        if (vers=="reweighted"){
          colnames(coefficients) <- paste0("class", 1:(length(object$classize))) 
          rownames(coefficients) <- 1:nrow(coefficients)
          coefficients=coefficients
        } else if (vers=="raw"){
          colnames(raw.coefficients) <- paste0("class", 1:(length(object$classize))) 
          rownames(raw.coefficients) <- 1:nrow(raw.coefficients)
          coefficients=raw.coefficients
        } 
      } else if (object$inputs$family=="gaussian" | object$inputs$family=="binomial"){
        if (vers=="reweighted"){
          coefficients=coefficients
        } else if (vers=="raw"){
          coefficients=raw.coefficients
        } 
      }
      return(coefficients)
    }
    
    # nonzero
    if (type=="nonzero"){
      if (vers=="reweighted"){
        reweighted.nonzeroCoef <- nonzeroCoef.enetLTS(object,vers="reweighted")
        nonzeroCoef=reweighted.nonzeroCoef
      } else if (vers=="raw"){
        raw.nonzeroCoef <- nonzeroCoef.enetLTS(object,vers="raw")
        nonzeroCoef=raw.nonzeroCoef
      } 
      return(nonzeroCoef)
    }
    
    # class
    if (type=="class"){
      if (object$inputs$family=="multinomial"){
        
        if (vers=="reweighted"){
          yhat             <- newX %*% coefficients
          probs            <- exp(yhat)/apply(exp(yhat),1,sum)
          reweighted.class <- apply(probs,1,which.max)
          class            <-  reweighted.class
        } else if (vers=="raw"){
          raw.yhat  <- newX %*% raw.coefficients
          raw.probs <- exp(raw.yhat)/apply(exp(raw.yhat),1,sum)
          raw.class <- apply(raw.probs,1,which.max)
          class     <- raw.class
        } 
        return(class)
      } else if (object$inputs$family=="binomial"){
        if (vers=="reweighted"){
          res              <- newX %*% coefficients
          cnum             <- ifelse(res>0.5,2,1)
          reweighted.class <- object$classnames[cnum]
          class            <- reweighted.class
        } else if (vers=="raw"){
          res       <- newX %*% raw.coefficients
          cnum      <- ifelse(res>0.5,2,1)
          raw.class <- object$classnames[cnum]
          class     <- raw.class
        }
        return(class)
      }
    }
    
  }
