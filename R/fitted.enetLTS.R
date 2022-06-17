

fitted.enetLTS <-
  function(object,vers=c("reweighted","raw"),type=c("response","class"),...){
    
    vers <- match.arg(vers)
    type <- match.arg(type)
    
    if(type=="class" & object$inputs$family=="gaussian"){stop("'class' is only available for logistic regression")}
    
    reweighted.coefficients <- object$coefficients
    raw.coefficients        <- object$raw.coefficients
    
    xx <- object$inputs$xx
    yy <- object$inputs$yy

    if (object$inputs$family=="multinomial"){
      if (vers=="reweighted"){
        u <- xx %*% reweighted.coefficients
        if (type=="class"){
          probs         <- exp(u)/apply(exp(u),1,sum)
          fitted.values <- apply(probs,1,which.max)
        } else if (type=="response"){
          fitted.values           <- exp(u)/apply(exp(u),1,sum)
          colnames(fitted.values) <- paste0("class", 1:(length(object$classize))) 
          rownames(fitted.values) <- 1:nrow(xx)
        }
        nfit <- list(fitted.values=fitted.values)
      } else if (vers=="raw"){
        uu <- xx %*% raw.coefficients
        if (type=="class"){
          raw.probs         <- exp(uu)/apply(exp(uu),1,sum)
          raw.fitted.values <- apply(raw.probs,1,which.max)
        } else if (type=="response"){
          fitted.values           <- exp(uu)/apply(exp(uu),1,sum)
          colnames(fitted.values) <- paste0("class", 1:(length(object$classize))) 
          rownames(fitted.values) <- 1:nrow(xx)
        }
        nfit <- fitted.values
      } 
    } else if (object$inputs$family=="binomial"){
        if (vers=="reweighted"){
          u <- xx %*% reweighted.coefficients
          if (type=="class"){
            fitted.values <- ifelse(u>0.5,0,1)
          } else if (type=="response"){
            fitted.values <-  1/(1+exp(-u))
          }
        } else if (vers=="raw"){
          uu <- xx %*% raw.coefficients
          if (type=="class"){
            fitted.values <- ifelse(uu>0.5,0,1)
          } else if (type=="response"){
            fitted.values <-  1/(1+exp(-uu))
          }
          nfit <- fitted.values
        } 
      } else if (object$inputs$family=="gaussian"){
        if (vers=="reweighted"){
          res=as.matrix(xx %*% reweighted.coefficients)
          fitted.values <- res
        } else if (vers=="raw"){
          res=as.matrix(xx %*% raw.coefficients)
          fitted.values <- res
        } 
        nfit <- fitted.values
      }
    return(nfit)
  }



