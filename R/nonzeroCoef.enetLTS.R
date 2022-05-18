
nonzeroCoef.enetLTS  <- function (object,vers=c("reweighted","raw"))
{
  vers <- match.arg(vers)
  
  if(vers=="reweighted"){
    beta <- object$coefficients
  } else if (vers=="raw"){
    beta <- object$raw.coefficients
  }
  
  family <- object$inputs$family
  
  if (family=="multinomial"){
    nr    <- nrow(beta)
    beta  <- abs(beta)>0 # this is sparse
    which <- seq(nr)
    ones  <- rep(1,ncol(beta))
    nz    <- as.vector((beta%*%ones)>0)
    which <- which[nz]
    beta  <- as.matrix(beta[which,,drop=FALSE])
    nzel  <- function(x, which) if (any(x)) 
      which[x]
    else NULL
    which <- apply(beta, 2, nzel, which)
    if(!is.list(which)) which = data.frame(which)# apply can return a matrix!!
    beta  <- which
    
  } else if (family=="binomial" | family=="gaussian"){
    
    beta <- as.matrix(drop(beta))
    beta <- abs(beta)>0      # this is sparse
    beta <- which(beta)
    # names(beta) <- 1:length(beta)
  }
  return(beta)
}

