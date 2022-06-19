

coef.enetLTS <-
  function(object ,vers=c("reweighted","raw"), zeros=TRUE,...)
  {
    vers=match.arg(vers)
    nbeta <- predict.enetLTS(object,newX=object$inputs$xx,vers=vers,type="coefficients",...)
    if (object$inputs$family=="gaussian" | object$inputs$family=="binomial") {
      nbeta <- as.numeric(unlist(nbeta))
    }
    if (isTRUE(zeros)) {
      nbeta <- nbeta
      names(nbeta) <- 1:length(nbeta)
    } else if (!isTRUE(zeros)) {
      nbeta <- nonzeroCoef.enetLTS(object,vers)
      # namesbeta <- which(nbeta != 0)
      # nbeta <- nbeta[nbeta != 0]
      # names(nbeta) <- namesbeta
    }
    if (object$inputs$family=="multinomial") {
      nbeta <- matrix(unlist(nbeta),ncol=length(object$classize))
      colnames(nbeta) <- paste0("class", 1:(length(object$classize))) 
      rownames(nbeta) <- 1:nrow(nbeta)
    }
    return(nbeta)
  }




