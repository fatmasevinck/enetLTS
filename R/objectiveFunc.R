

### objective function for multinomial family
objMultinom <- function(x, y, fit, coef, ind, alpha, lambda){
  pnlty <- lambda * ((1/2) * (1 - alpha) * norm(coef, type="F")^2 +
                       alpha * sum(apply(coef, 1, function(zz){sum(abs(zz))})))
  if (is.null(ind)){
    probj <- drop(predict(fit,newx=x,type="response"))
    deviances <- (-apply(y*log(probj),1,sum))
    obj <- mean(deviances) + pnlty
  } else {
    probj <- drop(predict(fit,newx=x[ind,],type="response"))
    deviances <- (-apply(y[ind,]*log(probj),1,sum))
    obj <- mean(deviances) + pnlty
  }
  return(obj)
}



### objective function for binomial family
objBinom <- function(x, y, coef, ind, alpha, lambda){
   if (is.null(ind)){
      obj <- mean(phiBY3(x %*% coef,ifelse(y==0,1,0),c3=0.5))
   } else {
      obj <- mean(phiBY3(x[ind,] %*% coef,ifelse(y[ind]==0,1,0),c3=0.5))
   }
   return(obj)
}


### objective function for gaussian family
objGaus <- function(x, y, coef, ind, alpha, lambda){
   if (is.null(ind)){
      obj <- (1/2) * mean((y - x %*% coef)^2) + 
         lambda * sum(1/2 * (1-alpha) * coef^2 + alpha * abs(coef))
   } else {
      obj <- (1/2) * mean((y[ind] - x[ind,] %*% coef)^2) +
         lambda * sum(1/2 * (1-alpha) * coef^2 + alpha * abs(coef))
   }
   return(obj)
}



