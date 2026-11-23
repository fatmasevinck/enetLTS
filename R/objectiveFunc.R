

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

### objective function for multinomial family, from the linear predictors
### `eta` (n x k) instead of a fitted glmnet object; saves a second prediction
### in every C-step
objMultinomEta <- function(eta, y, coef, ind, alpha, lambda){
  pnlty <- lambda * ((1/2) * (1 - alpha) * norm(coef, type="F")^2 +
                       alpha * sum(abs(coef)))
  eta   <- eta[ind, , drop = FALSE]
  prob  <- exp(eta - apply(eta, 1, max))
  prob  <- prob / rowSums(prob)
  mean(-rowSums(y[ind, , drop = FALSE] * log(prob))) + pnlty
}


### Elastic-net penalty P_alpha(b) = sum(1/2 (1 - alpha) b^2 + alpha |b|) of a
### vector or matrix of coefficients (all entries, intercept excluded by the caller)
penaltyEN <- function(b, alpha) sum(0.5 * (1 - alpha) * b^2 + alpha * abs(b))

### Value reported in `$objective`: the penalized criterion of the raw fit on the
### final subset H, on the scale on which the C-steps fit and penalize,
###     (1/h) sum_{i in H} loss_i  +  lambda * P_alpha(s_H * beta),
### where `slopes` are the original-scale coefficients without the intercept and
### `sH` the scale used to standardize x for the fit (so s_H * beta are the
### standardized coefficients). `loss` is the per-observation loss on H.
reportObjective <- function(loss, slopes, sH, alpha, lambda) {
  mean(loss) + lambda * penaltyEN(sH * slopes, alpha)
}
