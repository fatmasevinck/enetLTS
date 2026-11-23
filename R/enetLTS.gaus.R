
enetLTS.gaus <- function(xx, yy, alphas, lambdas, lambdaw, h, hsize, nobs, nvars, intercept, nsamp, 
                         s1, nfold, repl, scal, ncores, nCsteps, tol, seed, del, plot,
                         type.response, crossfit = FALSE){
   
   family <- "gaussian"

   if (is.null(lambdas)){
     l0 <- robustHD::lambda0(xx, yy, normalize=scal, intercept=intercept)
     lambdas <- seq(l0, 0, by=-0.025*l0)
   }
   
   if (any(lambdas<0)) stop ("lambdas should be non-negative")
   if (any(lambdaw<0)) stop ("lambdaw should be non-negative")
   
   x  <- standardize.x(xx, centerFun=median, scaleFun=mad)
   y  <- center.y(yy, centerFun=median)

   indexall  <- warmCsteps.mod(x, y, h, hsize, alphas, lambdas, nsamp, s1, 
                              scal, ncores, nCsteps, tol, seed, family)
   
   if ((length(alphas)==1) & (length(lambdas)==1)){
      if (plot==TRUE) warning("There is no meaning to see plot for a single 
                              combination of alpha and lambda")
      indexbest   <- drop(indexall)
      alphabest   <- alphas
      lambdabest  <- lambdas
   } else {
      CVresults   <- cv.gaussian.enetLTS(indexall, x, y, alphas, lambdas, 
                                   nfold, repl, ncores, plot)
      indexbest   <- CVresults$indexbest
      alphabest   <- CVresults$alphaopt
      lambdabest  <- CVresults$lambdaopt
      evalCritCV  <- CVresults$evalCrit
   }
   if (isTRUE(crossfit)) {
     if (!scal) stop("'crossfit' requires scal=TRUE")
     cf <- crossfitTrim(x, y, indexbest, h, alphas, nfold, del)
     indexbest  <- cf$H
     alphabest  <- cf$alpha
     lambdabest <- cf$lambda
   }
   alpha0 <- alphabest; lambda0 <- lambdabest      # (alpha, lambda) of the raw fit, for $objective
   if (scal){
     xs     <- standardize.x(xx, index=indexbest, centerFun=mean, scaleFun=sd)
     ys     <- center.y(yy, index=indexbest, centerFun=mean)
     
     fit    <- glmnet(xs[indexbest,], ys[indexbest], alpha=alphabest, lambda=lambdabest,
                    standardize=FALSE, intercept=FALSE)
     a00                   <- if (intercept==FALSE) 0 else drop(attr(ys,"center") + fit$a0 - 
                                                as.vector(as.matrix(fit$beta)) %*% (attr(xs,"center") / attr(xs,"scale")))
     raw.coefficients      <- drop(as.matrix(fit$beta) / attr(xs,"scale"))
     sH                    <- attr(xs,"scale")
     
     raw.residuals         <- yy - cbind(1,xx) %*% c(a00,raw.coefficients)
     # final reweighting:
     raw.wt <- weight.gaussian(raw.residuals, indexbest, del)$we
     if (isTRUE(crossfit)) raw.wt <- as.integer(seq_len(nrow(xx)) %in% cf$W)
     
     xss <- standardize.x(xx, index=which(raw.wt==1), centerFun=mean, scaleFun=sd)
     yss <- center.y(yy, index=which(raw.wt==1), centerFun=mean)
     
     if (isTRUE(crossfit)) {
       wi  <- which(raw.wt==1)
       sel <- crossfitFinal(xss[wi,,drop=FALSE], yss[wi], nfold,
                           alphas = finalAlphaGrid(alphas, seq(0, 1, length.out = 21)))
       alphabest <- sel$alpha
       if (is.null(lambdaw) || length(lambdaw) != 1) lambdaw <- sel$lambda
     } else if (is.null(lambdaw)){
       lambdaw <- cv.glmnet(xss[which(raw.wt==1),], yss[which(raw.wt==1)], family=family, nfolds=nfold,
                            alpha=alphabest, standardize=FALSE, intercept=FALSE, type.measure="mse")$lambda.min
     } else if (!is.null(lambdaw) & length(lambdaw)==1){
       lambdaw <- lambdaw
     } else if (!is.null(lambdaw) & length(lambdaw)>1){
       lambdaw <- cv.glmnet(xss[which(raw.wt==1),], yss[which(raw.wt==1)], family=family, lambda=lambdaw, nfolds=nfold,
                            alpha=alphabest, standardize=FALSE, intercept=FALSE, type.measure="mse")$lambda.min
     }
     
     fitw <- glmnet(xss[which(raw.wt==1),], yss[which(raw.wt==1)], alpha=alphabest, lambda=lambdaw,
                    standardize=FALSE, intercept=FALSE)
     a0                   <- if (intercept==FALSE) 0 else drop(attr(yss,"center") + fitw$a0 - 
                                               as.vector(as.matrix(fitw$beta)) %*% (attr(xss,"center") / attr(xss,"scale")))
     coefficients         <- drop(as.matrix(fitw$beta)/attr(xss,"scale"))
     
     reweighted.residuals <- yy - cbind(1,xx) %*% c(a0,coefficients)
     wgt <- weight.gaussian(reweighted.residuals, raw.wt==1, del)$we
     
   } else {                # scal = FALSE
     fit <- glmnet(x[indexbest,], y[indexbest], alpha=alphabest, lambda=lambdabest,
                   standardize=FALSE, intercept=FALSE)
     a00                  <- if (intercept==FALSE) 0 else drop(attr(y,"center") + fit$a0 - 
                                                as.vector(as.matrix(fit$beta)) %*% (attr(x,"center") / attr(x,"scale")))
     raw.coefficients     <- drop(as.matrix(fit$beta)/attr(x,"scale"))
     sH                   <- attr(x,"scale")
     raw.residuals        <- yy - cbind(1,xx) %*% c(a00,raw.coefficients)
     # final reweighting:
     raw.wt               <- weight.gaussian(raw.residuals, indexbest, del)$we 
     if (is.null(lambdaw)){
       lambdaw <- cv.glmnet(x[which(raw.wt==1),], y[which(raw.wt==1)], nfolds=nfold, alpha=alphabest,
                            standardize=FALSE, intercept=FALSE, type.measure="mse")$lambda.min
     } else if (!is.null(lambdaw) & length(lambdaw)==1){
       lambdaw <- lambdaw
     } else if (!is.null(lambdaw) & length(lambdaw)>1){
       lambdaw <- cv.glmnet(x[which(raw.wt==1),], y[which(raw.wt==1)], lambda=lambdaw, nfolds=nfold, alpha=alphabest,
                            standardize=FALSE, intercept=FALSE, type.measure="mse")$lambda.min
     }
     fitw <- glmnet(x[which(raw.wt==1),], y[which(raw.wt==1)], alpha=alphabest, lambda=lambdaw,
                    standardize=FALSE, intercept=FALSE)
     
     a0                   <- if (intercept==FALSE) 0 else drop(attr(y,"center") + 
                                               fitw$a0-as.vector(as.matrix(fitw$beta))%*%(attr(x,"center")/attr(x,"scale")))
     coefficients         <- drop(as.matrix(fitw$beta)/attr(x,"scale"))
     reweighted.residuals <-  yy - cbind(1,xx) %*% c(a0,coefficients)
     wgt                  <- weight.gaussian(reweighted.residuals, raw.wt==1, del)$we
     ## back-transformed to the original scale
   }
   num.nonzerocoef  <- sum(coefficients!=0)
   
   intercept        <- isTRUE(intercept)
   if(intercept) xx <- addIntercept(xx)
   
   if (intercept){
      coefficients     <- c(a0,coefficients)
      raw.coefficients <- c(a00,raw.coefficients)
   } else {
      coefficients     <- coefficients
      raw.coefficients <- raw.coefficients
   }
   
   raw.yhat <- xx %*% raw.coefficients
   yhat     <- xx %*% coefficients
   
   raw.fitted.values <-  if (type.response=="link" | type.response=="response"){
     raw.yhat
   } 
   
   fitted.values     <-  if (type.response=="link" | type.response=="response"){
     yhat
   }
   
   raw.residuals <- yy - raw.yhat
   residuals     <- yy - yhat
   
   raw.rmse <- sqrt(mean((yy - raw.yhat)^2))
   raw.mae  <- mean(abs(yy - raw.yhat))
   rmse     <- sqrt(mean((yy - yhat)^2))
   mae      <- mean(abs(yy - yhat))
   
   # $objective: penalized criterion of the raw fit on the subset 'indexbest' (see reportObjective())
   objective <- reportObjective(0.5 * (yy[indexbest] - xx[indexbest,,drop=FALSE] %*% raw.coefficients)^2,
                                if (intercept) raw.coefficients[-1] else raw.coefficients, sH, alpha0, lambda0)
   
   
   inputs <- list(xx=xx, yy=yy, family=family, crossfit=crossfit, alphas=alphas, lambdas=lambdas, lambdaw=lambdaw, hsize=hsize, 
                  intercept=intercept, nsamp=nsamp, s1=s1, nCsteps=nCsteps, nfold=nfold,
                  repl=repl, ncores=ncores, del=del, tol=tol, scal=scal, seed=seed, plot=plot)
   outlist <- list(
      objective         = objective,
      raw.rmse          = raw.rmse,
      rmse              = rmse,
      raw.mae           = raw.mae,
      mae               = mae,
      best              = sort(indexbest),
      raw.wt            = raw.wt,
      wt                = wgt,
      raw.coefficients  = raw.coefficients,
      coefficients      = coefficients,
      raw.fitted.values = raw.fitted.values,
      fitted.values     = fitted.values,
      raw.residuals     = drop(raw.residuals),
      residuals         = drop(residuals),
      alpha             = alphabest,
      alpha.raw         = alpha0,
      lambda            = lambdabest,
      lambdaw           = lambdaw,
      num.nonzerocoef   = num.nonzerocoef,
      n                 = nrow(x),
      p                 = ncol(x),
      h                 = h,
      inputs            = inputs
   )
   class(outlist) <- "gaussian"
   return(outlist)
}