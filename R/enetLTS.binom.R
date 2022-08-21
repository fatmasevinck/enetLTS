


enetLTS.binom <- function(xx, yy, alphas, lambdas, lambdaw, h, hsize, nobs, nvars, intercept, nsamp,
                          s1, nfold, repl, scal, iniscal, ncores, nCsteps, tol, seed, del, plot,
                          type.response){

   family <- "binomial"

   k <- dim(yy)
   if (!is.null(k)){
      nc <- as.integer(k[2])
      if (nc>2) stop ("More than two classes; not available for binomial family; use multinomial family",call.=FALSE)
      if (nc==1){ k <- NULL }
      if (nc==2){
        yy <- drop(yy[,2])
        k <- NULL
        }
   }
   if (is.null(k)){
      class_sizes <- as.factor(yy)
      classize <- table(class_sizes)
      minclass <- min(classize)
      if (minclass<8) warning ("one binomial class has fewer than 8  observations; dangerous ground")
      if (minclass<=1) stop ("one binomial class has 1 or 0 observations; not allowed")
      classnames <- names(classize)
      }

   if (length(yy)!=nobs) stop ("x and y have different number of rows in call to glmnet",call.=FALSE)

   if (is.null(lambdas)){
    l0 <- lambda00(xx, yy, normalize=scal, intercept=intercept)
    lambdas <-  seq(l0, 0, by=-0.025*l0)
   }
   if (any(lambdas<0)) stop ("lambdas should be non-negative")
   if (any(lambdaw<0)) stop ("lambdaw should be non-negative")

   if (isTRUE(iniscal)){
     x <- standardize.x(xx, centerFun=median, scaleFun=mad)
     y <- yy   # not centered for binomial
   } else {  # not initial robust scale!!!!! Careful!!!
     x <- xx
     y <- yy
   }

   indexall <- warmCsteps.mod(x, y, h, hsize, alphas, lambdas, nsamp, s1,
                              scal, ncores, nCsteps, tol, seed, family)

   if ((length(alphas)==1) & (length(lambdas)==1)){
     if (plot==TRUE) warning("There is no meaning to see plot for a single
                              combination of alpha and lambda")
     indexbest   <- drop(indexall)
     alphabest   <- alphas
     lambdabest  <- lambdas
   } else {
     CVresults   <- cv.binomial.enetLTS(indexall, x, y, alphas, lambdas, nfold, repl, ncores, plot)
     indexbest   <- CVresults$indexbest
     alphabest   <- CVresults$alphaopt
     lambdabest  <- CVresults$lambdaopt
     evalCritCV  <- CVresults$evalCrit
   }
   if (scal){
     xs <- standardize.x(xx, index=indexbest, centerFun=mean, scaleFun=sd)
     ys <- yy # # not centered for binomail

     fit <- glmnet(xs[indexbest,], ys[indexbest], family, alpha=alphabest, lambda=lambdabest,
                   standardize=FALSE, intercept=FALSE)

     a00 <- if (intercept==FALSE) 0 else drop(fit$a0 -
                                                as.vector(as.matrix(fit$beta)) %*% (attr(xs,"center") / attr(xs,"scale")))
     raw.coefficients <- drop(as.matrix(fit$beta) / attr(xs,"scale"))
     # final reweighting:
     # raw.residuals <- -(ys * xs %*% as.matrix(fit$beta)) + log(1 + exp(xs %*% as.matrix(fit$beta)))

     raw.wt <- weight.binomial(xx, yy, c(a00, raw.coefficients), intercept=intercept, del)
     xss <- standardize.x(xx, index=which(raw.wt==1), centerFun=mean, scaleFun=sd)
     yss <- yy

     if (is.null(lambdaw)){
       lambdaw <- cv.glmnet(xss[which(raw.wt==1),],yss[which(raw.wt==1)],family=family,nfolds=5,
                            alpha=alphabest,standardize=FALSE,intercept=FALSE,type.measure="mse")$lambda.min
     } else if (!is.null(lambdaw) & length(lambdaw)==1){
       lambdaw <- lambdaw
     } else if (!is.null(lambdaw) & length(lambdaw)>1){
       lambdaw <- cv.glmnet(xss[which(raw.wt==1),],yss[which(raw.wt==1)],family=family,lambda=lambdaw,nfolds=5,
                            alpha=alphabest,standardize=FALSE,intercept=FALSE,type.measure="mse")$lambda.min
     }

     fitw <- glmnet(xss[which(raw.wt==1),],yss[which(raw.wt==1)],family,alpha=alphabest,lambda=lambdaw,
                    standardize=FALSE,intercept=FALSE)  ## now take raw.wt instead of index

     a0 <- if (intercept==FALSE) 0 else drop(fitw$a0 -
                                               as.vector(as.matrix(fitw$beta)) %*% (attr(xss,"center") / attr(xss,"scale")))
     coefficients <- drop(as.matrix(fitw$beta)/attr(xss,"scale"))

     wgt <- weight.binomial(xx,yy,c(a0,coefficients),intercept,del)

  #   reweighted.residuals  <- -(yy * cbind(1,xx) %*% c(a0,coefficients)) + log(1+exp(cbind(1,xx) %*% c(a0,coefficients)))

   } else { # nonscaled
       fit <- glmnet(xx[indexbest,], yy[indexbest], family, alpha=alphabest, lambda=lambdabest,
                     standardize=FALSE, intercept=FALSE)

       a00 <- if (intercept==FALSE) 0 else drop(fit$a0)
       raw.coefficients <- drop(as.matrix(fit$beta))
       # final reweighting:
       # raw.residuals <- -(yy * xx %*% as.matrix(fit$beta)) + log(1+exp(xx %*% as.matrix(fit$beta)))
       raw.wt <- weight.binomial(xx,yy,c(a00,raw.coefficients),intercept,del)
       if (is.null(lambdaw)){
         lambdaw <- cv.glmnet(xx[which(raw.wt==1),],yy[which(raw.wt==1)],family=family,nfolds=5,
                              alpha=alphabest,standardize=FALSE,intercept=FALSE,type.measure="mse")$lambda.min
       } else if (!is.null(lambdaw) & length(lambdaw)==1){
         lambdaw <- lambdaw
       } else if (!is.null(lambdaw) & length(lambdaw)>1){
         lambdaw <- cv.glmnet(xx[which(raw.wt==1),],yy[which(raw.wt==1)],family=family,lambda=lambdaw,nfolds=5,
                              alpha=alphabest,standardize=FALSE,intercept=FALSE,type.measure="mse")$lambda.min
       }
       fitw <- glmnet(xx[which(raw.wt==1),],yy[which(raw.wt==1)],family,alpha=alphabest,lambda=lambdaw,
                      standardize=FALSE,intercept=FALSE)  ## now we take raw.wt instead of index

       a0 <- if (intercept==FALSE) 0 else drop(fitw$a0)
       coefficients <- drop(as.matrix(fitw$beta))

       wgt <- weight.binomial(xx, yy, c(a0,coefficients), intercept, del)
    }

   num.nonzerocoef <- sum(coefficients!=0)

   intercept <- isTRUE(intercept)
   if(intercept) xx <- addIntercept(xx)

   if (intercept){
      coefficients            <- c(a0,coefficients)
      names(coefficients)     <- 1:length(coefficients)
      raw.coefficients        <- c(a00,raw.coefficients)
      names(raw.coefficients) <- 1:length(raw.coefficients)
   } else {
      coefficients            <- coefficients
      names(coefficients)     <- 1:length(coefficients)
      raw.coefficients        <- raw.coefficients
      names(raw.coefficients) <- 1:length(raw.coefficients)
   }

   raw.yhat <- xx %*% raw.coefficients
   raw.probs <- 1/(1+exp(-raw.yhat))

   yhat     <- xx %*% coefficients
   probs <- 1/(1+exp(-yhat))

   # deviances
   raw.deviances  <- -(yy * xx %*% raw.coefficients) + log(1+exp(xx %*% raw.coefficients))
   deviances      <- -(yy * xx %*% coefficients) + log(1+exp(xx %*% coefficients))

   raw.fitted.values <-  if (type.response=="link"){
     raw.yhat
   } else if (type.response=="class"){
     ifelse(raw.yhat <= 0.5,0,1)
   } else if (type.response=="response"){
     raw.probs
   }

   fitted.values     <-  if (type.response=="link"){
     yhat
   } else if (type.response=="class"){
     ifelse(yhat <= 0.5,0,1)
   } else if (type.response=="response"){
     probs
   }

   # raw.rmse <- sqrt(mean((yy - raw.yhat)^2))
   # rmse     <- sqrt(mean((yy - yhat)^2))
   raw.mae  <- mean(abs(yy - raw.yhat))
   mae      <- mean(abs(yy - yhat))


   objective <- h * (mean((-yy[indexbest] * (xx[indexbest,] %*% coefficients)) +
                            log(1+exp(xx[indexbest,] %*% coefficients))) +
                       lambdabest * sum(1/2 * (1-alphabest) * coefficients^2 +
                                          alphabest*abs(coefficients)))

   inputs <- list(xx=xx, yy=yy, family=family, alphas=alphas, lambdas=lambdas, lambdaw=lambdaw,
                  hsize=hsize, intercept=intercept, nsamp=nsamp, s1=s1, nCsteps=nCsteps, nfold=nfold,
                  repl=repl, ncores=ncores, del=del, tol=tol, scal=scal, seed=seed, plot=plot,
                  type.response=type.response)
   outlist <- list(
      objective         = objective,
      raw.mae           = raw.mae,
      mae               = mae,
      best              = sort(indexbest),
      raw.wt            = raw.wt,
      wt                = wgt,
      raw.coefficients  = raw.coefficients,
      coefficients      = coefficients,
      raw.fitted.values = drop(raw.fitted.values),
      fitted.values     = drop(fitted.values),
      raw.residuals     = drop(raw.deviances),
      residuals         = drop(deviances),
      alpha             = alphabest,
      lambda            = lambdabest,
      lambdaw           = lambdaw,
      num.nonzerocoef   = num.nonzerocoef,
      n                 = nrow(x),
      p                 = ncol(x),
      h                 = h,
      classnames        = classnames,
      classize          = classize,
      inputs            = inputs
      # indexall = indexall
   )
   class(outlist) <- "binomial"
   return(outlist)
}
