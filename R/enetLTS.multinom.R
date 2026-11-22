
enetLTS.multinom <- function(xx, yy, alphas, lambdas, lambdaw, h, hsize, nobs, nvars, intercept, 
                             nsamp, s1, nfold, repl, scal, ncores, nCsteps, tol, seed, plot, typegrouped,
                             type.response, screen = TRUE, crossfit = FALSE, del = 0.0125){
  
  if (typegrouped) {
    typ <- "grouped"
  } else {
      typ <- "ungrouped"
      }
  
  family <- "multinomial"
  
  k <- dim(yy)
  if (!is.null(k)){
    nc <- as.integer(k[2])
    if (nc==1){
      k <- NULL
    }
  }
  y_factor <- yy
  if (is.null(k)){
    class_sizes <- as.factor(yy)
    classize <- table(class_sizes)
    minclass <- min(classize)
    ny <- as.integer(length(classize))
    classnames <- names(classize)
    if (minclass<=1) stop ( "one multinomial has 1 or 0 observations; not allowed" )
    k <- as.integer(length(classize))
    yy <- diag(ny)[as.integer(class_sizes), , drop = FALSE]
  } else if (!is.null(k) & k[2]>=2){   # if yy is a matrix
    nk <- k[1]
    if (nk!=nobs) stop ( "x and y have different number of rows in call to enetLTS", call.=FALSE )
    k <- as.integer(k[2])
    if (is.null(colnames(yy))) colnames(yy) <- 1:k
    classnames <- colnames(yy)
    classize <- rep(NA,k)
    for (ii in 1:k){
      classize[ii] <- sum(yy[,ii][yy[,ii]==1])
    }
  }
  
  ## class-free screening for leverage points: they are left out of the fit
  xxF <- xx; yyF <- yy; keep <- seq_len(nrow(xx))
  if (isTRUE(screen)) {
    ## the robust covariance (MCD) in the screening uses random subsamples;
    ## fix the RNG so that results are reproducible for a given seed
    if (!is.null(seed)) set.seed(seed)
    flag <- leverageScreen(xx)
    if (any(flag)) {
      kp <- which(!flag)
      if (all(colSums(yy[kp, , drop = FALSE]) >= 5)) {
        keep <- kp; xx <- xx[keep, , drop = FALSE]; yy <- yy[keep, , drop = FALSE]
      } else flag[] <- FALSE
    }
  } else flag <- rep(FALSE, nrow(xx))
  nrF <- nrow(xxF)

  if (any(lambdas<0)) stop ("lambdas should be non-negative")
  if (any(lambdaw<0)) stop ("lambdaw should be non-negative")

   x <- standardize.x(xx, centerFun=median, scaleFun=mad)
   y <- yy                # not centered for multinomial

   if (is.null(lambdas)) lambdas <- multinomLambdas(x, y)

   class_sizes <- apply(y,2,sum)
   hk          <- floor(class_sizes*hsize+1)
   h           <- sum(hk)

   indexall <- warmCsteps.mod(x, y, h, hsize, alphas, lambdas, nsamp, s1, 
                              scal, ncores, nCsteps, tol, seed, family)
   
   if ((length(alphas)==1) & (length(lambdas)==1)){
      if (plot==TRUE) warning("There is no meaning to see plot for a single 
                              combination of alpha and lambda")
      indexbest  <- drop(indexall)
      alphabest  <- alphas
      lambdabest <- lambdas
   } else {
      CVresults  <- cv.multinom.enetLTS(indexall, x, y, alphas, lambdas, 
                                       nfold, repl, ncores, plot)
      indexbest  <- CVresults$indexbest
      alphabest  <- CVresults$alphaopt
      lambdabest <- CVresults$lambdaopt
      evalCritCV <- CVresults$evalCrit
      if (length(lambdas) > 1 && isTRUE(all.equal(lambdabest, min(lambdas))))
        warning("the optimal lambda is the smallest value of the grid; ",
                "consider a grid that extends to smaller values", call. = FALSE)
   }
   if (isTRUE(crossfit)) {
     if (!scal) stop("'crossfit' requires scal=TRUE")
     cf <- multinomCrossfitTrim(x, max.col(y, ties.method="first"), ncol(y), indexbest, hsize,
                                alphas, nfold, del)
     indexbest  <- cf$H
     alphabest  <- cf$alpha
     lambdabest <- cf$lambda
   }
   alpha0 <- alphabest; lambda0 <- lambdabest      # (alpha, lambda) of the raw fit, for $objective
   if (scal){
     xs      <- standardize.x(xx, index=indexbest, centerFun=mean, scaleFun=sd)
     ys      <- yy           # not centered for multinomial
     fit     <- glmnet(xs[indexbest,], ys[indexbest,], family, alpha=alphabest, lambda=lambdabest,
                   standardize=FALSE, intercept=FALSE,type.multinomial=typ)
     a00 <- if (intercept==FALSE) rep(0,k) else drop(t(fit$a0) - 
       t(matrix(do.call(rbind, fit$beta), byrow=TRUE, nrow = k) %*% (attr(xs,"center") / attr(xs,"scale"))))
     raw.coefficients <- drop(matrix(do.call(rbind, fit$beta), byrow=FALSE, ncol = k ) / attr(xs,"scale"))
     sH <- attr(xs,"scale")
     # final reweighting:
     zj      <- drop(predict(fit,newx=xs,type="link"))
     raw.wt  <- multinomReweight(zj, ys)
     if (isTRUE(crossfit)) raw.wt <- as.integer(seq_len(nrow(xx)) %in% cf$W)

     xss     <- standardize.x(xx,index=which(raw.wt==1),centerFun=mean, scaleFun=sd)
     yss     <- yy
     
     if (isTRUE(crossfit)) {
       wi  <- which(raw.wt==1)
       sel <- multinomCrossfitFinal(xss[wi,,drop=FALSE], yss[wi,,drop=FALSE], nfold, typ,
                                   alphas = finalAlphaGrid(alphas, c(0, 0.1, 0.3, 0.5, 0.7, 1)))
       alphabest <- sel$alpha
       if (is.null(lambdaw) || length(lambdaw) != 1) lambdaw <- sel$lambda
     } else if (is.null(lambdaw)){
       lambdaw <- cv.glmnet(xss[which(raw.wt==1),], yss[which(raw.wt==1),], family=family, alpha=alphabest,nfolds=5,
                            standardize=FALSE, intercept=FALSE,type.multinomial=typ)$lambda.1se  
     } else if (!is.null(lambdaw) & length(lambdaw)==1){
       lambdaw <- lambdaw
     } else if (!is.null(lambdaw) & length(lambdaw)>1){
       lambdaw <- cv.glmnet(xss[which(raw.wt==1),], yss[which(raw.wt==1),], family=family, lambda=lambdaw, alpha=alphabest,nfolds=5,
                            standardize=FALSE, intercept=FALSE,type.multinomial=typ)$lambda.1se  
     }
     
     fitw    <- glmnet(xss[which(raw.wt==1),], yss[which(raw.wt==1),], family=family, alpha=alphabest, lambda=lambdaw,
                          standardize=FALSE, intercept=isTRUE(crossfit), type.multinomial=typ)
     
     a0      <- if (intercept==FALSE) rep(0,k) else drop(t(fitw$a0) - 
       t(matrix(do.call(rbind, fitw$beta), byrow=TRUE, nrow = k) %*% (attr(xss,"center") / attr(xss,"scale"))))
     coefficients <- drop(matrix(do.call(rbind, fitw$beta), byrow=FALSE, ncol = k ) / attr(xss,"scale"))
  
     # reweighting:
     zjw     <- drop(predict(fitw,newx=xss,type="link"))
     wt      <- if (isTRUE(crossfit)) raw.wt else multinomReweight(zjw, yss)
     
   } else { # nonscaled
     fit <- glmnet(xx[indexbest,], yy[indexbest,], family, alpha=alphabest, lambda=lambdabest,
                   standardize=FALSE, intercept=FALSE,type.multinomial=typ)
     a00 <- if (intercept==FALSE) rep(0,k) else drop(t(fit$a0)) 
     raw.coefficients <- drop(matrix(do.call(rbind, fit$beta), byrow=FALSE, ncol = k))
     sH <- rep(1, ncol(xx))          # this branch fits on the unscaled xx
     
     # final reweighting:
     zj      <- drop(predict(fit,newx=xx,type="link"))
     raw.wt  <- multinomReweight(zj, yy)
     
     if (is.null(lambdaw)){
       lambdaw <- cv.glmnet(xx[which(raw.wt==1),], yy[which(raw.wt==1),], family=family, alpha=alphabest,nfolds=5,
                            standardize=FALSE, intercept=FALSE,type.multinomial=typ)$lambda.1se  
     } else if (!is.null(lambdaw) & length(lambdaw)==1){
       lambdaw <- lambdaw
     } else if (!is.null(lambdaw) & length(lambdaw)>1){
       lambdaw <- cv.glmnet(xx[which(raw.wt==1),], yy[which(raw.wt==1),], family=family, alpha=alphabest,lambda=lambdaw,nfolds=5,
                            standardize=FALSE, intercept=FALSE,type.multinomial=typ)$lambda.1se  
     }
     
     fitw    <- glmnet(xx[which(raw.wt==1),], yy[which(raw.wt==1),], family=family, alpha=alphabest, lambda=lambdaw,
                    standardize=FALSE, intercept=FALSE, type.multinomial=typ)
     
     a0      <- if (intercept==FALSE) rep(0,k) else drop(t(fitw$a0))
     coefficients <- drop(matrix(do.call(rbind, fitw$beta), byrow=FALSE, ncol = k))
     
     # reweighting:
     zjw     <- drop(predict(fitw,newx=xx,type="link"))
     wt      <- multinomReweight(zjw, yy)
   }
   
   ## back to all observations: screened observations get weight zero
   if (length(keep) < nrF) {
     rw <- integer(nrF); rw[keep] <- raw.wt; raw.wt <- rw
     ww <- integer(nrF); ww[keep] <- wt;     wt     <- ww
     indexbest <- keep[indexbest]
     xx <- xxF; yy <- yyF
   }

   num.nonzerocoef <- sum(coefficients!=0)
   
   intercept <- isTRUE(intercept)
   if(intercept) xx <- addIntercept(xx)
   
   if (intercept){
     raw.coefficients <- rbind(a00,raw.coefficients)
     coefficients     <- rbind(a0,coefficients)
   } else {
     coefficients     <- coefficients
     raw.coefficients <- raw.coefficients
   }
   
   # add fitted values 
   raw.yhat  <- xx %*% raw.coefficients
   raw.probs <- multinomProbs(raw.yhat)
   
   yhat      <- xx %*% coefficients
   probs     <- multinomProbs(yhat)
   
   #  deviances 
   deviances     <- multinomDeviance(yhat, yy)
   raw.deviances <- multinomDeviance(raw.yhat, yy)
   
   # $objective: penalized criterion of the raw fit on the subset 'indexbest' (see reportObjective());
   # the intercept row (if any) is not penalized
   objective <- reportObjective(raw.deviances[indexbest],
                                if (intercept) raw.coefficients[-1,,drop=FALSE] else raw.coefficients,
                                sH, alpha0, lambda0)
   
   raw.fitted.values <- if (type.response=="link"){
     raw.yhat
   } else if (type.response=="class"){
     classnames[apply(raw.probs,1,which.max)]
   } else if (type.response=="response"){
     raw.probs
   }
   
   fitted.values     <- if (type.response=="link"){
     yhat
   } else if (type.response=="class"){
     classnames[apply(probs,1,which.max)]
   } else if (type.response=="response"){
     probs
   }
   
   ## mean absolute difference between the 0/1 response and the fitted probabilities
   raw.mae  <- mean(abs(yy - raw.probs))
   mae      <- mean(abs(yy - probs))

   
   inputs <- list(xx=xx, yy=yy, y_factor=y_factor, family=family, alphas=alphas, lambdas=lambdas, lambdaw=lambdaw, hsize=hsize, 
                  intercept=intercept, screen=screen, crossfit=crossfit, nsamp=nsamp, s1=s1, nCsteps=nCsteps, nfold=nfold,
                  repl=repl, ncores=ncores, tol=tol, scal=scal, seed=seed, plot=plot)
   outlist <- list(
     objective         = objective,
     raw.mae           = raw.mae,
     mae               = mae,
     best              = sort(indexbest),
     raw.wt            = raw.wt,
     wt                = wt,
     raw.coefficients  = raw.coefficients,
     coefficients      = coefficients,
     raw.fitted.values = raw.fitted.values,
     fitted.values     = fitted.values,
     raw.residuals     = drop(raw.deviances),
     residuals         = drop(deviances),
     alpha             = alphabest,
     alpha.raw         = alpha0,
     lambda            = lambdabest,
     lambdaw           = lambdaw,
     num.nonzerocoef   = num.nonzerocoef,
     n                 = nrF,
     p                 = ncol(x),
     h                 = h,
     classnames        = classnames,
     classize          = classize,
     screened          = which(flag),
     inputs            = inputs
   )
   class(outlist) <- "multinomial"
   return(outlist)
}




