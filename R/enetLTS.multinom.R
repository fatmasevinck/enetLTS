

enetLTS.multinom <- function(xx, yy, alphas, lambdas, lambdaw, h, hsize, nobs, nvars, intercept, 
                             nsamp, s1, nfold, repl, scal, ncores, nCsteps, tol, seed, plot, typegrouped,
                             type.response){
  
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
    yy <- diag(ny)[as.numeric(yy),]
  } else if (!is.null(k) & k[2]>2){    # if yy is a matrix
    nk <- k[1]
    if (nk!=nobs) stop ( "x and y have different number of rows in call to enetLTS", call.=FALSE )
    k <- as.integer(k[2])
    colnames(yy) <- 1:k
    classnames <- colnames(yy)
    classize <- rep(NA,k)
    for (ii in 1:k){
      classize[ii] <- sum(yy[,ii][yy[,ii]==1])
    }
  }
  
  if (is.null(lambdas)){
    lambdas <- seq(from=0.95,to=0.05,by=-0.05)
  } 
  if (any(lambdas<0)) stop ("lambdas should be non-negative")
  if (any(lambdaw<0)) stop ("lambdaw should be non-negative")

   x <- standardize.x(xx, centerFun=median, scaleFun=mad)
   y <- yy                # not centered for multinomial
   
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
   }
   if (scal){
     xs      <- standardize.x(xx, index=indexbest, centerFun=mean, scaleFun=sd)
     ys      <- yy           # not centered for multinomial
     fit     <- glmnet(xs[indexbest,], ys[indexbest,], family, alpha=alphabest, lambda=lambdabest,
                   standardize=FALSE, intercept=FALSE,type.multinomial=typ)
     a00 <- if (intercept==FALSE) rep(0,k) else drop(t(fit$a0) - 
       t(matrix(do.call(rbind, fit$beta), byrow=T, nrow = k) %*% (attr(xs,"center") / attr(xs,"scale"))))
     raw.coefficients <- drop(matrix(do.call(rbind, fit$beta), byrow=F, ncol = k ) / attr(xs,"scale"))
     # final reweighting:
     md      <- rep(NA,nrow(xs))
     zj      <- drop(predict(fit,newx=xs,type="link"))
     zj.svd  <- svd(scale(zj,TRUE,FALSE))
     svdrank <- sum(zj.svd$d>1e-4)
     if (svdrank<k-1){
        return(list(fit=fit,a00=a00,raw.coefficients=raw.coefficients))
     }
     # compute MDs:
     zj12    <- zj.svd$u[,1:svdrank] # first svdrank PCs
     for (j in 1:k){
        mcdj <- covMcd(zj12[ys[,j]==1,])
        md[y[,j]==1] <- sqrt(mahalanobis(zj12[ys[,j]==1,],mcdj$center,mcdj$cov))
        md[y[,j]==1] <- md[ys[,j]==1]*sqrt(qchisq(0.5,svdrank))/median(md[ys[,j]==1])
     }
     wmd     <-  weight.multinomial(md)
     raw.wt  <- ifelse(wmd != 0, 1, 0)

     xss     <- standardize.x(xx,index=which(raw.wt==1),centerFun=mean, scaleFun=sd)
     yss     <- yy
     
     if (is.null(lambdaw)){
       lambdaw <- cv.glmnet(xss[which(raw.wt==1),], yss[which(raw.wt==1),], family=family, alpha=alphabest,nfolds=5,
                            standardize=FALSE, intercept=FALSE,type.multinomial=typ)$lambda.1se  
     } else if (!is.null(lambdaw) & length(lambdaw)==1){
       lambdaw <- lambdaw
     } else if (!is.null(lambdaw) & length(lambdaw)>1){
       lambdaw <- cv.glmnet(xss[which(raw.wt==1),], yss[which(raw.wt==1),], family=family, lambda=lambdaw, alpha=alphabest,nfolds=5,
                            standardize=FALSE, intercept=FALSE,type.multinomial=typ)$lambda.1se  
     }
     
     fitw    <- glmnet(xss[which(raw.wt==1),], yss[which(raw.wt==1),], family=family, alpha=alphabest, lambda=lambdaw,
                          standardize=FALSE, intercept=FALSE,type.multinomial=typ)
     
     a0      <- if (intercept==FALSE) rep(0,k) else drop(t(fitw$a0) - 
       t(matrix(do.call(rbind, fitw$beta), byrow=T, nrow = k) %*% (attr(xss,"center") / attr(xss,"scale"))))
     coefficients <- drop(matrix(do.call(rbind, fitw$beta), byrow=F, ncol = k ) / attr(xss,"scale"))
  
     # reweighting:
     mdw     <- rep(NA,nrow(xss))
     zjw     <- drop(predict(fitw,newx=xss,type="link"))
     zj.svdw <- svd(scale(zjw,TRUE,FALSE))
     svdrank <- sum(zj.svdw$d>1e-4)
     if (svdrank<k-1){
        return(list(fit=fitw,a0=a0,coefficients=coefficients))
     }
     zj12w   <- zj.svdw$u[,1:svdrank] # first svdrank PCs
     for (j in 1:k){
       mcdjw <- covMcd(zj12w[yss[,j]==1,])
       mdw[yss[,j]==1] <- sqrt(mahalanobis(zj12w[yss[,j]==1,],mcdjw$center,mcdjw$cov))
       mdw[yss[,j]==1] <- mdw[yss[,j]==1]*sqrt(qchisq(0.5,svdrank))/median(mdw[yss[,j]==1])
     }
     wgt <- weight.multinomial(mdw)
     wt  <- ifelse(wgt != 0, 1, 0)
     
   } else { # nonscaled
     fit <- glmnet(xx[indexbest,], yy[indexbest,], family, alpha=alphabest, lambda=lambdabest,
                   standardize=FALSE, intercept=FALSE,type.multinomial=typ)
     a00 <- if (intercept==FALSE) rep(0,k) else drop(t(fit$a0)) 
     raw.coefficients <- drop(matrix(do.call(rbind, fit$beta), byrow=F, ncol = k))
     
     # final reweighting:
     md      <- rep(NA,nrow(xx))
     zj      <- drop(predict(fit,newx=xx,type="link"))
     zj.svd  <- svd(scale(zj,TRUE,FALSE))
     svdrank <- sum(zj.svd$d>1e-4)
     if (svdrank<k-1){
        return(list(fit=fit,a00=a00,raw.coefficients=raw.coefficients))
     }
     zj12    <- zj.svd$u[,1:svdrank] # first svdrank PCs
     for (j in 1:k){
       mcdj  <- covMcd(zj12[yy[,j]==1,])
       md[yy[,j]==1] <- sqrt(mahalanobis(zj12[yy[,j]==1,],mcdj$center,mcdj$cov))
       md[yy[,j]==1] <- md[yy[,j]==1]*sqrt(qchisq(0.5,svdrank))/median(md[yy[,j]==1])
     }
     wmd    <- weight.multinomial(md)
     raw.wt <- ifelse(wmd != 0, 1, 0)
     
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
     coefficients <- drop(matrix(do.call(rbind, fitw$beta), byrow=F, ncol = k))
     
     # reweighting:
     mdw     <- rep(NA,nrow(xx))
     zjw     <- drop(predict(fitw,newx=xx,type="link"))
     zj.svdw <- svd(scale(zjw,TRUE,FALSE))
     svdrank <- sum(zj.svdw$d>1e-4)
     if (svdrank<k-1){
        return(list(fit=fitw,a0=a0,coefficients=coefficients))
     }
     zj12w   <- zj.svdw$u[,1:svdrank] # first svdrank PCs
     for (j in 1:k){
       mcdjw <- covMcd(zj12w[yy[,j]==1,])
       mdw[yy[,j]==1] <- sqrt(mahalanobis(zj12w[yy[,j]==1,],mcdjw$center,mcdjw$cov))
       mdw[yy[,j]==1] <- mdw[yy[,j]==1]*sqrt(qchisq(0.5,svdrank))/median(mdw[yy[,j]==1])
     }
     wgt    <- weight.multinomial(mdw)
     wt     <- ifelse(wgt != 0, 1, 0)
   }
   
   num.nonzerocoef <- sum(coefficients!=0)
   
   # objective based on indexbest
   probj     <- drop(predict(fitw,newx=xx[indexbest,],type="response"))
   deviances <- (-apply(yy[indexbest,]*log(probj),1,sum))
   deviances[is.nan(deviances)] <- 0
   pnlty     <- lambdabest * ((1/2) * (1 - alphabest) * norm(coefficients, type="F")^2 +
                         alphabest * sum(apply(coefficients, 1, function(zz){sum(abs(zz))})))
   objective <- h*mean(deviances) + pnlty
   
   intercept <- isTRUE(intercept)
   if(intercept) xx <- addIntercept(xx)
   
   if (intercept){
     raw.coefficients <- rbind(a00,raw.coefficients)
     colnames(raw.coefficients) <- paste0("class", 1:(length(classize))) 
     rownames(raw.coefficients) <- 1:nrow(raw.coefficients)
     coefficients     <- rbind(a0,coefficients)
     colnames(coefficients) <- paste0("class", 1:(length(classize))) 
     rownames(coefficients) <- 1:nrow(coefficients)
   } else {
     coefficients     <- coefficients
     colnames(coefficients) <- paste0("class", 1:(nrow(classize))) 
     rownames(coefficients) <- 1:length(coefficients)
     raw.coefficients <- raw.coefficients
     colnames(raw.coefficients) <- paste0("class", 1:(length(classize))) 
     rownames(raw.coefficients) <- 1:nrow(raw.coefficients)
   }
   
   # add fitted values 
   raw.yhat  <- xx %*% raw.coefficients
   raw.probs <- exp(raw.yhat)/apply(exp(raw.yhat),1,sum)
   
   yhat      <- xx %*% coefficients
   probs     <- exp(yhat)/apply(exp(yhat),1,sum)
   
   #  deviances 
   deviances <- (-apply(yy*log(probs),1,sum))
   deviances[is.nan(deviances)] <- 0
   raw.deviances <- (-apply(yy*log(raw.probs),1,sum))
   raw.deviances[is.nan(raw.deviances)] <- 0
   
    if (type.response=="link"){
      raw.fitted.values <- raw.yhat
     colnames(raw.fitted.values) <- paste0("class", 1:(length(classize))) 
     rownames(raw.fitted.values) <- 1:nrow(xx)
   } else if (type.response=="class"){
     raw.fitted.values <- apply(raw.probs,1,which.max)
   } else if (type.response=="response"){
     raw.fitted.values <- raw.probs
     colnames(raw.fitted.values) <- paste0("class", 1:(length(classize))) 
     rownames(raw.fitted.values) <- 1:nrow(xx)
   }
   
   if (type.response=="link"){
     fitted.values <- yhat
     colnames(fitted.values) <- paste0("class", 1:(length(classize))) 
     rownames(fitted.values) <- 1:nrow(xx)
   } else if (type.response=="class"){
     fitted.values <- apply(probs,1,which.max)
   } else if (type.response=="response"){
     fitted.values <- probs
     colnames(fitted.values) <- paste0("class", 1:(length(classize))) 
     rownames(fitted.values) <- 1:nrow(xx)
   }
   
   # raw.rmse <- sqrt(mean((yy - raw.yhat)^2))
   # rmse     <- sqrt(mean((yy - yhat)^2))
   raw.mae  <- mean(abs(yy - raw.yhat))
   mae      <- mean(abs(yy - yhat))

   
   inputs <- list(xx=xx, yy=yy, y_factor=y_factor, family=family, alphas=alphas, lambdas=lambdas, lambdaw=lambdaw, hsize=hsize, 
                  intercept=intercept, nsamp=nsamp, s1=s1, nCsteps=nCsteps, nfold=nfold,
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
     lambda            = lambdabest,
     lambdaw           = lambdaw,
     num.nonzerocoef   = num.nonzerocoef,
     n                 = nrow(x),
     p                 = ncol(x),
     h                 = h,
     classnames        = classnames,
     classize          = classize,
     inputs            = inputs
   )
   class(outlist) <- "multinomial"
   return(outlist)
}





