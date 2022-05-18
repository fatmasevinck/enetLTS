## internal functions
# require(glmnet)
# source("utilities.R")
# source("objectiveFunc.R")

CStep.multinom <- 
  function(x, y, indx, h, hsize, alpha, lambda, scal)
  {
    k <- ncol(y)
    class_sizes <- apply(y,2,sum)
    hk <- floor(class_sizes*hsize+1)
    if (scal){
      xs <- standardize.x(x, index=indx, centerFun=mean, scaleFun=sd)
      ys <- y        # not centered response variable for multinomial case
      
      fit <- glmnet(xs[indx,], ys[indx,], family="multinomial", alpha=alpha, lambda=lambda, 
                    standardize=FALSE, intercept=FALSE)
      ### standardize=FALSE, intercept=FALSE,type.multinomial="grouped")
      beta <- matrix(do.call(rbind,lapply(fit$beta,matrix)),ncol=k) # to convert a matrix 
      md <- rep(NA,nrow(xs))
      zj <- drop(predict(fit,newx=xs,type="link"))
      zj.svd <- svd(scale(zj,TRUE,FALSE))
      svdrank <- sum(zj.svd$d>1e-4)
      # buraya svdrank hakkinda bir uyari eklenebilir. paket icin uyari olarak. 
      # cunku bunun full rank saglandigindan cok da emin degil.  
      if (svdrank<k-1){
        return(list(object=-Inf,index=indx,md=md,beta=beta))
      }
      
      # compute MDs:
      zj12 <- zj.svd$u[,1:svdrank] # first svdrank PCs
      for (j in 1:k){
        mcdj <- covMcd(zj12[ys[,j]==1,])
        md[ys[,j]==1] <- sqrt(mahalanobis(zj12[ys[,j]==1,],mcdj$center,mcdj$cov))
        md[ys[,j]==1] <- md[ys[,j]==1]*sqrt(qchisq(0.5,svdrank))/median(md[ys[,j]==1])
      }
      # P wmd <- weightmd(md)
      md.sort <- sort(md, decreasing=FALSE, index.return=TRUE)
      indxnew <- NULL
      for (l in 1:k){
        indxnew <- c(indxnew,md.sort$ix[y[,l][md.sort$ix]==1][1:hk[l]])  
      }
      obj <- objMultinom(xs, ys, fit, beta, indxnew, alpha, lambda)
      
    } else {
      
      fit <- glmnet(x[indx,], y[indx,], family="multinomial", alpha=alpha, lambda=lambda, 
                    standardize=FALSE, intercept=FALSE)
      ###                  standardize=FALSE, intercept=FALSE,type.multinomial="grouped")
      beta <- matrix(do.call(rbind,lapply(fit$beta,matrix)),ncol=k) # to convert a matrix 
      md <- rep(NA,nrow(x))
      zj <- drop(predict(fit,newx=x,type="link"))
      zj.svd <- svd(scale(zj,TRUE,FALSE))
      svdrank <- sum(zj.svd$d>1e-4)
      if (svdrank<k-1){
        return(list(object=-Inf,index=indx,md=md,beta=beta))
      }
      # compute MDs:
      zj12 <- zj.svd$u[,1:svdrank] # first svdrank PCs
      for (j in 1:k){
        mcdj <- covMcd(zj12[y[,j]==1,])
        md[y[,j]==1] <- sqrt(mahalanobis(zj12[y[,j]==1,],mcdj$center,mcdj$cov))
        md[y[,j]==1] <- md[y[,j]==1]*sqrt(qchisq(0.5,svdrank))/median(md[y[,j]==1])
      }
      # wmd <- weightmd(md)
      md.sort <- sort(md, decreasing=FALSE, index.return=TRUE)
      class_sizes <- apply(y,2,sum)
      hk <- floor(class_sizes*hsize+1)
      indxnew <- NULL
      for (l in 1:k){
        indxnew <- c(indxnew,md.sort$ix[y[,l][md.sort$ix]==1][1:hk[l]])  
      }
      obj <- objMultinom(x, y, fit, beta, indxnew, alpha, lambda) 
    }
    return(list(object=obj, index=indxnew, md=md, beta=beta))
  }


#####################################################################################################################

CStep.binom <- 
   function(x, y, indx, h, hsize, alpha, lambda, scal)
   {
     
      if (scal){
         xs <- standardize.x(x, index=indx, centerFun=mean, scaleFun=sd)
         ys <- y        # not centered response variable for binomial case
         
         fit <- glmnet(xs[indx,], ys[indx], family="binomial", alpha=alpha, lambda=lambda, 
                       standardize=FALSE, intercept=FALSE)
         beta <- matrix(fit$beta)
         resid <- -(ys * xs %*% beta) + log(1+exp(xs %*% beta))
         if(all(beta==0)){return(list(object=-Inf,index=indx,residu=resid,beta=beta))} 
         resid.sort <- sort(resid,decreasing=FALSE,index.return=TRUE) 
         h0 <- floor((length(y[y==0])+1)*hsize)
         h1 <- h - h0
         index0 <- resid.sort$ix[y[resid.sort$ix]==0][1:h0]
         index1 <- resid.sort$ix[y[resid.sort$ix]==1][1:h1]
         indxnew <- c(index0,index1)
         obj <- objBinom(xs, ys, beta, indxnew, alpha, lambda)
      } else {
         
         fit <- glmnet(x[indx,],y[indx],family="binomial",alpha=alpha,lambda=lambda,standardize=FALSE,intercept=FALSE)
         beta <- matrix(fit$beta)
         resid <- -(y * x %*% beta) + log(1+exp(x %*% beta))
         if(all(beta==0)){return(list(object=-Inf,index=indx,residu=resid,beta=beta))}
         resid.sort <- sort(resid,decreasing=FALSE,index.return=TRUE) 
         h0 <- floor((length(y[y==0])+1)*hsize)
         h1 <- h-h0
         index0 <- resid.sort$ix[y[resid.sort$ix]==0][1:h0]
         index1 <- resid.sort$ix[y[resid.sort$ix]==1][1:h1]
         indxnew <- c(index0,index1)
         obj <- objBinom(x, y, beta, indxnew, alpha, lambda)
      }
      return(list(object=obj, index=indxnew, residu=resid, beta=beta))
   }
#######################################################################################################
CStep.gaus <- 
   function(x, y, indx, h, hsize, alpha, lambda, scal)
   {
    
      if (scal){
         xs <- standardize.x(x, index=indx, centerFun=mean, scaleFun=sd)
         ys <- center.y(y, index=indx, centerFun=mean)
         
         fit <- glmnet(xs[indx,], ys[indx], alpha=alpha, lambda=lambda, 
                       standardize=FALSE, intercept=FALSE)
         beta <- matrix(fit$beta)
         resid <- ys - predict(fit, xs, exact=TRUE)
         resid.sort <- sort(abs(resid), index.return=TRUE)
         indxnew <- resid.sort$ix[1:h]
         obj <- objGaus(xs, ys, beta, indxnew, alpha, lambda)
      } else {
         
         fit <- glmnet(x[indx,], y[indx], alpha=alpha, lambda=lambda, 
                       standardize=FALSE, intercept=FALSE)
         beta <- matrix(fit$beta)
         resid <- y - predict(fit, x, exact=TRUE)
         resid.sort <- sort(abs(resid), index.return=TRUE)
         indxnew <- resid.sort$ix[1:h]
         obj <- objGaus(x, y, beta, indxnew, alpha, lambda)
      }
      return(list(object=obj, index=indxnew, residu=resid, beta=beta))
   }

