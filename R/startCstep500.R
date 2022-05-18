
startCstep500.multinom <- 
  function(xx, yy, h, hsize, alpha, lambda, nsamp, s1, scal, ncores, nCsteps, tol, seed)
  {
    ## internal function warmCsteps
    H2 <- selectbest10.multinom(xx, yy, h, hsize, alpha, lambda, nsamp, s1, scal, ncores, seed)
    best10subsets <- H2$idxbest
    s1 <- H2$s1_new
    lastbestindex <- mclapply(1:s1, function(zz, xx, yy, h, hsize, alpha, lambda, nCsteps, tol, best10subsets) {
      indexsubbest <- best10subsets[[zz]]
      objbest <- tol
      cstep.mod <- CStep.multinom(xx, yy, indexsubbest, h, hsize, alpha, lambda, scal)
      countloop <- 0
      while ((cstep.mod$object>objbest) & (countloop<nCsteps)){ 
        objbest <- cstep.mod$object 
        newindex <- cstep.mod$index  
        cstep.mod <- CStep.multinom(xx, yy, newindex, h, hsize, alpha, lambda, scal)
      }
      return(list(lastindex=newindex, objbest=objbest))
    }, 
    xx = xx, 
    yy = yy, 
    h = h, 
    hsize = hsize,
    alpha = alpha, 
    lambda = lambda, 
    nCsteps = nCsteps, 
    tol = tol, 
    best10subsets = best10subsets, 
    mc.cores = ncores, 
    mc.allow.recursive=FALSE) 
    
    obj <- NULL
    for (i in 1:s1){
      obj <- c(obj, lastbestindex[[i]]$objbest)
    }
    whichbestindex <- sort(obj, decreasing=TRUE, index.return=TRUE)$ix[1]   # NOTE !!!!!!!!!!!!
    index <- lastbestindex[[whichbestindex]]$lastindex
    
    return(index=index) 
  }

############################################################################################################
selectbest10.multinom <- function(x, y, h, hsize, alpha, lambda, nsamp, s1, scal, ncores, seed) {
  obj <- NULL
  subsets <- InitialSubset.multinom(x, y, h, hsize, alpha, lambda, nsamp, scal, ncores, seed)
  
  obj <- unlist(mclapply(1:nsamp, function(ob, sub){
    ob_val <- subsets[[ob]]$obj
  }, subsets, mc.cores=ncores, mc.allow.recursive=FALSE))
  
  obj_sorted <- sort(obj,decreasing=FALSE,index.return=TRUE)   # NOTE !!!!!
  
  obj <- obj_sorted$x[1:s1]
  s1_new <- length(obj[!is.infinite(obj)])
  idx <- obj_sorted$ix[1:s1_new]
  if (s1_new==0){
    stop(paste("Model is not suitable for alpha",alpha,"lambda",lambda,"for this data set. Choose another lambda."))
  }
  bestindex <- mclapply(1:s1_new, function(c, idx, subsets) {
    indx <- subsets[[idx[c]]]$indx
  }, idx, subsets, mc.cores = ncores)
  
  return(list(idxbest=bestindex, s1_new=s1_new, subsets=subsets))
}
################################################################################################

################################################################################################

startCstep500.binom <- 
   function(xx, yy, h, hsize, alpha, lambda, nsamp, s1, scal, ncores, nCsteps, tol, seed)
   {
      ## internal function warmCsteps

      H2 <- selectbest10.binom(xx, yy, h, hsize, alpha, lambda, nsamp, s1, scal, ncores, seed)
      best10subsets <- H2$idxbest
      s1 <- H2$s1_new
      
      lastbestindex <- mclapply(1:s1, function(zz, xx, yy, h, hsize, alpha, lambda, nCsteps, tol, best10subsets) {
         indexsubbest <- best10subsets[[zz]]
         objbest <- tol
         cstep.mod <- CStep.binom(xx, yy, indexsubbest, h, hsize, alpha, lambda/h, scal)
         countloop <- 0
         while ((cstep.mod$object>objbest) & (countloop<nCsteps)){ 
            objbest <- cstep.mod$object 
            newindex <- cstep.mod$index  
            cstep.mod <- CStep.gaus(xx, yy, newindex, h, hsize, alpha, lambda/h, scal)
         }
         return(list(lastindex=newindex, objbest=objbest))
      }, 
      xx = xx, 
      yy = yy, 
      h = h, 
      hsize = hsize,
      alpha = alpha, 
      lambda = lambda, 
      nCsteps = nCsteps, 
      tol = tol, 
      best10subsets = best10subsets, 
      mc.cores = ncores, 
      mc.allow.recursive=FALSE) 
      
      obj <- NULL
      for (i in 1:s1){
         obj <- c(obj, lastbestindex[[i]]$objbest)
      }
      whichbestindex <- sort(obj, decreasing=TRUE, index.return=TRUE)$ix[1]   # NOTE !!!!!!!!!!!!
      index <- lastbestindex[[whichbestindex]]$lastindex
      
      return(index=index) 
   }

############################################################################################################
selectbest10.binom <- function(x, y, h, hsize, alpha, lambda, nsamp, s1, scal, ncores, seed) {
   obj <- NULL
   subsets <- InitialSubset.binom(x, y, h, hsize, alpha, lambda, nsamp, scal, ncores, seed)
   
   obj <- unlist(mclapply(1:nsamp, function(ob, sub){
      ob_val <- subsets[[ob]]$obj
   }, subsets, mc.cores=ncores, mc.allow.recursive=FALSE))
   
   obj_sorted <- sort(obj,decreasing=TRUE,index.return=TRUE)   # NOTE !!!!!
   
   obj <- obj_sorted$x[1:s1]
   s1_new <- length(obj[!is.infinite(obj)])
   idx <- obj_sorted$ix[1:s1_new]
   if (s1_new==0){
      stop(paste("Model is not suitable for alpha",alpha,"lambda",lambda,"for this data set. Choose another lambda."))
   }
   bestindex <- mclapply(1:s1_new, function(c, idx, subsets) {
      indx <- subsets[[idx[c]]]$indx
   }, idx, subsets, mc.cores = ncores)
   
   return(list(idxbest=bestindex, s1_new=s1_new, subsets=subsets))
   # return(idxbest=bestindex)
}

################################################################################################


startCstep500.gaus <- 
   function(xx, yy, h, hsize, alpha, lambda, nsamp, s1, scal, ncores, nCsteps, tol, seed)
   {
      ## internal function warmCsteps
     
      #  source("objectiveFunc.R")
      #  source("InitialSubsets.R")
      #  source("Csteps.R")
      #  source("utilities.R")
      
      H2 <- selectbest10.gaus(xx, yy, h, hsize, alpha, lambda, nsamp, s1, scal, ncores, seed)
      best10subsets <- H2$idxbest
      s1 <- H2$s1_new
    
         lastbestindex <- mclapply(1:s1, function(zz, xx, yy, h, hsize, alpha, lambda, nCsteps, tol, best10subsets) {
            indexsubbest <- best10subsets[[zz]]
            objbest <- tol
            cstep.mod <- CStep.gaus(xx, yy, indexsubbest, h, hsize, alpha, lambda/h, scal)
            countloop <- 0
            while ((cstep.mod$object>objbest) & (countloop<nCsteps)){ 
               objbest <- cstep.mod$object 
               newindex <- cstep.mod$index  
               cstep.mod <- CStep.gaus(xx, yy, newindex, h, hsize, alpha, lambda/h, scal)
            }
            return(list(lastindex=newindex, objbest=objbest))
         }, 
         xx=xx, 
         yy=yy, 
         h=h, 
         hsize=hsize,
         alpha=alpha, 
         lambda=lambda, 
         nCsteps=nCsteps, 
         tol=tol, 
         best10subsets=best10subsets, 
         mc.cores=ncores, 
         mc.allow.recursive=FALSE) 
      
      obj <- NULL
      for (i in 1:s1){
         obj <- c(obj, lastbestindex[[i]]$objbest)
      }
      whichbestindex <- sort(obj, decreasing=TRUE, index.return=TRUE)$ix[1]
      index <- lastbestindex[[whichbestindex]]$lastindex
      
      return(index=index) 
   }

############################################################################################################
selectbest10.gaus <- function(x, y, h, hsize, alpha, lambda, nsamp, s1, scal, ncores, seed) {
   obj <- NULL
   subsets <- InitialSubset.gaus(x, y, h, hsize, alpha, lambda, nsamp, scal, ncores, seed)

   obj <- unlist(mclapply(1:nsamp, function(ob, sub){
      ob_val <- subsets[[ob]]$obj
   }, subsets, mc.cores=ncores, mc.allow.recursive=FALSE))
   
   obj_sorted <- sort(obj, decreasing=FALSE, index.return=TRUE)
   
   obj <- obj_sorted$x[1:s1]
   s1_new <- length(obj[!is.infinite(obj)])
   idx <- obj_sorted$ix[1:s1_new]
   if (s1_new==0){
      stop(paste("Model is not suitable for alpha",alpha,"lambda",lambda,"for this data set. Choose another lambda."))
   }
   bestindex <- mclapply(1:s1_new, function(c, idx, subsets) {
      indx <- subsets[[idx[c]]]$indx
   }, idx, subsets, mc.cores = ncores)
   
   return(list(idxbest=bestindex, s1_new=s1_new, subsets=subsets))
   # return(idxbest=bestindex)
}
