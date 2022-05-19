
index.subsets.gaus <- function(x, nsamp){
   is <- replicate(nsamp,sample.int(nrow(x), 3))
   return(is)
}

index.subsets.binom <- function(x, y, nsamp){
   is <- replicate(nsamp,c(sample(which(y==1),2),sample(which(y==0),2)))
   return(is)
}


## starting with 500 index subsets sometimes problematic. That is the reason we made it optional
index.subsets.multinom <- function(x, y, nsamp){
  k <- ncol(y)
  # ysave <<- y
  # xsave <<- x
  is <- replicate(nsamp,
                  drop(matrix(do.call(rbind,lapply(lapply(1:k,function(c,y){ 
                    sample(which(y[,c]==1),2)},y),matrix)),ncol=1)))
  return(is)
}

###############################################################################################################
InitialSubset.multinom <- function(x, y, h, hsize, alpha, lambda, nsamp, scal, ncores, seed) {
  # gives initial 500 subsamples after Two C Steps
  if (!is.null(seed)) set.seed(seed)
  
  index.subsets <- index.subsets.multinom(x, y, nsamp)
  
  twoCstep <- function(c, x, y, index.subsets, h, hsize, alpha, lambda, scal) {
    ## C step 1 
    Cstep1 <- CStep.multinom(x, y, index.subsets[,c], h, hsize, alpha, lambda, scal=FALSE)
    indx1 <- Cstep1$index
    ## C step 2
    Cstep2 <- CStep.multinom(x, y, indx1, h, hsize, alpha, lambda, scal) # h observations
    indx2 <- Cstep2$index
    object <- Cstep2$object
    return(list(obj=object, indx=indx2)) 
  }
  subsets <- mclapply(1:nsamp,
                      FUN = twoCstep,
                      x = x,
                      y = y,
                      index.subsets = index.subsets,
                      h = h,
                      hsize = hsize,
                      alpha = alpha,
                      lambda = lambda,
                      scal =scal,
                      mc.cores = ncores,
                      mc.cleanup = TRUE,
                      mc.allow.recursive = FALSE)
  return(subsets=subsets)
}


###############################################################################################################
InitialSubset.binom <- function(x, y, h, hsize, alpha, lambda, nsamp, scal, ncores, seed) {
   # gives initial 500 subsamples after Two C Steps
   if (!is.null(seed)) set.seed(seed)
   
   index.subsets <- index.subsets.binom(x, y, nsamp)
   
   twoCstep <- function(c, x, y, index.subsets, h, hsize, alpha, lambda, scal) {
      ## C step 1 
      Cstep1 <- CStep.binom(x, y, index.subsets[,c], h, hsize, alpha, lambda/4, scal=FALSE)
      indx1 <- Cstep1$index
      ## C step 2
      Cstep2 <- CStep.binom(x, y, indx1, h, hsize, alpha, lambda/h, scal) # h observations
      indx2 <- Cstep2$index
      object <- Cstep2$object
      return(list(obj=object, indx=indx2)) }
   
   subsets <- mclapply(1:nsamp,
                       FUN = twoCstep,
                       x = x,
                       y = y,
                       index.subsets = index.subsets,
                       h = h,
                       hsize = hsize,
                       alpha = alpha,
                       lambda = lambda,
                       scal =scal,
                       mc.cores = ncores,
                       mc.cleanup = TRUE,
                       mc.allow.recursive = FALSE)
   # return(subsets=subsets)
   return(subsets=subsets)
}

###############################################################################################################
InitialSubset.gaus <- function(x, y, h, hsize, alpha, lambda, nsamp, scal, ncores, seed) {
   # gives initial 500 subsamples after Two C Steps
   if (!is.null(seed)) set.seed(seed)
   
   index.subsets <- index.subsets.gaus(x, nsamp)
   
   twoCstep <- function(c, x, y, index.subsets, h, hsize, alpha, lambda, scal) {
      ## C step 1 
      Cstep1 <- CStep.gaus(x, y, index.subsets[,c], h, hsize, alpha, lambda/3, scal=FALSE)
      indx1 <- Cstep1$index
      ## C step 2
      Cstep2 <- CStep.gaus(x, y, indx1, h, hsize, alpha, lambda/h, scal) # h observations
      indx2 <- Cstep2$index
      object <- Cstep2$object
      return(list(obj=object, indx=indx2)) }
   
   subsets <- mclapply(1:nsamp,
                       FUN = twoCstep,
                       x = x,
                       y = y,
                       index.subsets = index.subsets,
                       h = h,
                       hsize = hsize,
                       alpha = alpha,
                       lambda = lambda,
                       scal =scal,
                       mc.cores = ncores,
                       mc.cleanup = TRUE,
                       mc.allow.recursive = FALSE)
   # return(subsets=subsets)
   return(subsets=subsets)
}

