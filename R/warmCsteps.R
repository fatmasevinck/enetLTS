

warmCsteps.mod <- 
   function(x, y, h, hsize, alphas, lambdas, nsamp, s1, 
            scal, ncores, nCsteps, tol, seed, family){
      
      alpha <- alphas[1]
     
      if (family == "gaussian"){
        lambda <- lambdas[1] 
        startCstep500.fun <- startCstep500.gaus
        CStep.fun <- CStep.gaus
      } else if (family == "binomial"){
        lambda <- lambdas[1] 
        startCstep500.fun <- startCstep500.binom
        CStep.fun <- CStep.binom
      } else if (family == "multinomial"){
        lambda <- lambdas[1] 
        # ll <- length(lambdas) ### PF changed
        # lambda <- lambdas[ceiling(ll/2)] ### PF changed from lambdas[1] 
        startCstep500.fun <- startCstep500.multinom
        CStep.fun <- CStep.multinom
      } 
      
      startCstep500 <- startCstep500.fun(x, y, h, hsize, alpha, lambda, nsamp, s1,
                                         scal, ncores, nCsteps, tol, seed)
      
      if (length(alphas)==1 & length(lambdas)==1) {
         indexall[,1,1] <- startCstep500
         return(indexall=indexall)
      } else if (length(alphas)>1 & length(lambdas)==1) {
         firstindex <- startCstep500
         indexall <- run.lambda.alphas(x, y, h, hsize, alphas, lambdas, firstindex, scal, nCsteps, tol, CStep.fun)
         return(indexall=indexall)
      } else if (length(alphas)==1 & length(lambdas)>1) {
         firstindex <- startCstep500
         indexall <- run.alpha.lambdas(x, y, h, hsize, alphas, lambdas, firstindex, scal, nCsteps, tol, CStep.fun)
         return(indexall=indexall)
      } else if (length(alphas)>1 & length(lambdas)>1) {
         firstindex <- startCstep500
         idx.alphas.lambda <- run.lambda.alphas(x, y, h, hsize, alphas, lambdas[1], 
                                                firstindex, scal, nCsteps, tol, CStep.fun)
         temp_result <- mclapply(1:length(alphas),
                                 FUN = run.sngl.alpha,
                                 x = x,
                                 y = y,
                                 h = h,
                                 hsize = hsize,
                                 alphas = alphas,
                                 lambdas = lambdas,
                                 idx.alphas.lambda = idx.alphas.lambda,
                                 scal = scal, 
                                 nCsteps = nCsteps,
                                 tol = tol,
                                 CStep.fun = CStep.fun,
                                 mc.cores = ncores,
                                 mc.cleanup = TRUE,
                                 mc.allow.recursive = FALSE)
         
         indexall <- temp_result
         return(indexall)
      }
   }

run.sngl.alpha <- function(cc, x, y, h, hsize, alphas, lambdas, idx.alphas.lambda, 
                           scal, nCsteps, tol, CStep.fun) {
   
   alpha <- alphas[cc]
   index.alphas.lambda <- idx.alphas.lambda[,cc]
   IndexMatrix <- matrix(NA, nrow=h, ncol=(length(lambdas)-1))
   
   for (ii in 1:(length(lambdas)-1)) {
      lambda <- lambdas[ii+1] 
      newindex_la <- index.alphas.lambda
      objbest <- tol 
      cstep.mod <- CStep.fun(x, y, newindex_la, h, hsize, alpha, lambda, scal)
      countloop <- 0
      # Make sure at least one additional run has been done in order to ensure existence of all objects
      countloop <- countloop + 1
      objbest <- cstep.mod$object 
      newindex_la <- cstep.mod$index  
      cstep.mod <- CStep.fun(x, y, newindex_la, h, hsize, alpha, lambda, scal)
      index1_la <- newindex_la
      while ((cstep.mod$object > objbest) & (countloop < nCsteps)) { 
         countloop <- countloop + 1
         objbest <- cstep.mod$object 
         newindex_la <- cstep.mod$index  
         cstep.mod <- CStep.fun(x, y, newindex_la, h, hsize, alpha, lambda, scal)
         index1_la <- newindex_la
      }
      IndexMatrix[,ii] <- newindex_la[1:h]   
   }
   IndexMatrix <- cbind(index.alphas.lambda, IndexMatrix)
   return(IndexMatrix)
}

run.lambda.alphas <- function(x, y, h, hsize, alphas, lambdas, index1_al, scal, nCsteps, tol, CStep.fun) {
   lambda <- lambdas
   index1_la <- index1_al
   IndexMatrix <- matrix(NA,nrow=h,ncol=(length(alphas)-1))
   for (ii in 1:(length(alphas)-1)){
      alpha <- alphas[ii+1] 
      newindex_la <- index1_la
      objbest <- tol 
      cstep.mod <- CStep.fun(x, y, newindex_la, h, hsize, alpha, lambda, scal)
      countloop <- 0
      countloop <- countloop + 1
      objbest <- cstep.mod$object 
      newindex_la <- cstep.mod$index  
      cstep.mod <- CStep.fun(x, y, newindex_la, h, hsize, alpha, lambda, scal)
      index1_la <- newindex_la
      while ((cstep.mod$object > objbest) & (countloop < nCsteps)){ 
         countloop <- countloop + 1
         objbest <- cstep.mod$object 
         newindex_la <- cstep.mod$index  
         cstep.mod <- CStep.fun(x, y, newindex_la, h, hsize, alpha, lambda, scal)
         index1_la <- newindex_la
      } 
      IndexMatrix[,ii] <- newindex_la[1:h]   
   }
   IndexMatrix <- cbind(index1_al, IndexMatrix)
   return(IndexMatrix)
}

run.alpha.lambdas <- function(x, y, h, hsize, alphas, lambdas, index1_al, scal, nCsteps, tol, CStep.fun) {
   alpha <- alphas
   index1_la <- index1_al
   IndexMatrix <- matrix(NA,nrow=h,ncol=(length(lambdas)-1))
   for (ii in 1:(length(lambdas)-1)){
      lambda <- lambdas[ii+1] 
      newindex_la <- index1_la
      objbest <- tol 
      cstep.mod <- CStep.fun(x, y, newindex_la, h, hsize, alpha, lambda, scal)
      countloop <- 0
      countloop <- countloop + 1
      objbest <- cstep.mod$object 
      newindex_la <- cstep.mod$index  
      cstep.mod <- CStep.fun(x, y, newindex_la, h, hsize, alpha, lambda, scal)
      index1_la <- newindex_la
      while ((cstep.mod$object > objbest) & (countloop < nCsteps)){ 
         countloop <- countloop + 1
         objbest <- cstep.mod$object 
         newindex_la <- cstep.mod$index  
         cstep.mod <- CStep.fun(x, y, newindex_la, h, hsize, alpha, lambda, scal)
         index1_la <- newindex_la
      }
      IndexMatrix[,ii] <- newindex_la[1:h]   
   }
   IndexMatrix <- cbind(index1_al,IndexMatrix)
   return(IndexMatrix)
}

