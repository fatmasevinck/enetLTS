
cv.gaussian.enetLTS <- function(index=NULL, xx, yy, alphas, lambdas, 
                            nfold, repl, ncores, plot=TRUE){

   family <- "gaussian"
   RTMSPE <- RMSPE <- NULL
   wh <- (alphas<0 | alphas>1)
   if (sum(wh)>0) stop ("alphas can take the values only between 0 and 1")
   if (missing(alphas)) stop ("provide an alphas sequence")
   if (missing(lambdas)) stop ("provide an lambdas sequence")
   
   combis_ind <- expand.grid(1:length(alphas), 1:length(lambdas))
   indcombi <- 1:nrow(combis_ind)
   
   calc_evalCrit <- function(rowind, combis_ind, alphas, lambdas, 
                             index, xx, yy, nfold, repl){
      
      i <- combis_ind[rowind, 1] # alphas
      j <- combis_ind[rowind, 2] # lambdas
      alpha <- alphas[i]
      lambda <- lambdas[j]
      
      if (is.null(index)) {
         x <- xx
         y <- yy } else {
            x <- xx[index[[i]][,j],]
            y <- yy[index[[i]][,j]]
         }
      evalCritl <- rep(NA,repl)
      for (l in 1:repl){ # repeate CV
         folds <- cvFolds(length(y), K = nfold, R = 1, type = "random")
         loss <- rep(NA,nrow(x))
         for (f in 1:nfold) {
            xtrain <- x[folds$subsets[folds$which != f, 1], ]
            ytrain <- y[folds$subsets[folds$which != f, 1] ]
            xtest <- x[folds$subsets[folds$which == f, 1], ]
            ytest <- y[folds$subsets[folds$which == f, 1] ]
            res <- tryCatch({
               hpen <- length(ytrain)
               trainmod <- glmnet(xtrain, ytrain, alpha=alpha, lambda=lambda/hpen,
                                  standardize=FALSE, intercept=FALSE)}, error=function(err){
                                     error <- TRUE
                                     return(error)
                                  })
            if (is.logical(res)){
               print(paste("CV broke off for alpha=", alpha, "and lambda=", lambda))
            } else {
               trainmod <- res
               loss[folds$which == f ] <- ytest - xtest %*% matrix(trainmod$beta)
            }
         }
         if (is.null(index)){
            evalCritl[l] <- sqrt(uptrimCrit((loss)^2))
         } else {
            evalCritl[l] <- sqrt(mean((loss)^2))
         }
      }
      evalCrit <- mean(evalCritl,na.rm=TRUE)
      return(list(evalCrit = evalCrit))
   }
   temp_result <- mclapply(1:nrow(combis_ind),
                           FUN = calc_evalCrit,
                           combis_ind = combis_ind,
                           alphas = alphas,
                           lambdas = lambdas,
                           index = index,
                           xx = xx,
                           yy = yy,
                           nfold = nfold,
                           repl = repl,
                           mc.cores = ncores,
                           mc.allow.recursive = FALSE)
   
   evalCrit <- matrix(unlist(temp_result), ncol = length(lambdas), byrow = FALSE) # row alphas, col lambdas
   dimnames(evalCrit) <- list(paste("alpha", alphas), paste("lambda", lambdas))
   
   optind <- which(evalCrit == min(evalCrit, na.rm = TRUE), arr.ind = TRUE)[1, ]
   alpha_optind <- unlist(optind[1])  # row
   lambda_optind <- unlist(optind[2])  # col
   minevalCrit <- evalCrit[alpha_optind,lambda_optind]
   indexbest <- index[[alpha_optind]][,lambda_optind]
   alphas <- round(alphas,4)
   lambdas <- round(lambdas,4)
   alpha <- alphas[alpha_optind]
   lambda <- lambdas[lambda_optind]
   
   if (plot==TRUE){
      print(paste("optimal model: lambda =", lambda, "alpha =", alpha))
      lenCol <- length(alphas)*length(lambdas)
      mycol.b <- colorRampPalette(c("black", "blue2", "purple", "orange", "yellow"))(lenCol)
      
      ggmspe <- evalCrit
      rownames(ggmspe) <- alphas
      colnames(ggmspe) <- lambdas
      ggmspe <- melt(ggmspe)
      if (is.null(index)){
         names(ggmspe) <- c("lambda", "alpha", "RTMSPE")
         mspeplot <- ggplot(ggmspe, aes(x=as.factor(alpha), y=as.factor(lambda), fill=RTMSPE)) +
            geom_tile() +  scale_fill_gradientn(colours=mycol.b) + theme(axis.text.x=element_text(angle=-90))
         mspeplot <- mspeplot + ggtitle(paste0("RTMSPE (minimum at lambda=", lambda,",alpha=",alpha,", ",family,")"))
      } else {
         names(ggmspe) <- c("lambda", "alpha", "RMSPE")
         mspeplot <- ggplot(ggmspe, aes(x=as.factor(alpha), y=as.factor(lambda), fill=RMSPE)) +
            geom_tile() +  scale_fill_gradientn(colours=mycol.b) + theme(axis.text.x=element_text(angle=-90))
         mspeplot <- mspeplot + ggtitle(paste0("RMSPE (minimum at lambda=", lambda, ",alpha=", alpha, ", ", family,")"))
      }
      mspeplot <- mspeplot + xlab("lambda") +  ylab("alpha")
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(1,1)))
      print(mspeplot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
   }
   
   return(list(indexbest=indexbest,evalCrit=evalCrit,minevalCrit=minevalCrit,lambdaopt=lambda,alphaopt=alpha))
}



