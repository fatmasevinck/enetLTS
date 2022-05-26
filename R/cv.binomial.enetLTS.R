
cv.binomial.enetLTS <- function(index=NULL, xx, yy, alphas, lambdas,
                             nfold, repl, ncores, plot=TRUE){
   family <- "binomial"
   RTMSPE <- RMSPE <- TMNLL <- MNLL <- NULL
   wh <- (alphas<0 | alphas>1)
   if (sum(wh)>0) stop("alphas can take the values only between 0 and 1")
   if (missing(alphas)) stop("provide an alphas sequence")
   if (missing(lambdas)) stop("provide an lambdas sequence")

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
         folds0 <- cvFolds(length(y[y==0]), K = nfold, R = 1, type = "random")
         folds1 <- cvFolds(length(y[y==1]), K = nfold, R = 1, type = "random")
         loss0 <- rep(NA,sum(y==0))
         loss1 <- rep(NA,sum(y==1))
         for (f in 1:nfold) {
            xtrain0 <- x[y==0,][folds0$subsets[folds0$which != f,1], ]
            ytrain0 <- y[y==0][folds0$subsets[folds0$which != f,1] ]
            xtest0 <- x[y==0,][folds0$subsets[folds0$which == f,1], ]
            ytest0 <- y[y==0][folds0$subsets[folds0$which == f,1] ]
            xtrain1 <- x[y==1,][folds1$subsets[folds1$which != f,1], ]
            ytrain1 <- y[y==1][folds1$subsets[folds1$which != f,1] ]
            xtest1 <- x[y==1,][folds1$subsets[folds1$which == f,1], ]
            ytest1 <- y[y==1][folds1$subsets[folds1$which == f,1] ]
            xtrain <- rbind(xtrain0,xtrain1); ytrain <- c(ytrain0,ytrain1)
            xtest <- rbind(xtest0,xtest1); ytest <- c(ytest0,ytest1)
            res <- tryCatch({
               hpen <- length(ytrain)
               trainmod <- glmnet(xtrain, ytrain, family=family, alpha=alpha, lambda=lambda/hpen,
                                  standardize=FALSE, intercept=FALSE)}, error=function(err){
                                     error <- TRUE
                                     return(error)
                                  })
            if (is.logical(res)){
               print(paste("CV broke off for alpha=", alpha, "and lambda=", lambda))
            } else {
               trainmod <- res
               loss0[folds0$which == f ] <- -(ytest0 * xtest0 %*% matrix(trainmod$beta)) +
                  log(1+exp(xtest0 %*% matrix(trainmod$beta)))
               loss1[folds1$which == f ] <- -(ytest1 * xtest1 %*% matrix(trainmod$beta)) +
                  log(1+exp(xtest1 %*% matrix(trainmod$beta)))
            }
         }
         loss <- c(loss0,loss1)
         if (is.null(index)){
            evalCritl[l] <- uptrimCrit(loss)
         } else {
            evalCritl[l] <- mean(loss,na.rm=TRUE)
         }
      }
      evalCrit <- mean(evalCritl, na.rm=TRUE)
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
      mycol.b <- colorRampPalette(c("black","blue2", "purple", "orange", "yellow"))(lenCol)

      ggmspe <- evalCrit
      rownames(ggmspe) <- alphas
      colnames(ggmspe) <- lambdas
      ggmspe <- melt(ggmspe)
      if (is.null(index)){

         names(ggmspe) <- c("lambda","alpha","TMNLL")
         mspeplot <- ggplot(ggmspe,aes(x=as.factor(alpha),y=as.factor(lambda),fill=TMNLL)) +
            geom_tile() +  scale_fill_gradientn(colours=mycol.b) + theme(axis.text.x=element_text(angle=-90))
         mspeplot <- mspeplot + ggtitle(paste0("TMNLL (optimal at lambda=",lambda,",alpha=",alpha,",  ",family,")"))
      } else {
         names(ggmspe) <- c("lambda","alpha","MNLL")
         mspeplot <- ggplot(ggmspe,aes(x=as.factor(alpha),y=as.factor(lambda),fill=MNLL)) +
            geom_tile() +  scale_fill_gradientn(colours=mycol.b) + theme(axis.text.x=element_text(angle=-90))
         mspeplot <- mspeplot + ggtitle(paste0("MNLL (optimal at lambda=",lambda,",alpha=",alpha,",  ",family,")"))
      }
      mspeplot <- mspeplot + xlab("lambda") +  ylab("alpha")
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(1,1)))
      print(mspeplot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
   }
   return(list(indexbest=indexbest,
               evalCrit=evalCrit,
               minevalCrit=minevalCrit,
               lambdaopt=lambda,
               alphaopt=alpha)
          )
}

