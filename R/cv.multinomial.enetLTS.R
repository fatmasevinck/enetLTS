

cv.multinom.enetLTS <- function(index=NULL, xx, yy, alphas, lambdas,
                                nfold, repl, ncores, plot=TRUE){

  RTMSPE <- RMSPE <- TMNLL <- MNLL <- NULL
   # return only smallest fraction 1-trim of values:
   uptrim <- function(x,trim=0.1) {
      return(sort(x)[1:(length(x)*(1-trim))])
   }
   family <- "multinomial"
   k <- dim(yy)[2]
   n <- dim(xx)[1]

   wh <- (alphas<0 | alphas>1)
   if (sum(wh)>0)        stop("alphas can take the values only between 0 and 1")
   if (missing(alphas))  stop("provide an alphas sequence")
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
            y <- yy[index[[i]][,j],]
         }
      evalCritl <- rep(NA,repl)
      for (l in 1:repl){ # repeate CV
         folds_k <- lapply(1:k, function(c,y,nfold){cvFolds(length(y[,c][y[,c]==1]),K=nfold,R=1,type="random")},y,nfold)
        # print(folds_k)
         loss <- lapply(1:k, function(c,y){rep(NA,length(y[,c][y[,c]==1]))},y)
         for (f in 1:nfold) {

            xtrain_k <- lapply(1:k, function(c,x,y,f) {
               x[which(y[,c]==1),][folds_k[[c]]$subsets[folds_k[[c]]$which != f,1], ]
            },x,y,f)
            xtrain <- do.call(rbind,xtrain_k)

            xtest_k <- lapply(1:k, function(c,x,y,f) {
               x[which(y[,c]==1),][folds_k[[c]]$subsets[folds_k[[c]]$which == f,1], ]
            },x,y,f)
            xtest <- do.call(rbind,xtest_k)

            ytrain_k <- lapply(1:k, function(c,y,f) {
               y[which(y[,c]==1),][folds_k[[c]]$subsets[folds_k[[c]]$which != f,1], ]
            },y,f)
            ytrain <- do.call(rbind,ytrain_k)

            ytest_k <- lapply(1:k, function(c,y,f) {
               y[which(y[,c]==1),][folds_k[[c]]$subsets[folds_k[[c]]$which == f,1], ]
            },y,f)
            ytest <- do.call(rbind,ytest_k)

            res <- tryCatch({
               trainmod <- glmnet(xtrain,ytrain,family,alpha=alpha,lambda=lambda,
                                  standardize=FALSE,intercept=FALSE)},error=function(err){
                                     error <- TRUE
                                     return(error)
                                  })
            if (is.logical(res)){
               print(paste("CV broke off for alpha=",alpha ,"and lambda=", lambda))
            } else {
               trainmod <- res

               loss_k <- lapply(1:k, function(c,xtest_k,ytest_k,trainmod,f){
                  (-apply(ytest_k[[c]]*log(drop(predict(trainmod,newx=xtest_k[[c]],type="response"))),1,sum))
               },xtest_k,ytest_k,trainmod,f)

               for (ii in 1:k){
                  loss[[ii]][folds_k[[ii]]$which == f ] <- loss_k[[ii]]
               }
            }
         }

         evalgroups <- vector("list",k)
         for (i in 1:k){
            evalgroups[[i]] <- uptrim(loss[[i]])
         }

         evalCritl[l] <- mean(unlist(evalgroups),na.rm=TRUE) # that works if we have repl more than 1, now only mean
      }
      evalCrit <- mean(evalCritl,na.rm=TRUE)
      return(evalCrit = evalCrit)
   }
   temp_result <- mclapply(indcombi,
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

   evalCrit <- matrix(unlist(temp_result), ncol = length(lambdas) , byrow = FALSE) # row alphas, col lambdas
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
            geom_tile() + scale_fill_gradientn(colours=mycol.b) + theme(axis.text.x=element_text(angle=-90))
         mspeplot <- mspeplot + ggtitle(paste0("TMNLL (optimal at lambda=",lambda,",alpha=",alpha,",  ",family,")"))
      } else {
         names(ggmspe) <- c("lambda","alpha","MNLL")
         mspeplot <- ggplot(ggmspe,aes(x=as.factor(alpha),y=as.factor(lambda),fill=MNLL)) +
            geom_tile() + scale_fill_gradientn(colours=mycol.b) + theme(axis.text.x=element_text(angle=-90))
         mspeplot <- mspeplot + ggtitle(paste0("MNLL (optimal at lambda=",lambda,",alpha=",alpha,",",family,")"))
      }
      mspeplot <- mspeplot + xlab("lambda") + ylab("alpha")
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(1,1)))
      print(mspeplot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
   }
   return(list(indexbest=indexbest,evalCrit=evalCrit,minevalCrit=minevalCrit,lambdaopt=lambda,alphaopt=alpha))
}

