## Author: F. Sevinc KURNAZ 
## YTU & TUwien 
enetLTS <-
   function(xx,
            yy,
            family=c("gaussian","binomial","multinomial"),
            alphas=seq(0,1,length=41),
            lambdas=NULL,
            lambdaw=NULL,
            intercept=TRUE,
            scal=TRUE,
            hsize=0.75,
            nsamp=c(500,10),
            nCsteps=20,
            nfold=5,
            repl=1,
            ncores=1,
            tol=-1e6,
            seed=NULL,
            del=0.0125,
            crit.plot=FALSE,
            typegrouped=FALSE,
            type.response=c("link","response","class"))
   {

      matchedCall      <- match.call()
      matchedCall[[1]] <- as.name("enetLTS")
      family           <- match.arg(family)
      #type <- match.arg(type)

      type.response    <- match.arg(type.response)
      if(family=="gaussian" && type.response=="class") stop("'class' is not available for gaussian")

      scal      <- isTRUE(scal)
      intercept <- isTRUE(intercept)
      plot      <- isTRUE(crit.plot)

      alphas    <- sort(alphas)
      wh        <- (alphas<0 | alphas>1)
      if (sum(wh)>0) stop("alphas can take the values only between 0 and 1")
      alphas    <- as.double(alphas)

      ncx       <- dim(xx)
      if (is.null(ncx) | (ncx[2]<=1)) stop ("X should be a matrix with 2 or more columns")

      xx        <- addColnames(as.matrix(xx))
      nobs      <- as.integer(ncx[1])  # number of observatıons
      nvars     <- as.integer(ncx[2]) # number of variables

      hsize     <- rep(hsize, length.out=1)
      if(!isTRUE(is.numeric(hsize) && 0.5 <= hsize && hsize <= 1)) {
        stop("'hsize' must be between 0.5 and 1")
      }
      h         <- floor((nobs+1)*hsize)

      yy        <- drop(yy)   # dont use matrix
      dimy      <- dim(yy)
      nrowy     <- ifelse(is.null(dimy), length(yy), dimy[1])

      if (nrowy!=nobs) stop (paste("number of observations in y (",nrowy,") not equal to the number of rows of x (",nobs,")",sep=""))

      if (repl<=0) stop ("repl has to be a positive number")
      if (nCsteps<=0) stop ("nCsteps has to be a positive number")

      nsamp <- rep(nsamp, length.out=2)
      if(!is.numeric(nsamp) || any(!is.finite(nsamp))) {
        nsamp <- formals$nsamp()
        warning("missing or infinite values in 'nsamp'; using default values")
      } else nsamp <- as.integer(nsamp)
      s1 <- nsamp[2]  # the number of subsamples to keep after the first phase
      nsamp <- nsamp[1] # the number of initial subsamples

      ncores <- rep(ncores, length.out=1)
      if (is.na(ncores)) ncores <- detectCores()  # use all available cores
      if (!is.numeric(ncores) || is.infinite(ncores) || ncores < 1) {
         ncores <- 1  # use default value
         warning ("invalid value of 'ncores'; using default value")
      } else ncores <- as.integer(ncores)
      # check whether parallel computing should be used
      haveNCores <- ncores > 1
      if (haveNCores) {
         if (.Platform$OS.type == "windows") {
            cl <- makePSOCKcluster(rep.int("localhost", ncores))
         } else cl <- makeForkCluster(ncores)
         on.exit(stopCluster(cl))
      }

      tol <- rep(tol, length.out=1)
      if(!is.numeric(tol) || !is.finite(tol)) {
        tol <- formals()$tol
        warning("missing or infinite value of 'tol'; using default value")
      }

      fit <- switch(family,
                    "gaussian"    = enetLTS.gaus(xx, yy, alphas, lambdas, lambdaw, h, hsize, nobs, nvars, intercept,
                                              nsamp, s1, nfold, repl, scal, ncores, nCsteps, tol, seed, del, plot,
                                              type.response),
                    "binomial"    = enetLTS.binom(xx, yy, alphas, lambdas, lambdaw, h, hsize, nobs, nvars, intercept,
                                               nsamp, s1, nfold, repl, scal, ncores, nCsteps, tol, seed, del, plot,
                                               type.response),
                    "multinomial" = enetLTS.multinom(xx, yy, alphas, lambdas, lambdaw, h, hsize, nobs, nvars, intercept,
                                                     nsamp, s1, nfold, repl, scal, ncores, nCsteps, tol, seed, plot, typegrouped,
                                                     type.response)

                     )
      output <- fit
      class(output) <- "enetLTS"

      output$call   <- matchedCall
      output
   }


