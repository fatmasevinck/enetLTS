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
            type.response=c("link","response","class"),
            hgrid=c(0.65,0.75,0.85,0.95),
            screen=TRUE,
            crossfit=NULL)
   {

      matchedCall      <- match.call()
      matchedCall[[1]] <- as.name("enetLTS")
      family           <- match.arg(family)

      # hsize = "optimal": fit for every value of 'hgrid' and keep the best one
      if (is.character(hsize)) {
        match.arg(hsize, "optimal")
        return(enetLTSoptimal(matchedCall, parent.frame(), hgrid))
      }


      type.response    <- match.arg(type.response)
      if (isTRUE(crossfit) && family == "binomial") {
        warning("'crossfit' is available for the gaussian and multinomial families only and is ignored")
        crossfit <- FALSE
      }
      if(family=="gaussian" && type.response=="class") stop("'class' is not available for gaussian")

      scal      <- isTRUE(scal)
      # default: cross-fitted trimming for the multinomial family with scal = TRUE, off otherwise
      if (is.null(crossfit)) crossfit <- (family == "multinomial" && scal)
      intercept <- isTRUE(intercept)
      plot      <- isTRUE(crit.plot)

      alphas    <- sort(alphas)
      wh        <- (alphas<0 | alphas>1)
      if (sum(wh)>0) stop("alphas can take the values only between 0 and 1")
      alphas    <- as.double(alphas)

      ncx       <- dim(xx)
      if (is.null(ncx) || (ncx[2]<=1)) stop ("X should be a matrix with 2 or more columns")
      if (anyNA(xx) || anyNA(yy)) stop("missing values in 'xx' or 'yy' are not allowed")

      xx        <- addColnames(as.matrix(xx))
      if (family == "binomial" && is.factor(yy)) {
        if (!all(levels(yy) %in% c("0", "1")))
          stop("for family = \"binomial\", 'yy' must be coded 0/1")
        yy <- as.numeric(as.character(yy))
      }
      nobs      <- as.integer(ncx[1])  # number of observations
      nvars     <- as.integer(ncx[2]) # number of variables

      hsize     <- rep(hsize, length.out=1)
      if(!isTRUE(is.numeric(hsize) && 0.5 <= hsize && hsize <= 1)) {
        stop("'hsize' must be between 0.5 and 1")
      }
      h         <- floor((nobs+1)*hsize)

      yy        <- drop(yy)   # a one-column matrix becomes a vector
      dimy      <- dim(yy)
      nrowy     <- ifelse(is.null(dimy), length(yy), dimy[1])

      if (nrowy!=nobs) stop (paste("number of observations in y (",nrowy,") not equal to the number of rows of x (",nobs,")",sep=""))

      if (repl<=0) stop ("repl has to be a positive number")
      if (nCsteps<=0) stop ("nCsteps has to be a positive number")

      nsamp <- rep(nsamp, length.out=2)
      if(!is.numeric(nsamp) || any(!is.finite(nsamp))) {
        nsamp <- c(500, 10)
        warning("missing or infinite values in 'nsamp'; using default values")
      } else nsamp <- as.integer(nsamp)
      s1 <- nsamp[2]  # the number of subsamples to keep after the first phase
      nsamp <- nsamp[1] # the number of initial subsamples

      ncores <- checkNcores(ncores)
      startWorkers(ncores)
      on.exit(stopWorkers())

      tol <- rep(tol, length.out=1)
      if(!is.numeric(tol) || !is.finite(tol)) {
        tol <- -1e6
        warning("missing or infinite value of 'tol'; using default value")
      }

      fit <- switch(family,
                    "gaussian"    = enetLTS.gaus(xx, yy, alphas, lambdas, lambdaw, h, hsize, nobs, nvars, intercept,
                                              nsamp, s1, nfold, repl, scal, ncores, nCsteps, tol, seed, del, plot,
                                              type.response, crossfit),
                    "binomial"    = enetLTS.binom(xx, yy, alphas, lambdas, lambdaw, h, hsize, nobs, nvars, intercept,
                                               nsamp, s1, nfold, repl, scal, ncores, nCsteps, tol, seed, del, plot,
                                               type.response),
                    "multinomial" = enetLTS.multinom(xx, yy, alphas, lambdas, lambdaw, h, hsize, nobs, nvars, intercept,
                                                     nsamp, s1, nfold, repl, scal, ncores, nCsteps, tol, seed, plot, typegrouped,
                                                     type.response, screen, crossfit, del)

                     )
      output <- fit
      class(output) <- "enetLTS"

      output$call   <- matchedCall
      output
   }

