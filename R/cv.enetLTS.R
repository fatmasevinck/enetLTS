

cv.enetLTS <- function(index=NULL, family=c("gaussian","binomial","multinomial"), xx, yy, alphas, lambdas,
                       nfold=5, repl=1, ncores=NA, plot=TRUE)
{
  #  require(ggplot2)
  family <- match.arg(family)

  if (!is.null(index) & !is.array(index)){
    stop("index has to be an array")
  }

  alphas <- sort(alphas)
  wh <- (alphas<0 | alphas>1)
  if (sum(wh)>0) stop("alphas can take the values only between 0 and 1")
  alphas <- as.double(alphas)

  if (missing(lambdas) & family=="gaussian"){
    l0 <- robustHD::lambda0(xx, yy, normalize=TRUE, intercept=TRUE)
    lambdas <- seq(l0, 0, by=-0.025*l0)
  } else if (missing(lambdas) & family=="binomial"){
    l0 <- lambda00(xx, yy, normalize=TRUE, intercept=TRUE)
    lambdas <- seq(l0, 0, by=-0.025*l0)
  } else if (missing(lambdas) & family=="multinomial"){
    lambdas <- seq(from=0.005,to=0.3,by=0.005)
  }

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

  if (repl<=0) stop ("repl has to be a positive number")

  output <- switch(family,
                   "gaussian"    = cv.gaussian.enetLTS(index, xx, yy, alphas, lambdas,
                                                       nfold, repl, ncores, plot),
                   "binomial"    = cv.binomial.enetLTS(index, xx, yy, alphas, lambdas,
                                                       nfold, repl, ncores, plot),
                   "multinomial" = cv.multinom.enetLTS(index, xx, yy, alphas, lambdas,
                                                       nfold, repl, ncores, plot)
                   )

}

