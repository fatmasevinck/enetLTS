

cv.enetLTS <- function(index=NULL, family=c("gaussian","binomial","multinomial"), xx, yy, alphas, lambdas,
                       nfold=5, repl=1, ncores=NA, plot=TRUE)
{
  family <- match.arg(family)

  if (!is.null(index) && !is.list(index)) {
    stop("'index' has to be a list with one matrix of subsets (columns: lambdas) per alpha")
  }

  if (family=="multinomial" && is.null(dim(yy))){
    ## class labels -> indicator matrix, as in enetLTS()
    yy <- as.factor(yy)
    yy <- diag(nlevels(yy))[as.integer(yy), , drop=FALSE]
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
    lambdas <- multinomLambdas(xx, yy)
  }

  ncores <- checkNcores(ncores)
  startWorkers(ncores)
  on.exit(stopWorkers())

  if (repl<=0) stop ("repl has to be a positive number")

  switch(family,
                   "gaussian"    = cv.gaussian.enetLTS(index, xx, yy, alphas, lambdas,
                                                       nfold, repl, ncores, plot),
                   "binomial"    = cv.binomial.enetLTS(index, xx, yy, alphas, lambdas,
                                                       nfold, repl, ncores, plot),
                   "multinomial" = cv.multinom.enetLTS(index, xx, yy, alphas, lambdas,
                                                       nfold, repl, ncores, plot)
                   )

}

