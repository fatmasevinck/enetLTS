
#################################################################################################
# modified from lambda0 function in robustHD package

lambda00 <- function(x, y, normalize = TRUE, intercept = TRUE, const = 2,
                     prob = 0.95, tol = .Machine$double.eps^0.5,
                     eps = .Machine$double.eps, ...) {
   # initializations
   n <- length(y)
   x <- as.matrix(x)
   if(nrow(x) != n) stop(sprintf("'x' must have %d rows", n))
   normalize <- isTRUE(normalize)
   intercept <- isTRUE(intercept)

   x <- robStandardize(x, eps=eps, ...)
   centerX <- attr(x, "center")
   scaleX <- attr(x, "scale")
   # drop variables with too small a scale
   keep <- which(scaleX > eps)
   if(length(keep) == 0) stop("scale of all predictors is too small")
   x <- x[, keep, drop=FALSE]
   centerX <- centerX[keep]
   scaleX <- scaleX[keep]
   # compute largest lambda
   corY <- sapply(seq_len(ncol(x)), function(j) {
      tmp <- winsorize.default(x[, j], const=const, prob=prob, standardized=TRUE, tol=tol)
      rpbr(y,tmp)
   })
   max(abs(corY)) * 2 / nrow(x)
}

rpbr <- function(y,x){
   nj <- unlist(table(y))
   n <- sum(nj)
   M0 <- median(x[y==0])
   M1 <- median(x[y==1])
   return((M1-M0)/mad(x)*sqrt(prod(nj)/(n*(n-1))))
}

