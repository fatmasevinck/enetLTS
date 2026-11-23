\name{coef.enetLTS}
\alias{coef.enetLTS}
\title{Coefficients of an \code{"enetLTS"} fit}
\description{
Extracts the estimated coefficients of the reweighted or the raw fit.
}
\usage{
\method{coef}{enetLTS}(object, vers = c("reweighted", "raw"), zeros = TRUE, ...)
}
\arguments{
\item{object}{a fitted object of class \code{"enetLTS"}.}
\item{vers}{a character string: \code{"reweighted"} (the default) for the
  reweighted fit or \code{"raw"} for the raw fit.}
\item{zeros}{logical; if \code{TRUE} (the default) all coefficients are returned,
  if \code{FALSE} only the nonzero ones (for the multinomial family: the rows that
  are nonzero in at least one class).}
\item{\dots}{not used.}
}
\value{
For \code{family = "gaussian"} and \code{"binomial"}, a named numeric vector: the
intercept \code{"(Intercept)"} (if \code{intercept = TRUE}) followed by the
slopes, named by the column names of \code{xx} (\code{x1}, \code{x2}, \dots if
\code{xx} has no column names).
For \code{family = "multinomial"}, a matrix with one row per coefficient (named as
above) and one column per class.
}
\note{
In versions before 1.2.0, \code{zeros = FALSE} returned the indices of the nonzero
coefficients; these are now given by \code{\link{nonzeroCoef.enetLTS}}.
}
\seealso{
\code{\link{enetLTS}}, \code{\link{predict.enetLTS}},
\code{\link{nonzeroCoef.enetLTS}}, \code{\link{print.enetLTS}}
}
\examples{
## linear regression
set.seed(86)
n <- 100; p <- 25                              # observations and variables
beta <- c(rep(1, 6), rep(0, p - 6))            # six nonzero coefficients
x <- matrix(rnorm(n * p), n, p)
y <- drop(x \%*\% beta) + rnorm(n)
m <- 10                                        # contaminated observations
y[1:m] <- y[1:m] + 10                          # vertical outliers
x[1:m, ] <- x[1:m, ] + 10                      # bad leverage points

## logistic regression
set.seed(86)
n <- 100; p <- 25
beta <- c(rep(1, 6), rep(0, p - 6))
xb <- matrix(rnorm(n * p), n, p)
yb <- rbinom(n, 1, plogis(drop(xb \%*\% beta)))
mb <- 5
i0 <- which(yb == 0)[1:mb]                      # observations of class 0 ...
xb[i0, ] <- xb[i0, ] + 10                      # ... become bad leverage points

## multinomial regression
set.seed(86)
n <- 120; p <- 15; K <- 3                      # observations, variables, classes
Bm <- matrix(0, p, K)
Bm[1:2, 1] <- 3; Bm[3:4, 2] <- 3; Bm[5:6, 3] <- 3
xm <- matrix(rnorm(n * p), n, p)
pr <- exp(xm \%*\% Bm); pr <- pr / rowSums(pr)
ym <- apply(pr, 1, function(q) sample(K, 1, prob = q))
mm <- 10
xm[1:mm, 1:2] <- xm[1:mm, 1:2] + 5               # far out in the predictors of class 1 ...
ym[1:mm] <- 3                                   # ... but labeled as class 3

\donttest{
fit  <- enetLTS(x, y)
fitb <- enetLTS(xb, yb, family = "binomial")
fitm <- enetLTS(xm, ym, family = "multinomial")
coef(fit, zeros = FALSE)
coef(fit, vers = "raw", zeros = FALSE)
coef(fitb, zeros = FALSE)
coef(fitm, zeros = FALSE)
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{regression}
