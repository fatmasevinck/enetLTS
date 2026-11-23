\name{print.enetLTS}
\alias{print.enetLTS}
\title{Print an \code{"enetLTS"} fit}
\description{
Prints a short summary of a fitted \code{"enetLTS"} object.
}
\usage{
\method{print}{enetLTS}(x, vers = c("reweighted", "raw"),
      digits = max(3, getOption("digits") - 3), maxprint = 20, ...)
}
\arguments{
\item{x}{a fitted object of class \code{"enetLTS"}.}
\item{vers}{a character string: \code{"reweighted"} (the default) for the
  reweighted fit or \code{"raw"} for the raw fit.}
\item{digits}{the number of digits for the coefficients.}
\item{maxprint}{the largest number of nonzero coefficients (for the multinomial
  family: rows) that are printed; all are available with \code{\link{coef.enetLTS}}.}
\item{\dots}{not used.}
}
\details{
The family and the call are printed, followed by the nonzero coefficients with the
variable names (\code{x1}, \code{x2}, \dots if \code{xx} has no column names), the
number of nonzero slopes, the selected \code{alpha}, \code{lambda} and
\code{lambdaw}, and the number of observations flagged as outliers. For the
multinomial family the rows of the coefficient matrix that are nonzero in at least
one class are printed.
}
\value{
The object \code{x}, invisibly.
}
\seealso{
\code{\link{enetLTS}}, \code{\link{coef.enetLTS}}, \code{\link{outliers.enetLTS}}
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
fit
print(fit, vers = "raw")
fitb
fitm
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{models}
\keyword{regression}
