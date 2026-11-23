\name{weights.enetLTS}
\alias{weights.enetLTS}
\title{Binary outlier weights of an \code{"enetLTS"} fit}
\description{
Extracts the 0/1 weights that mark the outliers of the reweighted and/or the raw fit.
}
\usage{
\method{weights}{enetLTS}(object, vers = c("reweighted", "raw", "both"), index = FALSE, ...)
}
\arguments{
\item{object}{a fitted object of class \code{"enetLTS"}.}
\item{vers}{a character string: \code{"reweighted"} (the default) for the
  reweighted fit, \code{"raw"} for the raw fit, or \code{"both"} for both fits.}
\item{index}{logical; if \code{TRUE}, the weights are named by the observation
  numbers. The default is \code{FALSE}.}
\item{\dots}{not used.}
}
\details{
An observation has weight 1 if it is regular and 0 if it is flagged as an outlier.
The raw weights \code{raw.wt} come from the raw fit and define the observations
used for the reweighted fit; the weights \code{wt} come from the reweighted fit.
With \code{del} the argument of \code{\link{enetLTS}}, the rules are:
\describe{
\item{gaussian}{standardized residuals (location and scale estimated from the
  observations with raw weight 1, with a consistency factor) larger than
  \code{qnorm(1 - del)} in absolute value get weight 0.}
\item{binomial}{Pearson residuals larger than \code{qnorm(1 - del)} in absolute
  value get weight 0.}
\item{multinomial}{with \code{crossfit = TRUE} (the default) an observation gets
  weight 0 if its cross-fitted probability of the observed class is below a
  calibrated cut-off (see \code{crossfit} in \code{\link{enetLTS}}); with
  \code{crossfit = FALSE} it gets weight 0 if it is outlying in the space of the
  linear predictors (class-wise robust distances) or if the fitted probability of
  its observed class is below 1/6. Observations removed by \code{screen = TRUE}
  have weight 0 in both vectors.}
}
}
\value{
A list with the component \code{wt} (\code{vers = "reweighted"}), \code{raw.wt}
(\code{vers = "raw"}) or both (\code{vers = "both"}); each is an integer vector of
length \eqn{n}.
}
\seealso{
\code{\link{enetLTS}}, \code{\link{outliers.enetLTS}}
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
which(weights(fit)$wt == 0)
which(weights(fitb)$wt == 0)
w <- weights(fitm, vers = "both")
table(raw = w$raw.wt, reweighted = w$wt)
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{regression}
