\name{plot.enetLTS}
\alias{plot.enetLTS}
\title{Plots of an \code{"enetLTS"} fit}
\description{
Draws the coefficient, residual and diagnostic plots of a fitted \code{"enetLTS"} object.
}
\usage{
\method{plot}{enetLTS}(x, method = c("coefficients", "resid", "diagnostic"),
     vers = c("reweighted", "raw"), ...)
}
\arguments{
\item{x}{a fitted object of class \code{"enetLTS"}.}
\item{method}{the plot to draw: \code{"coefficients"} (\code{\link{plotCoef.enetLTS}}),
  \code{"resid"} (\code{\link{plotResid.enetLTS}}) or \code{"diagnostic"}
  (\code{\link{plotDiagnostic.enetLTS}}). If \code{method} is not given, all
  three are drawn; in an interactive session R asks before each new page.}
\item{vers}{a character string: \code{"reweighted"} (the default) for the
  reweighted fit or \code{"raw"} for the raw fit.}
\item{\dots}{further arguments passed to the plot function, for example \code{zoom}
  for \code{method = "resid"} and \code{method = "diagnostic"}.}
}
\value{
The object \code{x}, invisibly. The plots are drawn as a side effect.
}
\seealso{
\code{\link{enetLTS}}, \code{\link{plotCoef.enetLTS}},
\code{\link{plotResid.enetLTS}}, \code{\link{plotDiagnostic.enetLTS}}
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
plot(fit, method = "resid")
plot(fitb, method = "diagnostic", vers = "raw")
plot(fitm, method = "coefficients")
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{hplot}
