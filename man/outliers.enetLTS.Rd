\name{outliers.enetLTS}
\alias{outliers.enetLTS}
\alias{print.outliers.enetLTS}
\title{Observations flagged as outliers by an \code{"enetLTS"} fit}
\description{
Returns, for every training observation, whether it has weight zero, its loss, the
observed and the predicted response and a type label. The print method shows only
the flagged observations.
}
\usage{
outliers.enetLTS(object, vers = c("reweighted", "raw"))

\method{print}{outliers.enetLTS}(x, ...)
}
\arguments{
\item{object}{a fitted object of class \code{"enetLTS"}.}
\item{vers}{\code{"reweighted"} (the default; flagged means weight 0 in \code{wt})
  or \code{"raw"} (weight 0 in \code{raw.wt}).}
\item{x}{an object of class \code{"outliers.enetLTS"}.}
\item{\dots}{not used.}
}
\value{
A data frame of class \code{"outliers.enetLTS"} with one row per observation and
the columns \code{index}, \code{flagged}, \code{loss} (absolute residual for
\code{"gaussian"}, \eqn{-\log \hat{p}_{i,y_i}} otherwise), \code{observed},
\code{predicted} and \code{type}.

For \code{"gaussian"}, \code{type} is \code{"outlier"} for flagged observations.
For the classification families it is \code{"leverage point"} if a flagged
observation is predicted in its observed class, \code{"possible mislabel / outlier"}
if it is predicted in another class, and \code{"leverage point (screened)"} or
\code{"bad leverage point (screened)"} (predicted in the observed class or not) if
it was removed by the screening step of the multinomial family. Observations that
are not flagged are \code{"regular"}.
}
\seealso{
\code{\link{enetLTS}}, \code{\link{weights.enetLTS}}, \code{\link{plotResid.enetLTS}}
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
out <- outliers.enetLTS(fit)
out
## compare with the contaminated observations 1, ..., m
table(flagged = out$flagged, contaminated = out$index \%in\% 1:m)
outliers.enetLTS(fitb)
outliers.enetLTS(fitm)
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{robust}
\keyword{regression}
