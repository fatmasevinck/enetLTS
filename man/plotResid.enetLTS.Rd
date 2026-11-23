\name{plotResid.enetLTS}
\alias{plotResid.enetLTS}
\title{Residual plots of an \code{"enetLTS"} fit}
\description{
Plots the residuals of the reweighted or the raw fit against the observation index
and, for the gaussian and binomial families, against the linear predictor.
}
\usage{
plotResid.enetLTS(object, vers = c("reweighted", "raw"), zoom = TRUE, ...)
}
\arguments{
\item{object}{a fitted object of class \code{"enetLTS"}.}
\item{vers}{a character string: \code{"reweighted"} (the default) for the
  reweighted fit or \code{"raw"} for the raw fit.}
\item{zoom}{logical; if \code{TRUE} (the default) and the flagged observations lie
  far from the regular ones, a second panel shows only the range of the regular
  observations (and of the flagged observations close to them), so that the fit and
  the cut-off lines remain visible.}
\item{\dots}{not used.}
}
\details{
For \code{family = "gaussian"} the residuals are standardized as in the
reweighting step (location and scale from the observations with raw weight 1 for
the reweighted fit, from the best subset for the raw fit), and the horizontal
lines at \eqn{\pm}\code{qnorm(1 - del)} are the cut-offs of the reweighting step.
For \code{family = "binomial"} and \code{"multinomial"} the residuals are the
negative log-likelihood contributions (see \code{\link{residuals.enetLTS}}).
The multinomial family has only the index plot.

Observations are marked as flagged if they have weight 0 in \code{wt}
(\code{vers = "reweighted"}) or if they are not in the best subset \code{best}
of the raw fit (\code{vers = "raw"}). For the gaussian family regular
observations are blue crosses and flagged ones red triangles. For the binomial
and multinomial families the class is shown by color and shape, and flagged
observations are drawn filled, regular ones open.
}
\value{
A list of the drawn \code{"ggplot"} objects (two for the gaussian and binomial
families, one for the multinomial family), invisibly.
}
\seealso{
\code{\link{enetLTS}}, \code{\link{residuals.enetLTS}},
\code{\link{plotDiagnostic.enetLTS}}, \code{\link{outliers.enetLTS}}
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
plotResid.enetLTS(fit)
plotResid.enetLTS(fit, vers = "raw", zoom = FALSE)
plotResid.enetLTS(fitb)
plotResid.enetLTS(fitm)
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{hplot}
