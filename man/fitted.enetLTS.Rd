\name{fitted.enetLTS}
\alias{fitted.enetLTS}
\title{Fitted values of an \code{"enetLTS"} fit}
\description{
Fitted values of the reweighted and/or the raw fit for all observations of \code{xx}.
}
\usage{
\method{fitted}{enetLTS}(object, vers = c("reweighted", "raw", "both"),
       type = c("response", "class"), ...)
}
\arguments{
\item{object}{a fitted object of class \code{"enetLTS"}.}
\item{vers}{a character string: \code{"reweighted"} (the default) for the
  reweighted fit, \code{"raw"} for the raw fit, or \code{"both"} for both fits.}
\item{type}{\code{"response"} (the default) gives the fitted values for
  \code{family = "gaussian"} and the fitted probabilities for \code{"binomial"}
  and \code{"multinomial"}; \code{"class"} gives the class with the largest fitted
  probability (not available for \code{"gaussian"}).}
\item{\dots}{not used.}
}
\value{
A list with the component \code{fitted.values} (\code{vers = "reweighted"}),
\code{raw.fitted.values} (\code{vers = "raw"}) or both (\code{vers = "both"}).
Each is a one-column matrix (gaussian, binomial), a matrix with one column per
class (multinomial, \code{type = "response"}) or a vector of class labels
(\code{type = "class"}).
}
\seealso{
\code{\link{enetLTS}}, \code{\link{predict.enetLTS}}, \code{\link{residuals.enetLTS}}
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
head(fitted(fit)$fitted.values)
table(fitted(fitb, type = "class")$fitted.values, yb)
head(fitted(fitm, vers = "raw")$raw.fitted.values)
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{regression}
