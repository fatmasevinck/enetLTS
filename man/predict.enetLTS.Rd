\name{predict.enetLTS}
\alias{predict.enetLTS}
\title{Predictions from an \code{"enetLTS"} fit}
\description{
Linear predictors, fitted values or probabilities, class labels, coefficients or the
indices of the nonzero coefficients from a fitted \code{"enetLTS"} object.
}
\usage{
\method{predict}{enetLTS}(object, newX, vers = c("reweighted", "raw"),
        type = c("link", "response", "coefficients", "nonzero", "class"), ...)
}
\arguments{
\item{object}{a fitted object of class \code{"enetLTS"}.}
\item{newX}{a numeric matrix (or an object that \code{as.matrix} converts to one)
  with the same predictors, in the same order, as \code{xx} in
  \code{\link{enetLTS}}. Not needed for \code{type = "coefficients"} and
  \code{type = "nonzero"}.}
\item{vers}{a character string: \code{"reweighted"} (the default) for the
  reweighted fit or \code{"raw"} for the raw fit.}
\item{type}{\code{"link"} (the default) gives the linear predictors;
  \code{"response"} the fitted values (gaussian) or the probabilities (binomial,
  multinomial); \code{"class"} the class with the largest probability (binomial,
  multinomial); \code{"coefficients"} the coefficients; \code{"nonzero"} the
  indices of the nonzero coefficients, as \code{\link{nonzeroCoef.enetLTS}}.}
\item{\dots}{not used.}
}
\value{
A list with one component. Its name is \code{link}, \code{response},
\code{class}, \code{coefficients} or \code{nonzeroCoef}, with the prefix
\code{raw.} for \code{vers = "raw"} (for example \code{raw.response}). For the
multinomial family, \code{link} and \code{response} are matrices with one column
per class.
}
\seealso{
\code{\link{enetLTS}}, \code{\link{coef.enetLTS}}, \code{\link{fitted.enetLTS}},
\code{\link{assess.enetLTS}}
}
\examples{
## linear regression
set.seed(86)
n <- 100; p <- 25                              # observations and variables
beta <- c(rep(1, 6), rep(0, p - 6))            # six nonzero coefficients
x <- matrix(rnorm(n * p), n, p)
y <- drop(x \%*\% beta) + rnorm(n)
xt <- matrix(rnorm(200 * p), 200, p)           # clean test data
yt <- drop(xt \%*\% beta) + rnorm(200)
m <- 10                                        # contaminated observations
y[1:m] <- y[1:m] + 10                          # vertical outliers
x[1:m, ] <- x[1:m, ] + 10                      # bad leverage points

## logistic regression
set.seed(86)
n <- 100; p <- 25
beta <- c(rep(1, 6), rep(0, p - 6))
xb <- matrix(rnorm(n * p), n, p)
yb <- rbinom(n, 1, plogis(drop(xb \%*\% beta)))
xbt <- matrix(rnorm(200 * p), 200, p)          # clean test data
ybt <- rbinom(200, 1, plogis(drop(xbt \%*\% beta)))
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
xmt <- matrix(rnorm(200 * p), 200, p)         # clean test data
pr <- exp(xmt \%*\% Bm); pr <- pr / rowSums(pr)
ymt <- apply(pr, 1, function(q) sample(K, 1, prob = q))
mm <- 10
xm[1:mm, 1:2] <- xm[1:mm, 1:2] + 5               # far out in the predictors of class 1 ...
ym[1:mm] <- 3                                   # ... but labeled as class 3

\donttest{
fit  <- enetLTS(x, y)
fitb <- enetLTS(xb, yb, family = "binomial")
fitm <- enetLTS(xm, ym, family = "multinomial")
yhat <- predict(fit, xt, type = "response")$response
mean((yt - yhat)^2)
pb <- predict(fitb, xbt, type = "class")$class
table(predicted = pb, observed = ybt)
head(predict(fitm, xmt, type = "response")$response)
predict(fitm, xmt, vers = "raw", type = "class")$raw.class[1:10]
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{regression}
