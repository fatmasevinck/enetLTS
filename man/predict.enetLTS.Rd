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
data("GaussianExample")
data("BinomialExample")
data("MultinomialExample")
\donttest{
fit  <- enetLTS(GaussianExample$x, GaussianExample$y)
fitb <- enetLTS(BinomialExample$x, BinomialExample$y, family = "binomial")
fitm <- enetLTS(MultinomialExample$x, MultinomialExample$y, family = "multinomial")
yhat <- predict(fit, GaussianExample$xtest, type = "response")$response
mean((GaussianExample$ytest - yhat)^2)
pb <- predict(fitb, BinomialExample$xtest, type = "class")$class
table(predicted = pb, observed = BinomialExample$ytest)
head(predict(fitm, MultinomialExample$xtest, type = "response")$response)
predict(fitm, MultinomialExample$xtest, vers = "raw", type = "class")$raw.class[1:10]
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{regression}
