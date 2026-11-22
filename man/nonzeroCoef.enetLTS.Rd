\name{nonzeroCoef.enetLTS}
\alias{nonzeroCoef.enetLTS}
\title{Indices of the nonzero coefficients of an \code{"enetLTS"} fit}
\description{
Returns the positions of the nonzero coefficients of the reweighted or the raw fit.
}
\usage{
nonzeroCoef.enetLTS(object, vers = c("reweighted", "raw"))
}
\arguments{
\item{object}{a fitted object of class \code{"enetLTS"}.}
\item{vers}{a character string: \code{"reweighted"} (the default) for the
  reweighted fit or \code{"raw"} for the raw fit.}
}
\value{
For \code{family = "gaussian"} and \code{"binomial"}, an integer vector with the
positions of the nonzero entries of \code{coef(object, vers = vers)}. With
\code{intercept = TRUE} position 1 is the intercept, so the slope of predictor
\eqn{j} has position \eqn{j + 1}.
For \code{family = "multinomial"}, a list with one such vector per class (a data
frame if all classes have the same number of nonzero coefficients).
}
\seealso{
\code{\link{enetLTS}}, \code{\link{coef.enetLTS}}, \code{\link{predict.enetLTS}}
}
\examples{
data("GaussianExample")
data("BinomialExample")
data("MultinomialExample")
\donttest{
fit  <- enetLTS(GaussianExample$x, GaussianExample$y)
fitb <- enetLTS(BinomialExample$x, BinomialExample$y, family = "binomial")
fitm <- enetLTS(MultinomialExample$x, MultinomialExample$y, family = "multinomial")
nonzeroCoef.enetLTS(fit)
nonzeroCoef.enetLTS(fitb, vers = "raw")
nonzeroCoef.enetLTS(fitm)
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{regression}
