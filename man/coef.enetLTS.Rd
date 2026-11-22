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
data("GaussianExample")
data("BinomialExample")
data("MultinomialExample")
\donttest{
fit  <- enetLTS(GaussianExample$x, GaussianExample$y)
fitb <- enetLTS(BinomialExample$x, BinomialExample$y, family = "binomial")
fitm <- enetLTS(MultinomialExample$x, MultinomialExample$y, family = "multinomial")
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
