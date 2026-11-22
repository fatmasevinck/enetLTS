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
data("GaussianExample")
data("BinomialExample")
data("MultinomialExample")
\donttest{
fit  <- enetLTS(GaussianExample$x, GaussianExample$y)
fitb <- enetLTS(BinomialExample$x, BinomialExample$y, family = "binomial")
fitm <- enetLTS(MultinomialExample$x, MultinomialExample$y, family = "multinomial")
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
