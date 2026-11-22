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
data("GaussianExample")
data("BinomialExample")
data("MultinomialExample")
\donttest{
fit  <- enetLTS(GaussianExample$x, GaussianExample$y)
fitb <- enetLTS(BinomialExample$x, BinomialExample$y, family = "binomial")
fitm <- enetLTS(MultinomialExample$x, MultinomialExample$y, family = "multinomial")
head(fitted(fit)$fitted.values)
table(fitted(fitb, type = "class")$fitted.values, BinomialExample$y)
head(fitted(fitm, vers = "raw")$raw.fitted.values)
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{regression}
