\name{plotCoef.enetLTS}
\alias{plotCoef.enetLTS}
\title{Coefficient plot of an \code{"enetLTS"} fit}
\description{
Needle plot of the estimated slopes against the predictor index.
}
\usage{
plotCoef.enetLTS(object, vers = c("reweighted", "raw"), colors = NULL, ...)
}
\arguments{
\item{object}{a fitted object of class \code{"enetLTS"}.}
\item{vers}{a character string: \code{"reweighted"} (the default) for the
  reweighted fit or \code{"raw"} for the raw fit.}
\item{colors}{not used; kept for compatibility with earlier versions.}
\item{\dots}{not used.}
}
\details{
Every slope is drawn as a vertical segment from zero, and the nonzero slopes are
marked by a point. The intercept (one per class for the multinomial family) is
given in the subtitle, together with the number of nonzero slopes. For the
multinomial family there is one panel per class.
}
\value{
The \code{"ggplot"} object, invisibly; the plot is drawn as a side effect.
}
\seealso{
\code{\link{enetLTS}}, \code{\link{coef.enetLTS}}, \code{\link{plot.enetLTS}}
}
\examples{
data("GaussianExample")
data("BinomialExample")
data("MultinomialExample")
\donttest{
fit  <- enetLTS(GaussianExample$x, GaussianExample$y)
fitb <- enetLTS(BinomialExample$x, BinomialExample$y, family = "binomial")
fitm <- enetLTS(MultinomialExample$x, MultinomialExample$y, family = "multinomial")
plotCoef.enetLTS(fit)
plotCoef.enetLTS(fitb, vers = "raw")
plotCoef.enetLTS(fitm)
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{hplot}
