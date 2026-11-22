\name{plot.enetLTS}
\alias{plot.enetLTS}
\title{Plots of an \code{"enetLTS"} fit}
\description{
Draws the coefficient, residual and diagnostic plots of a fitted \code{"enetLTS"} object.
}
\usage{
\method{plot}{enetLTS}(x, method = c("coefficients", "resid", "diagnostic"),
     vers = c("reweighted", "raw"), ...)
}
\arguments{
\item{x}{a fitted object of class \code{"enetLTS"}.}
\item{method}{the plot to draw: \code{"coefficients"} (\code{\link{plotCoef.enetLTS}}),
  \code{"resid"} (\code{\link{plotResid.enetLTS}}) or \code{"diagnostic"}
  (\code{\link{plotDiagnostic.enetLTS}}). If \code{method} is not given, all
  three are drawn; in an interactive session R asks before each new page.}
\item{vers}{a character string: \code{"reweighted"} (the default) for the
  reweighted fit or \code{"raw"} for the raw fit.}
\item{\dots}{further arguments passed to the plot function, for example \code{zoom}
  for \code{method = "resid"} and \code{method = "diagnostic"}.}
}
\value{
The object \code{x}, invisibly. The plots are drawn as a side effect.
}
\seealso{
\code{\link{enetLTS}}, \code{\link{plotCoef.enetLTS}},
\code{\link{plotResid.enetLTS}}, \code{\link{plotDiagnostic.enetLTS}}
}
\examples{
data("GaussianExample")
data("BinomialExample")
data("MultinomialExample")
\donttest{
fit  <- enetLTS(GaussianExample$x, GaussianExample$y)
fitb <- enetLTS(BinomialExample$x, BinomialExample$y, family = "binomial")
fitm <- enetLTS(MultinomialExample$x, MultinomialExample$y, family = "multinomial")
plot(fit, method = "resid")
plot(fitb, method = "diagnostic", vers = "raw")
plot(fitm, method = "coefficients")
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{hplot}
