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
data("GaussianExample")
data("BinomialExample")
data("MultinomialExample")
\donttest{
fit  <- enetLTS(GaussianExample$x, GaussianExample$y)
fitb <- enetLTS(BinomialExample$x, BinomialExample$y, family = "binomial")
fitm <- enetLTS(MultinomialExample$x, MultinomialExample$y, family = "multinomial")
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
