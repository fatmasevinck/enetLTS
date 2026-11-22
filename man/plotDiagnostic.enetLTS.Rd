\name{plotDiagnostic.enetLTS}
\alias{plotDiagnostic.enetLTS}
\title{Diagnostic plot of an \code{"enetLTS"} fit}
\description{
Plots the response against the linear predictor (gaussian, binomial) or the first
two principal components of the linear predictors (multinomial).
}
\usage{
plotDiagnostic.enetLTS(object, vers = c("reweighted", "raw"), zoom = TRUE, ...)
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
For \code{family = "gaussian"} the plot shows \eqn{y} against
\eqn{X\hat{\beta}} with the line \eqn{y = X\hat{\beta}}. For
\code{family = "binomial"} it shows the 0/1 response against \eqn{X\hat{\beta}}
with the fitted logistic curve. For \code{family = "multinomial"} it shows the
scores of the first two principal components of the \eqn{K} linear predictors.

Observations are marked as flagged if they have weight 0 in \code{wt}
(\code{vers = "reweighted"}) or if they are not in the best subset \code{best}
of the raw fit (\code{vers = "raw"}). For the gaussian family regular
observations are blue crosses and flagged ones red triangles. For the binomial
and multinomial families the class is shown by color and shape, and flagged
observations are drawn filled, regular ones open.
}
\value{
The \code{"ggplot"} object, invisibly; the plot is drawn as a side effect.
}
\seealso{
\code{\link{enetLTS}}, \code{\link{plotResid.enetLTS}}, \code{\link{outliers.enetLTS}}
}
\examples{
data("GaussianExample")
data("BinomialExample")
data("MultinomialExample")
\donttest{
fit  <- enetLTS(GaussianExample$x, GaussianExample$y)
fitb <- enetLTS(BinomialExample$x, BinomialExample$y, family = "binomial")
fitm <- enetLTS(MultinomialExample$x, MultinomialExample$y, family = "multinomial")
plotDiagnostic.enetLTS(fit)
plotDiagnostic.enetLTS(fitb)
plotDiagnostic.enetLTS(fitm, zoom = FALSE)
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{hplot}
