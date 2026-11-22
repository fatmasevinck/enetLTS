\name{residuals.enetLTS}
\alias{residuals.enetLTS}
\title{Residuals of an \code{"enetLTS"} fit}
\description{
Residuals of the reweighted and/or the raw fit for all observations of \code{xx}.
}
\usage{
\method{residuals}{enetLTS}(object, vers = c("reweighted", "raw", "both"), ...)
}
\arguments{
\item{object}{a fitted object of class \code{"enetLTS"}.}
\item{vers}{a character string: \code{"reweighted"} (the default) for the
  reweighted fit, \code{"raw"} for the raw fit, or \code{"both"} for both fits.}
\item{\dots}{not used.}
}
\details{
For \code{family = "gaussian"} the residuals are \eqn{y_i - \hat{y}_i}.
For \code{family = "binomial"} and \code{"multinomial"} they are the negative
log-likelihood contributions of the observations,
\eqn{-y_i\eta_i + \log(1 + e^{\eta_i})} and \eqn{-\log \hat{p}_{i,y_i}},
respectively (\eqn{\eta_i} is the linear predictor and \eqn{\hat{p}_{i,y_i}} the
fitted probability of the observed class).
}
\value{
A list with the component \code{residuals} (\code{vers = "reweighted"}),
\code{raw.residuals} (\code{vers = "raw"}) or both (\code{vers = "both"}); each is
a numeric vector.
}
\seealso{
\code{\link{enetLTS}}, \code{\link{fitted.enetLTS}}, \code{\link{plotResid.enetLTS}},
\code{\link{outliers.enetLTS}}
}
\examples{
data("GaussianExample")
data("BinomialExample")
data("MultinomialExample")
\donttest{
fit  <- enetLTS(GaussianExample$x, GaussianExample$y)
fitb <- enetLTS(BinomialExample$x, BinomialExample$y, family = "binomial")
fitm <- enetLTS(MultinomialExample$x, MultinomialExample$y, family = "multinomial")
head(residuals(fit)$residuals)
str(residuals(fitb, vers = "both"))
summary(residuals(fitm)$residuals)
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{regression}
