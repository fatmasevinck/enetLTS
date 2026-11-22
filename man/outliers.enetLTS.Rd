\name{outliers.enetLTS}
\alias{outliers.enetLTS}
\alias{print.outliers.enetLTS}
\title{Observations flagged as outliers by an \code{"enetLTS"} fit}
\description{
Returns, for every training observation, whether it has weight zero, its loss, the
observed and the predicted response and a type label. The print method shows only
the flagged observations.
}
\usage{
outliers.enetLTS(object, vers = c("reweighted", "raw"))

\method{print}{outliers.enetLTS}(x, ...)
}
\arguments{
\item{object}{a fitted object of class \code{"enetLTS"}.}
\item{vers}{\code{"reweighted"} (the default; flagged means weight 0 in \code{wt})
  or \code{"raw"} (weight 0 in \code{raw.wt}).}
\item{x}{an object of class \code{"outliers.enetLTS"}.}
\item{\dots}{not used.}
}
\value{
A data frame of class \code{"outliers.enetLTS"} with one row per observation and
the columns \code{index}, \code{flagged}, \code{loss} (absolute residual for
\code{"gaussian"}, \eqn{-\log \hat{p}_{i,y_i}} otherwise), \code{observed},
\code{predicted} and \code{type}.

For \code{"gaussian"}, \code{type} is \code{"outlier"} for flagged observations.
For the classification families it is \code{"leverage point"} if a flagged
observation is predicted in its observed class, \code{"possible mislabel / outlier"}
if it is predicted in another class, and \code{"leverage point (screened)"} or
\code{"bad leverage point (screened)"} (predicted in the observed class or not) if
it was removed by the screening step of the multinomial family. Observations that
are not flagged are \code{"regular"}.
}
\seealso{
\code{\link{enetLTS}}, \code{\link{weights.enetLTS}}, \code{\link{plotResid.enetLTS}}
}
\examples{
data("GaussianExample")
data("BinomialExample")
data("MultinomialExample")
\donttest{
fit  <- enetLTS(GaussianExample$x, GaussianExample$y)
fitb <- enetLTS(BinomialExample$x, BinomialExample$y, family = "binomial")
fitm <- enetLTS(MultinomialExample$x, MultinomialExample$y, family = "multinomial")
out <- outliers.enetLTS(fit)
out
## compare with the contaminated observations
table(flagged = out$flagged, contaminated = out$index \%in\% GaussianExample$outliers)
outliers.enetLTS(fitb)
outliers.enetLTS(fitm)
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{robust}
\keyword{regression}
