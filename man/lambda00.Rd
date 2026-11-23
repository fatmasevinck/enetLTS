\name{lambda00}
\alias{lambda00}
\title{Robust upper end of the penalty grid for logistic regression}
\description{
Robust estimate of the smallest value of the penalty parameter that sets all
slopes of the (lasso-penalized) logistic regression to zero. \code{\link{enetLTS}}
uses it as the upper end of the default \code{lambdas} grid for
\code{family = "binomial"}.
}
\usage{
lambda00(x, y, normalize = TRUE, intercept = TRUE, const = 2, prob = 0.95,
         tol = .Machine$double.eps^0.5, eps = .Machine$double.eps, ...)
}
\arguments{
\item{x}{a numeric matrix of predictors.}
\item{y}{a numeric 0/1 response vector.}
\item{normalize}{logical; whether the winsorized predictors are normalized. The
  default is \code{TRUE}.}
\item{intercept}{logical; whether the model has an intercept. The default is
  \code{TRUE}.}
\item{const}{tuning constant of the univariate winsorization. The default is 2.}
\item{prob}{probability of the \eqn{\chi^2}{chi-squared} quantile used in the
  bivariate winsorization. The default is 0.95.}
\item{tol}{a small positive value used to detect singularity in the computation of
  the correlations for the bivariate winsorization.}
\item{eps}{a small positive value below which a robust scale estimate is treated as
  zero.}
\item{\dots}{further arguments passed to \code{\link[robustHD]{robStandardize}}.}
}
\details{
The procedure follows \code{\link[robustHD]{lambda0}} (Alfons et al., 2013),
which is used for \code{family = "gaussian"}: the predictors are robustly
standardized and the correlations between the response and each predictor are
computed on winsorized data. For the 0/1 response of logistic regression the
Pearson correlation is replaced by a robustified point-biserial correlation
(Kurnaz et al., 2018).
}
\value{
A single positive number.
}
\references{
Kurnaz, F.S., Hoffmann, I. and Filzmoser, P. (2018). Robust and sparse
estimation methods for high-dimensional linear and logistic regression.
\emph{Chemometrics and Intelligent Laboratory Systems}, \bold{172}, 211--222.
\doi{10.1016/j.chemolab.2017.11.017}

Alfons, A., Croux, C. and Gelper, S. (2013). Sparse least trimmed squares
regression for analyzing high-dimensional large data sets.
\emph{The Annals of Applied Statistics}, \bold{7}(1), 226--248.
}
\seealso{
\code{\link{enetLTS}}, \code{\link[robustHD]{lambda0}},
\code{\link[robustHD]{sparseLTS}}
}
\examples{
data("BinomialExample")
lambda00(BinomialExample$x, BinomialExample$y)
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{robust}
\keyword{regression}
