\name{cv.enetLTS}
\alias{cv.enetLTS}
\title{Cross-validation over the \eqn{(\alpha, \lambda)} grid}
\description{
Repeated \eqn{k}-fold cross-validation over a grid of \eqn{\alpha} and
\eqn{\lambda} values, as used inside \code{\link{enetLTS}}. Returns the selected
pair and, optionally, draws a heat map of the criterion.
}
\usage{
cv.enetLTS(index = NULL, family = c("gaussian", "binomial", "multinomial"),
           xx, yy, alphas, lambdas, nfold = 5, repl = 1, ncores = NA, plot = TRUE)
}
\arguments{
\item{index}{\code{NULL} (the default) or the subsets found by the C-steps: a list
  with one element per value of \code{alphas}, each a matrix with \eqn{h} rows
  (observation numbers) and one column per value of \code{lambdas}. If given, the
  cross-validation for a pair \eqn{(\alpha_i, \lambda_j)} uses only the
  observations \code{index[[i]][, j]}.}
\item{family}{\code{"gaussian"}, \code{"binomial"} or \code{"multinomial"}.}
\item{xx}{a numeric matrix of predictors.}
\item{yy}{the response: numeric (gaussian), 0/1 (binomial), or class labels or
  an indicator matrix (multinomial).}
\item{alphas}{the values of the mixing parameter \eqn{\alpha \in [0, 1]}
  (\eqn{\alpha = 1}: lasso, \eqn{\alpha = 0}: ridge).}
\item{lambdas}{the values of the penalty parameter. If missing, the default grid of
  \code{\link{enetLTS}} for the family is used.}
\item{nfold}{the number of folds. The default is 5.}
\item{repl}{the number of repetitions of the cross-validation; the criterion is
  averaged over the repetitions. The default is 1.}
\item{ncores}{the number of processor cores; \code{NA} (the default) uses all
  available cores.}
\item{plot}{logical; if \code{TRUE} (the default) a heat map of the criterion over
  the grid is drawn.}
}
\details{
The folds are random for the gaussian family and stratified by class for the
binomial and multinomial families. The criterion is the root mean squared
prediction error (gaussian) or the mean negative log-likelihood (binomial,
multinomial) of the test folds. If \code{index = NULL}, the observations are not
trimmed beforehand, so the largest 10\% of the test-fold losses are left out of
the mean (labels \code{RTMSPE}, \code{TMNLL}). For the multinomial family the
largest 10\% are always left out, within every class.
}
\value{
A list with components
\item{indexbest}{the subset of the selected pair (\code{NULL} if \code{index = NULL}).}
\item{evalCrit}{the matrix of the criterion, rows for \code{alphas} and columns
  for \code{lambdas}.}
\item{minevalCrit}{the smallest value of the criterion.}
\item{lambdaopt}{the selected \eqn{\lambda} (rounded to 4 decimals).}
\item{alphaopt}{the selected \eqn{\alpha} (rounded to 4 decimals).}
}
\note{
This function is called by \code{\link{enetLTS}}; it is exported so that the
tuning parameters can also be chosen for a given set of subsets.
}
\seealso{
\code{\link{enetLTS}}
}
\examples{
data("GaussianExample")
cv <- cv.enetLTS(family = "gaussian", xx = GaussianExample$x,
                 yy = GaussianExample$y, alphas = c(0.5, 1), ncores = 1)
cv[c("alphaopt", "lambdaopt")]
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{models}
\keyword{regression}
