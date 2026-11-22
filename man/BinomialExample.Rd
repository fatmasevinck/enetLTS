\name{BinomialExample}
\alias{BinomialExample}
\docType{data}
\title{Simulated high-dimensional logistic regression data with outliers}
\description{
Simulated data for \code{family = "binomial"} with \eqn{n = 100} observations and
\eqn{p = 200} predictors (\eqn{p > n}), following the simulation design of
Kurnaz, Hoffmann and Filzmoser (2018).
}
\usage{data("BinomialExample")}
\format{A list with components
\describe{
  \item{\code{x}}{training predictors, a \eqn{100 \times 200} matrix.}
  \item{\code{y}}{training response, a 0/1 vector of length 100.}
  \item{\code{outliers}}{indices of the contaminated training observations (1, 3, 7, 13).}
  \item{\code{xtest}, \code{ytest}}{clean test data with 200 observations.}
  \item{\code{beta}}{the true slope vector (the true intercept is 1).}
}}
\details{
The first 20 predictors are informative (\eqn{\beta_j = 1}) and correlated with
\eqn{\mathrm{cor}(x_j, x_k) = 0.9^{|j-k|}}; the other 180 are uninformative and correlated with
\eqn{0.5^{|j-k|}}. The response is \eqn{y = 1} if \eqn{1 + x^\top\beta + \varepsilon > 0} with
\eqn{\varepsilon \sim N(0,1)}, and \eqn{y = 0} otherwise. The first 10\% of the observations of
class 0 (4 observations) are contaminated: their informative predictors are drawn from \eqn{N(20, 1)} while their
label stays 0, so they are bad leverage points. Values are rounded to three decimals. The
generating code is in the \code{data-raw} directory of the package sources.
}
\references{
Kurnaz, F.S., Hoffmann, I. and Filzmoser, P. (2018). Robust and sparse
estimation methods for high-dimensional linear and logistic regression.
\emph{Chemometrics and Intelligent Laboratory Systems}, \bold{172}, 211--222.
\doi{10.1016/j.chemolab.2017.11.017}
}
\examples{
data("BinomialExample")
table(BinomialExample$y)
BinomialExample$outliers
}
\keyword{datasets}
