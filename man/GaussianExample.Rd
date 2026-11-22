\name{GaussianExample}
\alias{GaussianExample}
\docType{data}
\title{Simulated high-dimensional linear regression data with outliers}
\description{
Simulated data for \code{family = "gaussian"} with \eqn{n = 100} observations and
\eqn{p = 200} predictors (\eqn{p > n}), following the simulation design of
Kurnaz, Hoffmann and Filzmoser (2018).
}
\usage{data("GaussianExample")}
\format{A list with components
\describe{
  \item{\code{x}}{training predictors, a \eqn{100 \times 200} matrix.}
  \item{\code{y}}{training response, a numeric vector of length 100.}
  \item{\code{outliers}}{indices of the contaminated training observations (1 to 10).}
  \item{\code{xtest}, \code{ytest}}{clean test data with 200 observations.}
  \item{\code{beta}}{the true slope vector (the true intercept is 1).}
}}
\details{
The predictors form three independent blocks: two informative blocks of 10 predictors each with
\eqn{\mathrm{cor}(x_j, x_k) = 0.9^{|j-k|}}, and 180 uninformative predictors with
\eqn{0.2^{|j-k|}}. The response is \eqn{y = 1 + x^\top\beta + \varepsilon} with
\eqn{\varepsilon \sim N(0,1)}, \eqn{\beta_j = 1} for the 20 informative predictors and 0 otherwise.
The first 10 observations are contaminated: their informative predictors are drawn from
\eqn{N(20, 1)} and their errors from \eqn{N(20\,\hat\sigma_y, 1)}, where \eqn{\hat\sigma_y} is the
standard deviation of the clean response, so they are both leverage points and vertical outliers.
Values are rounded to three decimals. The generating code is in the \code{data-raw} directory of
the package sources.
}
\references{
Kurnaz, F.S., Hoffmann, I. and Filzmoser, P. (2018). Robust and sparse
estimation methods for high-dimensional linear and logistic regression.
\emph{Chemometrics and Intelligent Laboratory Systems}, \bold{172}, 211--222.
\doi{10.1016/j.chemolab.2017.11.017}
}
\examples{
data("GaussianExample")
dim(GaussianExample$x)
GaussianExample$outliers
}
\keyword{datasets}
