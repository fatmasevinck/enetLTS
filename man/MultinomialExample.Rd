\name{MultinomialExample}
\alias{MultinomialExample}
\docType{data}
\title{Simulated high-dimensional multinomial regression data with outliers}
\description{
Simulated data for \code{family = "multinomial"} with \eqn{n = 150} observations,
\eqn{p = 200} predictors (\eqn{p > n}) and \eqn{K = 3} classes.
}
\usage{data("MultinomialExample")}
\format{A list with components
\describe{
  \item{\code{x}}{training predictors, a \eqn{150 \times 200} matrix.}
  \item{\code{y}}{training response, a factor with levels 1, 2, 3.}
  \item{\code{outliers}}{indices of the contaminated training observations (1 to 15).}
  \item{\code{xtest}, \code{ytest}}{clean test data with 300 observations.}
  \item{\code{beta}}{the true \eqn{200 \times 3} matrix of slopes (the true intercepts are 0).}
}}
\details{
The predictors are independent standard normal. Each class is driven by its own predictors:
\eqn{\beta_{1,1} = \beta_{2,1} = 3} (class 1), \eqn{\beta_{3,2} = \beta_{4,2} = 3} (class 2),
\eqn{\beta_{5,3} = \beta_{6,3} = 3} (class 3); the class probabilities follow the multinomial
logistic model. The first 15 observations are contaminated: predictors 1 and 2 are shifted by 5
(far out in the direction of class 1) and the label is set to 3, so they are bad leverage points.
Values are rounded to three decimals. The generating code is in the \code{data-raw} directory of
the package sources.
}
\examples{
data("MultinomialExample")
table(MultinomialExample$y)
}
\keyword{datasets}
