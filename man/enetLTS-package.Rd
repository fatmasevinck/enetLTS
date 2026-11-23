\name{enetLTS}
\alias{enetLTS}
\alias{enetLTS-package}
\title{Robust and Sparse Methods for High Dimensional Linear, Binary and
Multinomial Regression}
\description{
Computes the robust elastic net estimator enetLTS for linear regression
(\code{family = "gaussian"}), binary logistic regression (\code{"binomial"}) and
multinomial logistic regression (\code{"multinomial"}). The estimator searches for
an outlier-free subset of \eqn{h} observations on which the classical elastic net
is fitted, and a reweighting step improves its efficiency.
}
\usage{
enetLTS(xx, yy, family = c("gaussian", "binomial", "multinomial"),
        alphas = seq(0, 1, length = 41), lambdas = NULL, lambdaw = NULL,
        intercept = TRUE, scal = TRUE, hsize = 0.75, nsamp = c(500, 10),
        nCsteps = 20, nfold = 5, repl = 1, ncores = 1, tol = -1e6,
        seed = NULL, del = 0.0125, crit.plot = FALSE, typegrouped = FALSE,
        type.response = c("link", "response", "class"),
        hgrid = c(0.65, 0.75, 0.85, 0.95), screen = TRUE, crossfit = NULL)
}
\arguments{
\item{xx}{a numeric matrix of predictors, \eqn{n \times p}; \eqn{p > n} is allowed.}
\item{yy}{the response. Numeric for \code{family = "gaussian"}; a numeric 0/1
  vector (numeric, or a factor with levels 0 and 1) for
  \code{"binomial"}; for \code{"multinomial"} a factor (or a vector of
  class labels) with \eqn{K \ge 3} classes, or an \eqn{n \times K} indicator matrix.}
\item{family}{\code{"gaussian"} (the default), \code{"binomial"} or
  \code{"multinomial"}.}
\item{alphas}{the grid of the elastic net mixing parameter
  \eqn{\alpha \in [0, 1]} (\eqn{\alpha = 1}: lasso, \eqn{\alpha = 0}: ridge).
  The default is 41 equally spaced values.}
\item{lambdas}{the grid of the penalty parameter \eqn{\lambda}. By default the
  grid runs from \eqn{\lambda_0} down to 0 in steps of \eqn{0.025\lambda_0},
  where \eqn{\lambda_0} is a robust estimate of the smallest penalty that sets all
  slopes to zero: \code{\link[robustHD]{lambda0}} for \code{"gaussian"} and
  \code{\link{lambda00}} for \code{"binomial"}. For \code{"multinomial"} the
  default is 20 values, equally spaced on the log scale, from the
  \eqn{\lambda_{max}} of \pkg{glmnet} (\eqn{\alpha = 1}, robustly standardized
  data) down to 1\% of it.}
\item{lambdaw}{the penalty of the reweighting step. By default it is chosen by
  \code{\link[glmnet]{cv.glmnet}} on the observations with raw weight 1.}
\item{intercept}{logical; whether the model has an intercept (for
  \code{"multinomial"}: class intercepts). The default is \code{TRUE}.}
\item{scal}{logical; if \code{TRUE} (the default), the predictors are centered and
  scaled by the means and standard deviations of the current subset in every
  C-step, and for \code{"gaussian"} the response is centered as well.}
\item{hsize}{the proportion of observations in the subset, \eqn{h =
  \lfloor (n+1) \cdot}\code{hsize}\eqn{\rfloor}. The default is 0.75.
  \emph{Experimental:} with \code{hsize = "optimal"} the model is fitted for every
  value in \code{hgrid}, and the value with the smallest cross-fitted loss on all
  observations (ignoring the worst 25\%) is chosen, taking the largest value within
  10\% of the best criterion. This takes about \code{length(hgrid)} times as long
  as a single fit. The result then contains \code{hsize.search} with the criterion
  of every candidate, and a warning is given if the chosen value is below 0.70.}
\item{nsamp}{a vector of two integers. Two C-steps are carried out on each of
  \code{nsamp[1]} elemental subsets (default 500), and the \code{nsamp[2]} best of
  them (default 10) are iterated further.}
\item{nCsteps}{the largest number of C-steps in one iteration; the C-steps stop
  earlier when the criterion no longer improves. The default is 20.}
\item{nfold}{the number of folds of the cross-validation. The default is 5.}
\item{repl}{the number of repetitions of the cross-validation. The default is 1.}
\item{ncores}{the number of processor cores (package \pkg{parallel}); \code{NA}
  uses all available cores. The default is 1.}
\item{tol}{not used by the current algorithm; kept for compatibility with earlier
  versions.}
\item{seed}{an optional seed for the random number generator, set before the
  elemental subsets are drawn (and, for \code{"multinomial"}, before the screening
  step). The default \code{NULL} does not set a seed.}
\item{del}{the tail probability of the reweighting step. For \code{"gaussian"} and
  \code{"binomial"} an observation gets weight 0 if its standardized or Pearson
  residual is larger than \code{qnorm(1 - del)} in absolute value; with the
  default 0.0125 this flags about 2.5\% of the observations under the model. For
  \code{"multinomial"} with \code{crossfit = TRUE} the cut-off is calibrated at
  \code{2 * del} (see \code{crossfit}).}
\item{crit.plot}{logical; if \code{TRUE} a heat map of the cross-validation
  criterion over the \eqn{(\alpha, \lambda)} grid is drawn. The default is
  \code{FALSE}.}
\item{typegrouped}{logical; \code{"multinomial"} only. If \code{TRUE}, the grouped
  lasso penalty of \pkg{glmnet} (\code{type.multinomial = "grouped"}) is used in the
  fits on the best subset and in the reweighting step. The default is \code{FALSE}.}
\item{type.response}{the scale of \code{fitted.values} and
  \code{raw.fitted.values} in the result: \code{"link"} (the default) for the linear
  predictors, \code{"response"} for the fitted values (gaussian) or probabilities,
  \code{"class"} for class labels (not for \code{"gaussian"}).}
\item{hgrid}{the candidate values of \code{hsize} (between 0.5 and 1), used only if
  \code{hsize = "optimal"}.}
\item{screen}{logical; \code{"multinomial"} only. If \code{TRUE} (the default),
  observations that are outlying in the predictors (leverage points) are found
  beforehand without using the class labels (spatial-sign principal components and
  robust distances of their scores) and are left out of the fit; they get weight
  zero in \code{raw.wt} and \code{wt} and are listed in \code{screened}. This
  protects against contamination that falls mostly into one class, for which the
  class-wise robust distances of the reweighting step are unreliable.}
\item{crossfit}{logical or \code{NULL} (the default); \code{"gaussian"} and
  \code{"multinomial"} only. \code{NULL} means \code{TRUE} for
  \code{"multinomial"} with \code{scal = TRUE} and \code{FALSE} otherwise. If
  \code{TRUE}, the observations that are kept are chosen by \emph{cross-fitted}
  residuals: an observation is never judged by a model that has seen it. The
  ordinary trimming keeps the observations that the current fit describes best;
  with many predictors this partly selects observations that noise variables happen
  to explain, which raises the false positive rate and the prediction error even on
  clean data. After the trimming, \eqn{\alpha} and \eqn{\lambda} are chosen again on
  the reweighted observations (common folds), within the range of \code{alphas}.
  For \code{"multinomial"} the cross-fitted loss is \eqn{-\log \hat{p}(y_i)}, and the
  reweighting step keeps an observation if its cross-fitted probability of the
  observed class, from a model fitted on all screened observations, is at least a
  cut-off \eqn{c}. The cut-off is calibrated: \eqn{c} is the largest value for
  which, if these probabilities were correct, the expected share of observations
  with \eqn{\hat{p}(y_i) < c} is at most \code{2 * del}. An earlier rule, a robust
  median/MAD bound on the cross-fitted loss within every class, is used with
  \code{options(enetLTS.wrule = "mad")}. The final model contains class intercepts.
  Requires \code{scal = TRUE}. For \code{"binomial"} it is ignored with a warning.}
}
\details{
The estimator applies the classical elastic net (\pkg{glmnet}) to subsets of the
data only. For one pair \eqn{(\alpha, \lambda)} the search starts from
\code{nsamp[1]} elemental subsets and improves them by C-steps; the C-steps keep the
\eqn{h} observations with the smallest squared residuals (\code{"gaussian"}) or
negative log-likelihood contributions (\code{"binomial"}, within each class). For
\code{"multinomial"} they keep, within each class \eqn{k}, the
\eqn{\lfloor n_k \cdot}\code{hsize}\eqn{+ 1\rfloor} observations with the smallest
robust (MCD) distance in the space of the linear predictors.
For \code{"binomial"} the start subsets are ranked by the Bianco-Yohai criterion
(Croux and Haesbroeck, 2003), as in Kurnaz et al. (2018). The other grid points
are computed with warm starts. The pair \eqn{(\alpha, \lambda)} is chosen by
\code{repl} times repeated \code{nfold}-fold cross-validation on the best subsets
(see \code{\link{cv.enetLTS}}), and the raw fit is the elastic net on the best
subset of the chosen pair.

In the reweighting step the observations with large residuals of the raw fit get
weight 0 (see \code{del} and \code{\link{weights.enetLTS}}), and the elastic net
with penalty \code{lambdaw} is fitted to the observations with weight 1; this is
the reweighted fit. The weights \code{wt} are then computed from the reweighted
fit with the same rule. With \code{crossfit = TRUE} the raw weights come from the
cross-fitted trimming instead, and for \code{"multinomial"} \code{wt} then equals
\code{raw.wt}.
}
\value{
An object of class \code{"enetLTS"}, a list with components
\item{objective}{the penalized criterion of the raw fit on the final subset
  \eqn{H} (the \eqn{h} retained observations), on the scale on which the C-steps
  fit and penalize: \eqn{\frac{1}{h}\sum_{i\in H}\ell_i + \lambda P_\alpha(s_H \circ \beta)},
  where \eqn{\ell_i} is the squared residual divided by two (\code{"gaussian"}) or
  the negative log-likelihood contribution (\code{"binomial"}, \code{"multinomial"}), \eqn{\beta} the raw
  coefficients without the intercept, \eqn{s_H} the standard deviations of the
  columns of \code{xx} on \eqn{H} (one if \code{scal = FALSE} and
  \code{family = "multinomial"}), and \eqn{(\alpha, \lambda)} those of the raw fit.
  It is reported only and is not used by the algorithm. The search for
  \code{"binomial"} ranks subsets by the Bianco-Yohai criterion without penalty, so
  for that family the value is not what the C-steps optimize.}
\item{raw.rmse, rmse}{\code{"gaussian"} only: root mean squared residual of the
  raw and the reweighted fit on all observations.}
\item{raw.mae, mae}{mean absolute residual of the raw and the reweighted fit on all
  observations; for \code{"binomial"} and \code{"multinomial"} the mean absolute
  difference between the 0/1 indicators and the fitted probabilities.}
\item{best}{the indices of the best subset of \eqn{h} observations, used for the
  raw fit.}
\item{raw.wt}{integer 0/1 weights from the raw fit; the observations with weight 1
  are used for the reweighted fit.}
\item{wt}{integer 0/1 weights from the reweighted fit (0 = flagged as outlier).}
\item{raw.coefficients, coefficients}{the coefficients of the raw and the
  reweighted fit, intercept first (a matrix with one column per class for
  \code{"multinomial"}). Use \code{\link{coef.enetLTS}} to get them with names.}
\item{raw.fitted.values, fitted.values}{the fitted values of the raw and the
  reweighted fit, on the scale given by \code{type.response}.}
\item{raw.residuals, residuals}{the residuals (\code{"gaussian"}) or the negative
  log-likelihood contributions (\code{"binomial"}, \code{"multinomial"}) of the
  raw and the reweighted fit.}
\item{alpha}{the selected mixing parameter of the reweighted fit.}
\item{alpha.raw}{\code{"gaussian"} and \code{"multinomial"}: the mixing parameter
  of the raw fit. It equals \code{alpha} unless \code{crossfit} selects
  \eqn{\alpha} again for the reweighted fit.}
\item{lambda}{the selected penalty of the raw fit.}
\item{lambdaw}{the penalty of the reweighted fit.}
\item{num.nonzerocoef}{the number of nonzero coefficients of the reweighted fit.}
\item{n, p}{the number of observations and of predictors.}
\item{h}{the size of the subsets.}
\item{classnames, classize}{\code{"binomial"} and \code{"multinomial"}: the class
  labels and the class sizes.}
\item{screened}{\code{"multinomial"} only: the indices of the observations removed
  by the screening step (\code{screen = TRUE}).}
\item{hsize.search}{only if \code{hsize = "optimal"}: the criterion of every
  candidate value.}
\item{inputs}{a list with the data and the settings used by the fit.}
\item{call}{the matched call.}
}
\references{
Kurnaz, F.S. and Filzmoser, P. (2023). enetLTS: Robust and sparse methods
for high dimensional linear, binary, and multinomial regression.
\emph{Journal of Open Source Software}, \bold{8}(82), 4773.
\doi{10.21105/joss.04773}

Kurnaz, F.S. and Filzmoser, P. (2023). Robust and sparse multinomial
regression in high dimensions. \emph{Data Mining and Knowledge Discovery},
\bold{37}(4), 1609--1629. \doi{10.1007/s10618-023-00936-6}

Kurnaz, F.S., Hoffmann, I. and Filzmoser, P. (2018). Robust and sparse
estimation methods for high-dimensional linear and logistic regression.
\emph{Chemometrics and Intelligent Laboratory Systems}, \bold{172}, 211--222.
\doi{10.1016/j.chemolab.2017.11.017}

Croux, C. and Haesbroeck, G. (2003). Implementing the Bianco and Yohai estimator
for logistic regression. \emph{Computational Statistics & Data Analysis},
\bold{44}, 273--295.
}
\seealso{
\code{\link{print.enetLTS}}, \code{\link{coef.enetLTS}},
\code{\link{predict.enetLTS}}, \code{\link{outliers.enetLTS}},
\code{\link{assess.enetLTS}}, \code{\link{plot.enetLTS}},
\code{\link{cv.enetLTS}}, \code{\link{GaussianExample}},
\code{\link{BinomialExample}}, \code{\link{MultinomialExample}}
}
\examples{
## linear regression
set.seed(86)
n <- 100; p <- 25                              # observations and variables
beta <- c(rep(1, 6), rep(0, p - 6))            # six nonzero coefficients
x <- matrix(rnorm(n * p), n, p)
y <- drop(x \%*\% beta) + rnorm(n)
m <- 10                                        # contaminated observations
y[1:m] <- y[1:m] + 10                          # vertical outliers
x[1:m, ] <- x[1:m, ] + 10                      # bad leverage points

## logistic regression
set.seed(86)
n <- 100; p <- 25
beta <- c(rep(1, 6), rep(0, p - 6))
xb <- matrix(rnorm(n * p), n, p)
yb <- rbinom(n, 1, plogis(drop(xb \%*\% beta)))
mb <- 5
i0 <- which(yb == 0)[1:mb]                      # observations of class 0 ...
xb[i0, ] <- xb[i0, ] + 10                      # ... become bad leverage points

## multinomial regression
set.seed(86)
n <- 120; p <- 15; K <- 3                      # observations, variables, classes
Bm <- matrix(0, p, K)
Bm[1:2, 1] <- 3; Bm[3:4, 2] <- 3; Bm[5:6, 3] <- 3
xm <- matrix(rnorm(n * p), n, p)
pr <- exp(xm \%*\% Bm); pr <- pr / rowSums(pr)
ym <- apply(pr, 1, function(q) sample(K, 1, prob = q))
mm <- 10
xm[1:mm, 1:2] <- xm[1:mm, 1:2] + 5               # far out in the predictors of class 1 ...
ym[1:mm] <- 3                                   # ... but labeled as class 3

\donttest{
fit  <- enetLTS(x, y)
fitb <- enetLTS(xb, yb, family = "binomial")
fitm <- enetLTS(xm, ym, family = "multinomial")
fit
which(fit$wt == 0)
fitb
outliers.enetLTS(fitb)
fitm
fitm$screened
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{robust}
\keyword{regression}
\keyword{classif}
