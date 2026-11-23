# enetLTS 1.2.0

Changes since version 1.1.0.

## New functions and data
* `outliers.enetLTS()` lists the observations flagged by the raw or the reweighted
  fit with their loss, observed and predicted response and a type. For the
  classification families the type separates leverage points from possible
  mislabeled observations.
* `assess.enetLTS()`, `confusion.enetLTS()` and `roc.enetLTS()` assess a fit on
  test data. `assess.enetLTS()` returns the measures in the layout of
  `glmnet::assess.glmnet`; the deviance is twice the mean negative
  log-likelihood.
* Example data sets `GaussianExample`, `BinomialExample` and `MultinomialExample`:
  high-dimensional (p > n) data with known outliers and clean test data. The
  first two follow the simulation design of Kurnaz, Hoffmann and Filzmoser (2018).

## New arguments
* `crossfit` (gaussian and multinomial): the observations that are kept are chosen
  by cross-fitted residuals, so that an observation is never judged by a model that
  has seen it, and `alpha` and `lambda` are chosen again on the reweighted data.
  The default is `NULL`, which means `TRUE` for `family = "multinomial"` with
  `scal = TRUE` and `FALSE` otherwise. The multinomial fit therefore differs from
  version 1.1.0 and takes about 2-3 times as long; `crossfit = FALSE` switches the
  cross-fitted trimming off. For the multinomial family the reweighting step flags
  an observation when its cross-fitted probability of the observed class is below a
  cut-off calibrated to flag a share `2 * del` of the observations if the model is
  correct. The earlier rule, a median/MAD bound on the cross-fitted loss within each
  class, is available with `options(enetLTS.wrule = "mad")`.
* `screen` (multinomial, default `TRUE`): observations that are outlying in the
  predictors are found beforehand, without the class labels, and left out of the
  fit. They get weight zero and are listed in `$screened`.
* `hsize = "optimal"` with `hgrid` (experimental): the fraction of retained
  observations is chosen from `hgrid` by cross-validated trimmed loss; the result
  contains `hsize.search`.
* New component `alpha.raw` (gaussian and multinomial): the mixing parameter of the
  raw fit. With `crossfit = TRUE`, `alpha` is the value chosen again for the
  reweighted fit.
* `plotResid.enetLTS()` and `plotDiagnostic.enetLTS()` have a new argument `zoom`:
  when the flagged observations are far away, a second panel shows the regular ones.

## Changes
* `coef()` returns named coefficients (`(Intercept)`, then the column names of
  `xx` or `x1`, `x2`, ...); for the multinomial family a matrix with one column per
  class. With `zeros = FALSE` it returns the nonzero coefficients instead of their
  indices, which are available from `nonzeroCoef.enetLTS()`.
* `print()` shows only the nonzero coefficients, the selected tuning parameters and
  the number of flagged observations.
* `$objective` is the penalized criterion of the raw fit on the final subset, on
  the scale on which the C-steps fit and penalize, per observation.
* Plots: new style with a color-blind safe palette; the class is shown by color
  and shape, and flagged observations are filled. `plotCoef.enetLTS()` draws the
  slopes as needles. The cut-off lines of the gaussian residual plot are at
  +/- `qnorm(1 - del)`, and the residuals are standardized as in the reweighting
  step, so the points outside the lines are the flagged observations. The plots
  of `vers = "reweighted"` use the weights `wt`.
  `plotResid.enetLTS()` and `plotDiagnostic.enetLTS()` return the plots invisibly,
  and `plot.enetLTS()` does not change `par(ask)` or wait for input.
* The default of `crit.plot` is `FALSE`; the cross-validation heat map is built
  without the `reshape` package. The packages `reshape`, `grid` and `cvTools` are
  no longer imported.
* `predict()` accepts a data frame for `newX`. `cv.enetLTS()` accepts class labels
  for the multinomial family. `enetLTS()` stops with an error if `xx` or `yy`
  contain missing values.
* The vignette was rewritten, and all help pages were revised. The examples of the
  help pages use small simulated data; the example data sets are used in the vignette.
* Multinomial fits no longer print the `glmnet` warning that the smallest lambda
  values of a path did not converge.

## Bug fixes
* `weight.binomial()` no longer returns `NA` weights for observations with a very
  large linear predictor.
* `residuals(..., vers = "both")` for the binomial family and
  `fitted(..., vers = "both", type = "class")` for the multinomial family return
  the correct values.
* `mae` and `raw.mae` of the binomial and multinomial families are computed from the
  fitted probabilities.
* Class labels of the multinomial family that are not numbers (for example `"a"`,
  `"b"`) are coded correctly.
* A column without spread (MAD equal to zero, for example a rarely occurring dummy)
  no longer stops the fit.
* The random subsamples of the screening step and the final choice of `alpha` with
  `crossfit = TRUE` respect `seed` and the range of the user-supplied `alphas`.
* `plotDiagnostic.enetLTS()` no longer sets `options(warn = -1)`.
