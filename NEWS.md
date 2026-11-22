# enetLTS 1.2.0

Changes since version 1.1.0 (CRAN).

## New functions and fixes
* New `outliers.enetLTS()` lists the observations flagged by the raw or reweighted
  fit (loss, observed and predicted response, type; classification families mark
  possible mislabeled observations and screened leverage points).
* New `assess.enetLTS()`, `confusion.enetLTS()` and `roc.enetLTS()` assess a fit on
  test data without relying on glmnet functions.
* `plotResid.enetLTS()` and `plotDiagnostic.enetLTS()` now use the reweighted
  weights `wt` for `vers = "reweighted"` (previously `raw.wt`), so the flagged
  points agree with `weights(fit)$wt` and `outliers.enetLTS(fit)`.
* Multinomial family, `crossfit = TRUE` (default): the reweighting step now uses
  cross-fitted class probabilities from a model fitted on all screened observations,
  with a cut-off calibrated to flag a share `2 * del` of the observations if the model
  is correct. The previous rule (median/MAD bound on the loss of the trimmed fit) flagged
  on average 10-12% of clean observations, because the fit on the trimmed subset gives
  over-confident probabilities. In simulations (n = 120, p = 15, 3 classes) the new rule
  flags 2-4% of clean observations, detects the same bad leverage points (10% and 20%
  contamination) and lowers the test error (clean data 0.205 -> 0.174, 20% contamination
  0.226 -> 0.188). With p > n (n = 100, p = 200) both rules behave similarly. The old rule
  is available with `options(enetLTS.wrule = "mad")`.
* Multinomial family: the leverage screening (`screen = TRUE`) uses random subsamples
  and was not controlled by `seed`, so the same data and `seed` could give different
  screened observations (and results) depending on the RNG state before the call.
  The RNG is now set from `seed` before the screening.
* With `crossfit = TRUE` (default for the multinomial family) the final, reweighted
  fit chose `alpha` from a fixed internal grid and ignored the user-supplied `alphas`
  (e.g. `alphas = 1` could give a final ridge-like fit). The internal grid is now
  restricted to the range of `alphas`.
* `print()` now shows only the nonzero coefficients, with variable names, together
  with the selected tuning parameters and the number of flagged observations.
* New example data sets `GaussianExample`, `BinomialExample` and `MultinomialExample`
  (high-dimensional, with known outliers and clean test data); the first two follow
  the simulation design of Kurnaz, Hoffmann and Filzmoser (2018).
* `assess.enetLTS()` now returns the measures in the layout of `glmnet::assess.glmnet`
  (one list element per measure, reweighted and raw fit side by side), and the
  deviance is twice the mean negative log-likelihood, as in glmnet.
* `outliers.enetLTS()`: the loss of the classification families is computed from the
  linear predictor (it was capped at -log(1e-15) = 34.5); screened observations are
  labeled `bad leverage point (screened)` when they are misclassified.
* Plots: new style with a color-blind safe palette. `plotCoef.enetLTS()` draws the
  slopes as needles (one panel per class for the multinomial family) with the
  intercept in the subtitle. In the plots of the classification families the class
  is shown by color and shape and flagged observations are filled. The cut-off
  lines of the gaussian residual plot are at +/- `qnorm(1 - del)` (they were fixed
  at +/- 2.5) and the residuals are standardized as in the reweighting step, so
  the points outside the lines are exactly the flagged observations. The residual
  plots of the classification families are labeled "negative log-likelihood".
  `roc.enetLTS(..., plot = TRUE)` draws a ggplot ROC curve with the AUC.
  New argument `zoom` of `plotResid.enetLTS()` and `plotDiagnostic.enetLTS()`: when
  the flagged observations are far away, a second panel zooms in on the regular ones.
* `plotDiagnostic.enetLTS()` no longer sets `options(warn = -1)` for the session.
* The vignette was rewritten: the same structure for every family, the example data
  sets, interpretation of the output, tables of the main arguments.
* `weight.binomial()` no longer returns `NA` weights for observations with a very
  large linear predictor (`Inf/Inf`); such observations now get a proper 0/1 weight.

* `coef()` returns the coefficients with names (`(Intercept)`, then the column names
  of `xx`, or `x1`, `x2`, ...); for the multinomial family a matrix with one column
  per class. `zeros = FALSE` now returns the nonzero coefficients (previously their
  indices, which are available with `nonzeroCoef.enetLTS()`).
* `predict()` accepts a data frame for `newX`, and for the multinomial family a
  single new observation (`type = "link"` failed before).
* `fitted(..., vers = "both", type = "class")` returned class numbers instead of
  class labels for the reweighted fit of the multinomial family; fixed.
* `residuals(..., vers = "both")` failed for the binomial family (it referred to a
  component that does not exist); fixed.
* `cv.enetLTS()` accepts class labels for the multinomial family (an indicator matrix
  was required), uses the same default `lambdas` as `enetLTS()` for that family, and
  returns its result visibly.
* `plot.enetLTS()` no longer waits for <Enter> before the first plot and no longer
  changes `par(ask)` for the session; in an interactive session R asks before each
  new page when all three plots are drawn.
* `plotResid.enetLTS()` returns the drawn plots (a list), `plotDiagnostic.enetLTS()`
  the drawn plot, invisibly.
* `mae` and `raw.mae` of the binomial and multinomial families are now computed from
  the fitted probabilities (they used the linear predictor).
* Class labels of the multinomial family that are not numbers (e.g. `"a"`, `"b"`)
  are coded correctly.
* The cross-validation heat map (`crit.plot = TRUE`, `cv.enetLTS(plot = TRUE)`) is
  built without `reshape`; the packages `reshape`, `grid` and `cvTools` are no longer
  imported.
* All help pages were revised: every argument, default and return value was checked
  against the code, and the examples use the example data sets.

## Multinomial: cross-fitted trimming is now the default (`crossfit = NULL`)
* `crossfit` now defaults to `NULL`: `TRUE` for `family = "multinomial"` with
  `scal = TRUE`, `FALSE` for every other case (gaussian, `scal = FALSE`;
  binomial ignores it). `crossfit = FALSE` switches the cross-fitted trimming off.
  Code that calls `enetLTS(..., family = "multinomial")` therefore returns
  different estimates than version 1.1.0 and takes about 2-3 times as long.
* Evidence (paper design with 10 variants: signal strength, class imbalance,
  4 classes, correlation, n = 90/180/450, dense and high-dimensional; contamination
  rows placed at random; 20 paired replications per cell; test MCR of the reweighted
  fit; difference to classical elastic net on the same data):
  - systematic label errors (30% of class 1 relabeled, plus leverage shift):
    average -0.017 (95% CI -0.022, -0.013); 6 of 10 variants significantly better,
    none worse;
  - 20% leverage contamination: -0.007 (-0.011, -0.003); dense model -0.076;
  - clean data: +0.011 (0.008, 0.013); 5 of 10 variants significantly worse
    (correlated +0.038, large n +0.023, strong signal +0.016, base +0.010, weak +0.007);
  - `crossfit = TRUE` against `crossfit = FALSE`: lower error in all four
    scenarios (-0.009 to -0.022), and far fewer false positives under systematic
    label errors (the default `FALSE` selects many noise variables there).
* New output `alpha.raw` (gaussian, multinomial): the mixing parameter of the raw fit. With
  `crossfit = TRUE`, `alpha` is the value re-selected for the reweighted fit while
  `lambda` and `objective` belong to the raw fit.
* Known limits: on clean data the trimming costs accuracy relative to classical
  elastic net, most with highly correlated predictors and large n; with class
  imbalance plus leverage contamination (`lev3`) the classical fit was better by 0.035;
  the fit is sparser than classical (false negative rate about 0.6).

## `$objective` (all three families)
* The reported `objective` was neither the criterion of the C-steps nor the
  criterion of any single fit: it combined the coefficients of the *reweighted*
  fit with the subset and lambda of the *raw* fit, penalized the coefficients on
  the *original* scale (the C-steps penalize the standardized ones, i.e.
  `s_H * beta`), and was multiplied by `h`. In a gaussian example (columns with
  sd 1, 5, 0.3) the penalty was 0.265 instead of 0.168. It is an output only,
  nothing in the algorithm uses it, so estimates and predictions are unchanged.
* `objective` is now `mean(loss over H) + lambda * P_alpha(s_H * beta)` for the
  raw fit on the final subset `H`, with the (alpha, lambda) of that fit (also
  with `crossfit = TRUE`, where alpha is re-selected later). Loss is half the
  squared residual (gaussian) or the deviance (binomial, multinomial). Values
  are therefore per observation (previously times `h`) and not comparable with
  earlier versions.
* For the binomial family the C-step ranks subsets by a robust criterion
  without penalty (`objBinom`); `objective` is the penalized deviance, so it
  is not the number the search optimises. This is documented in `?enetLTS`.
* Not changed, noted: for `family = "multinomial"` and `scal = FALSE` the C-steps
  run on the robustly standardized `x` while the final raw fit uses the
  unscaled `xx`.

## Gaussian: optional cross-fitted trimming (`crossfit = TRUE`, default `FALSE`)
* On the linear design of the paper (n = 150, p = 60) the reweighted subset of
  the ordinary C-steps is selected partly *because* noise variables explain
  it: classical elastic net on the 116 retained observations had false
  positive rate 0.86 and RMSPE 1.46, against 0.14 and 1.09 on 116 randomly
  chosen observations. The choice of alpha/lambda was not the cause.
* New argument `crossfit` (gaussian and multinomial; for multinomial see below): residuals that decide which
  observations are kept are cross-fitted, the reweighting set is taken from
  them, and alpha and lambda are re-selected on the reweighted data with
  common folds and the one-standard-error rule.

## Multinomial: leverage points concentrated in one class
* The class-wise robust distances of the reweighting step (and of the C-steps)
  break down when most of one class consists of leverage points: the MCD of that
  class then describes the contamination, the clean observations of the class
  are removed and the fit is wrong. With 20% leverage points in the first
  rows, i.e. mostly in class 1, the misclassification rate of the paper design
  (n = 180, p = 55) was 0.237 (both blocks) and 0.166 (informative block).
* New argument `screen` (default `TRUE`, multinomial only): leverage points are
  found beforehand without the class labels (robust standardisation,
  spatial-sign PCA, robust distances of the scores; at most 25% of the
  observations) and left out of the fit. They get weight zero in `raw.wt` and
  `wt` and are listed in `$screened`. The misclassification rate of the two
  cases above drops to 0.117, the rate of the clean data. `screen = FALSE`
  gives the previous behavior.

