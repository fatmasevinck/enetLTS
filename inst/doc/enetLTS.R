## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 8, fig.height = 3.6,
                      out.width = "100%", message = FALSE, warning = FALSE)

## ----install, eval = FALSE----------------------------------------------------
#  install.packages("enetLTS")

## ----load---------------------------------------------------------------------
library(enetLTS)

## ----qs-data------------------------------------------------------------------
data(GaussianExample)
x <- GaussianExample$x
y <- GaussianExample$y
dim(x)

## ----qs-fit-------------------------------------------------------------------
fit <- enetLTS(x, y, seed = 1)

## ----qs-print-----------------------------------------------------------------
fit

## ----qs-outliers--------------------------------------------------------------
outliers.enetLTS(fit)

## ----qs-predict---------------------------------------------------------------
predict(fit, newX = GaussianExample$xtest[1:5, ], type = "response")

## ----gaus-nonzero-------------------------------------------------------------
b  <- coef(fit)                         # intercept first
br <- coef(fit, vers = "raw")
c(reweighted = sum(b[-1] != 0), raw = sum(br[-1] != 0))

## ----gaus-ids-----------------------------------------------------------------
o   <- outliers.enetLTS(fit)
fl  <- o$index[o$flagged]
out <- GaussianExample$outliers
list(flagged = fl,
     detected = intersect(fl, out),
     missed = setdiff(out, fl),
     clean_but_flagged = setdiff(fl, out))

## ----del, eval = FALSE--------------------------------------------------------
#  fit.del <- enetLTS(x, y, del = 0.005, seed = 1)

## ----crossfit-----------------------------------------------------------------
fit.cf <- enetLTS(x, y, crossfit = TRUE, seed = 1)
ocf <- outliers.enetLTS(fit.cf)
flcf <- ocf$index[ocf$flagged]
list(detected = intersect(flcf, out), clean_but_flagged = setdiff(flcf, out))
sapply(list(default = fit, crossfit = fit.cf), function(f)
  assess.enetLTS(f, GaussianExample$xtest, GaussianExample$ytest, vers = "reweighted")$rmse)

## ----hsize, eval = FALSE------------------------------------------------------
#  fit60 <- enetLTS(x, y, hsize = 0.6, seed = 1)

## ----parallel, eval = FALSE---------------------------------------------------
#  fit.par <- enetLTS(x, y, ncores = 4, seed = 1)

## ----gaus-pred----------------------------------------------------------------
head(predict(fit, newX = GaussianExample$xtest, type = "response")$response)

## ----gaus-coef----------------------------------------------------------------
plotCoef.enetLTS(fit)

## ----gaus-resid---------------------------------------------------------------
plotResid.enetLTS(fit)

## ----gaus-diag----------------------------------------------------------------
plotDiagnostic.enetLTS(fit, zoom = FALSE)

## ----gaus-assess--------------------------------------------------------------
assess.enetLTS(fit, GaussianExample$xtest, GaussianExample$ytest)

## ----bin-fit------------------------------------------------------------------
data(BinomialExample)
xb <- BinomialExample$x
yb <- BinomialExample$y
table(yb)
fitb <- enetLTS(xb, yb, family = "binomial", seed = 1)
fitb

## ----bin-out------------------------------------------------------------------
ob <- outliers.enetLTS(fitb)
ob
outb <- BinomialExample$outliers
list(contaminated = outb,
     detected = intersect(ob$index[ob$flagged], outb),
     clean_but_flagged = setdiff(ob$index[ob$flagged], outb))

## ----bin-resid----------------------------------------------------------------
plotResid.enetLTS(fitb, zoom = FALSE)

## ----bin-diag-----------------------------------------------------------------
plotDiagnostic.enetLTS(fitb)

## ----bin-assess---------------------------------------------------------------
assess.enetLTS(fitb, BinomialExample$xtest, BinomialExample$ytest)

## ----bin-conf-----------------------------------------------------------------
confusion.enetLTS(fitb, BinomialExample$xtest, BinomialExample$ytest)

## ----bin-roc, fig.height = 4.5------------------------------------------------
roc.enetLTS(fitb, BinomialExample$xtest, BinomialExample$ytest, plot = TRUE)

## ----mult-fit-----------------------------------------------------------------
data(MultinomialExample)
xm <- MultinomialExample$x
ym <- MultinomialExample$y
table(ym)
fitm <- enetLTS(xm, ym, family = "multinomial", seed = 1)
fitm

## ----mult-out-----------------------------------------------------------------
om <- outliers.enetLTS(fitm)
om
outm <- MultinomialExample$outliers
list(detected = intersect(om$index[om$flagged], outm),
     missed = setdiff(outm, om$index[om$flagged]),
     clean_but_flagged = setdiff(om$index[om$flagged], outm))

## ----mult-coef, fig.height = 5.5----------------------------------------------
plotCoef.enetLTS(fitm)

## ----mult-resid---------------------------------------------------------------
plotResid.enetLTS(fitm, zoom = FALSE)

## ----mult-diag----------------------------------------------------------------
plotDiagnostic.enetLTS(fitm)

## ----mult-assess--------------------------------------------------------------
assess.enetLTS(fitm, MultinomialExample$xtest, MultinomialExample$ytest)
confusion.enetLTS(fitm, MultinomialExample$xtest, MultinomialExample$ytest)

## ----mult-classwise-----------------------------------------------------------
pm <- predict(fitm, newX = MultinomialExample$xtest, type = "class")$class
round(tapply(pm == MultinomialExample$ytest, MultinomialExample$ytest, mean), 3)

## ----hsize-opt, eval = FALSE--------------------------------------------------
#  fit.opt <- enetLTS(x, y, hsize = "optimal", seed = 1)

