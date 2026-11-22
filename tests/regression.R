## Regression tests: edge cases that used to fail (constant columns, overflow, class coding, ...).
suppressMessages(library(enetLTS))
set.seed(1)
n <- 90; y3 <- rep(1:3, each = 30)
x <- matrix(rnorm(n * 6), n)
x[, 1] <- x[, 1] + c(2, 0, -2)[y3]; x[, 2] <- x[, 2] + c(0, 2, -2)[y3]
args <- list(alphas = c(0.5, 1), lambdas = c(0.2, 0.1), nsamp = c(20, 3), seed = 1, crit.plot = FALSE)
fit <- function(...) suppressWarnings(do.call(enetLTS, c(list(...), args)))

## a column without spread (MAD = 0) must not break the fit
xd <- cbind(x, dummy = rbinom(n, 1, 0.1))
stopifnot(mad(xd[, 7]) == 0)
stopifnot(is.finite(sum(fit(xd, factor(y3), family = "multinomial")$coefficients)))
stopifnot(is.finite(sum(fit(xd, drop(x[, 1:2] %*% c(1, -1)) + rnorm(n), family = "gaussian")$coefficients)))

## deviances without overflow
eta <- rbind(c(800, 0, 0), c(0, 900, 0)); ind <- rbind(c(0, 0, 1), c(1, 0, 0))
stopifnot(all.equal(enetLTS:::multinomDeviance(eta, ind), c(800, 900)))
stopifnot(all(is.finite(enetLTS:::multinomProbs(eta))))

## a response given as a two-column indicator matrix
stopifnot(!inherits(try(fit(x, cbind(y3 == 1, y3 != 1) * 1, family = "multinomial"), silent = TRUE), "try-error"))

## class labels instead of positions
lab <- factor(c("a", "b", "c")[y3])
f <- fit(x, lab, family = "multinomial")
stopifnot(all(predict(f, newX = x[1:3, ], type = "class")$class %in% c("a", "b", "c")))
stopifnot(all(fitted(f, type = "class")$fitted.values %in% c("a", "b", "c")))

## the multinomial result does not depend on the number of cores
stopifnot(isTRUE(all.equal(f$coefficients, fit(x, lab, family = "multinomial", ncores = 2)$coefficients)))

## cross-validation loss: finite even when probabilities underflow
set.seed(3); xx <- matrix(rnorm(60 * 100), 60); xx[, 1] <- xx[, 1] + rep(c(3, 0, -3), each = 20)
Y <- diag(3)[rep(1:3, each = 20), ]
g <- glmnet::glmnet(xx, Y, family = "multinomial", alpha = 0, lambda = 1e-8, intercept = FALSE, standardize = FALSE)
stopifnot(all(is.finite(enetLTS:::cvSpec("multinomial")$loss(g, xx[1:5, ] * 400, Y[1:5, ]))))

## binomial classes: fitted and predicted classes agree with the sign of the linear predictor
yb <- as.numeric(y3 == 1); xb <- x; xb[, 1] <- xb[, 1] + 2 * yb
fb <- fit(xb, yb, family = "binomial", type.response = "class")
link <- drop(cbind(1, xb) %*% fb$coefficients)
stopifnot(all(drop(fitted(fb, type = "class")$fitted.values) == as.numeric(link > 0)))
stopifnot(all(as.numeric(predict(fb, newX = xb, type = "class")$class) == as.numeric(link > 0)))
cat("ok\n")


## hsize = "optimal"
local({
  set.seed(1); x <- matrix(rnorm(120*8), 120); y <- factor(sample(1:3, 120, TRUE)); x[, 1] <- x[, 1] + 2 * as.numeric(y); x[, 2] <- x[, 2] - 2 * as.numeric(y)
  f <- enetLTS(x, y, family="multinomial", hsize="optimal", hgrid=c(0.75,0.9),
               lambdas=c(0.1,0.05), alphas=1, nsamp=c(40,4), crit.plot=FALSE, ncores=1, seed=1)
  stopifnot(is.data.frame(f$hsize.search), sum(f$hsize.search$selected) == 1)
  r <- try(enetLTS(x, y, family="multinomial", hsize="foo", crit.plot=FALSE), silent=TRUE)
  stopifnot(inherits(r, "try-error"))
})

## class-free leverage screening (multinomial)
local({
  set.seed(3); n <- 150; x <- matrix(rnorm(n*10), n); y <- factor(rep(1:3, each = 50)); x[, 5] <- x[, 5] + 3 * (as.numeric(y) - 2); x[, 6] <- x[, 6] + c(2, -2, 0)[as.numeric(y)]
  x[1:30, 1:4] <- x[1:30, 1:4] + 10                       # 60% of class 1 (20% of all) are leverage points
  fl <- enetLTS:::leverageScreen(x)
  stopifnot(all(fl[1:30]), sum(fl[-(1:30)]) <= 5)
  f <- enetLTS(x, y, family="multinomial", lambdas=c(0.2,0.1), alphas=1, nsamp=c(40,4),
               crit.plot=FALSE, ncores=1, seed=1)
  stopifnot(length(f$wt) == n, length(f$raw.wt) == n, all(f$wt[1:30] == 0), all(f$screened %in% which(fl)))
  stopifnot(nrow(f$fitted.values) == n || length(f$fitted.values) >= n)
  f0 <- enetLTS(x, y, family="multinomial", lambdas=c(0.2,0.1), alphas=1, nsamp=c(40,4),
                crit.plot=FALSE, ncores=1, seed=1, screen=FALSE)
  stopifnot(length(f0$screened) == 0)
})

## cross-fitted trimming (gaussian and multinomial); binomial ignores it with a warning
local({
  set.seed(5); n <- 90; p <- 12; x <- matrix(rnorm(n*p), n); b <- c(2, -2, rep(0, p-2))
  y <- as.numeric(x %*% b + rnorm(n)); y[1:6] <- y[1:6] + 15
  f <- enetLTS(x, y, family="gaussian", crossfit=TRUE, lambdas=seq(0.5,0.05,length.out=6), alphas=c(0.5,1),
               nsamp=c(30,3), crit.plot=FALSE, ncores=1, seed=1)
  stopifnot(length(f$raw.wt) == n, length(f$wt) == n, all(is.finite(f$coefficients)), isTRUE(f$inputs$crossfit))
  set.seed(6); n <- 120; xm <- matrix(rnorm(n*8), n); ym <- factor(rep(1:3, each = 40))
  xm[, 1] <- xm[, 1] + 2*(as.numeric(ym) - 2); xm[, 2] <- xm[, 2] + c(2, -2, 0)[as.numeric(ym)]
  g <- enetLTS(xm, ym, family="multinomial", crossfit=TRUE, lambdas=c(0.2,0.1), alphas=c(0.5,1),
               nsamp=c(30,3), crit.plot=FALSE, ncores=1, seed=1)
  stopifnot(length(g$wt) == n, length(g$raw.wt) == n, all(is.finite(as.matrix(g$coefficients))), isTRUE(g$inputs$crossfit))
  stopifnot(mean(max.col(cbind(1, xm) %*% g$coefficients) == as.numeric(ym)) > 0.6)
  ## default: crossfit = NULL -> TRUE for multinomial (scal=TRUE), FALSE for gaussian and scal=FALSE
  g0 <- enetLTS(xm, ym, family="multinomial", lambdas=c(0.2,0.1), alphas=c(0.5,1),
                nsamp=c(30,3), crit.plot=FALSE, ncores=1, seed=1)
  stopifnot(isTRUE(g0$inputs$crossfit), isTRUE(all.equal(g0$coefficients, g$coefficients)))
  g1 <- enetLTS(xm, ym, family="multinomial", crossfit=FALSE, lambdas=c(0.2,0.1), alphas=c(0.5,1),
                nsamp=c(30,3), crit.plot=FALSE, ncores=1, seed=1)
  stopifnot(!isTRUE(g1$inputs$crossfit))
  g2 <- enetLTS(xm, ym, family="multinomial", scal=FALSE, lambdas=c(0.2,0.1), alphas=c(0.5,1),
                nsamp=c(30,3), crit.plot=FALSE, ncores=1, seed=1)
  stopifnot(!isTRUE(g2$inputs$crossfit))
  yb <- as.numeric(y3 == 1); xb <- x; xb[, 1] <- xb[, 1] + 2 * yb
  w <- NULL; h <- withCallingHandlers(do.call(enetLTS, c(list(xb, yb, family="binomial", crossfit=TRUE), args)),
         warning=function(e){ w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  stopifnot(any(grepl("gaussian and multinomial", w)))
})

## $objective is the penalised criterion of the raw fit on the final subset, on the standardised scale
local({
  chk <- function(f, family) {
    H <- f$best; a <- if (is.null(f$alpha.raw)) f$alpha else f$alpha.raw; l <- f$lambda
    xs <- enetLTS:::standardize.x(f$inputs$xx, index = H, centerFun = mean, scaleFun = sd)
    ys <- f$inputs$yy
    if (family == "gaussian") ys <- enetLTS:::center.y(ys, index = H, centerFun = mean)
    if (family == "multinomial") {
      g <- glmnet::glmnet(xs[H, ], ys[H, ], family = "multinomial", alpha = a, lambda = l, standardize = FALSE, intercept = FALSE)
      b <- do.call(cbind, lapply(g$beta, as.matrix)); eta <- xs[H, ] %*% b
      eta <- eta - apply(eta, 1, max); lp <- eta - log(rowSums(exp(eta)))
      loss <- -rowSums(ys[H, ] * lp)
    } else {
      g <- glmnet::glmnet(xs[H, ], ys[H], family = family, alpha = a, lambda = l, standardize = FALSE, intercept = FALSE)
      b <- as.matrix(g$beta); eta <- xs[H, ] %*% b
      loss <- if (family == "gaussian") 0.5 * (ys[H] - eta)^2 else -(ys[H] * eta) + log(1 + exp(eta))
    }
    ref <- mean(loss) + l * sum(0.5 * (1 - a) * b^2 + a * abs(b))
    stopifnot(isTRUE(all.equal(unname(f$objective), ref, tolerance = 1e-6)))
  }
  a2 <- list(alphas = c(0.5, 1), lambdas = c(0.2, 0.1), nsamp = c(20, 3), seed = 1, crit.plot = FALSE, ncores = 1)
  set.seed(21); n <- 80; xg <- matrix(rnorm(n*6), n) %*% diag(c(1, 5, 0.3, 2, 1, 4)); yg <- as.numeric(xg[, 1] + 0.2 * xg[, 2] + rnorm(n)); yg[1:4] <- yg[1:4] + 15
  chk(suppressWarnings(do.call(enetLTS, c(list(xg, yg, family = "gaussian"), a2))), "gaussian")
  chk(suppressWarnings(do.call(enetLTS, c(list(xg, yg, family = "gaussian", crossfit = TRUE), a2))), "gaussian")
  xb <- matrix(rnorm(n*5), n) %*% diag(c(1, 3, 0.5, 2, 1)); yb <- as.numeric(rbinom(n, 1, plogis(1.2 * xb[, 1] + 0.3 * xb[, 2])))
  chk(suppressWarnings(do.call(enetLTS, c(list(xb, yb, family = "binomial"), a2))), "binomial")
  xm <- matrix(rnorm(n*5), n) %*% diag(c(1, 3, 0.5, 2, 1)); ym <- factor(rep(1:3, length.out = n)); xm[, 1] <- xm[, 1] + 2 * (as.numeric(ym) - 2)
  chk(suppressWarnings(do.call(enetLTS, c(list(xm, ym, family = "multinomial"), a2))), "multinomial")
  chk(suppressWarnings(do.call(enetLTS, c(list(xm, ym, family = "multinomial", crossfit = FALSE), a2))), "multinomial")
})
