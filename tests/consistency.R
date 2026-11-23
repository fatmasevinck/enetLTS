## Results must not depend on the number of cores or on the parallel back end.
suppressMessages(library(enetLTS))
set.seed(1)
n <- 60; p <- 8
x <- matrix(rnorm(n * p), n, p)
y <- drop(1 + x[, 1:2] %*% c(1, -1) + rnorm(n))
yb <- as.numeric(y > 1)

fit <- function(family, y, ncores, backend = NULL) {
  options(enetLTS.backend = backend)
  on.exit(options(enetLTS.backend = NULL))
  enetLTS(x, y, family = family, alphas = c(0.5, 1), lambdas = c(0.2, 0.1, 0.05),
          nsamp = c(20, 3), seed = 1, ncores = ncores, crit.plot = FALSE)
}

for (family in c("gaussian", "binomial")) {
  yy <- if (family == "gaussian") y else yb
  ref <- fit(family, yy, 1)
  stopifnot(isTRUE(all.equal(ref$coefficients, fit(family, yy, 2)$coefficients)))
  if (.Platform$OS.type != "windows")
    stopifnot(isTRUE(all.equal(ref$coefficients, fit(family, yy, 2, "psock")$coefficients)))
}
cat("ok\n")
