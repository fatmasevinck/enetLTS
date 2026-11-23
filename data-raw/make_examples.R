## Generates the example data sets GaussianExample, BinomialExample and
## MultinomialExample (run from the package root directory).
## Each data set is a list with the training data x, y, the indices of the
## contaminated observations (outliers), and clean test data xtest, ytest.

r3 <- function(x) round(x, 3)

## ---- GaussianExample: design of Kurnaz, Hoffmann and Filzmoser (2018) ----
## Two informative blocks (5% of p each, rho = 0.9) and an uninformative block
## (90%, rho = 0.2); y = 1 + x'beta + e. The first 10% of the observations get
## informative predictors from N(20, 1) and errors from N(20 * sd(y), 1)
## (bad leverage points that are also vertical outliers).
set.seed(1)
n <- 100; p <- 200; ntest <- 200
p1 <- 10; p2 <- 10; pb <- p - p1 - p2; pa <- p1 + p2
toep <- function(k, r) r^abs(outer(1:k, 1:k, "-"))
R1 <- chol(toep(p1, 0.9)); Rb <- chol(toep(pb, 0.2))
gen <- function(N) cbind(matrix(rnorm(N * p1), N, p1) %*% R1, matrix(rnorm(N * p2), N, p2) %*% R1,
                         matrix(rnorm(N * pb), N, pb) %*% Rb)
beta <- c(rep(1, pa), rep(0, pb))
x <- gen(n); y <- 1 + drop(x %*% beta) + rnorm(n); sy <- sd(y)
xtest <- gen(ntest); ytest <- 1 + drop(xtest %*% beta) + rnorm(ntest)
outliers <- 1:10
x[outliers, 1:pa] <- matrix(rnorm(length(outliers) * pa, 20, 1), length(outliers), pa)
y[outliers] <- 1 + drop(x[outliers, ] %*% beta) + rnorm(length(outliers), 20 * sy, 1)
GaussianExample <- list(x = r3(x), y = r3(y), outliers = outliers,
                        xtest = r3(xtest), ytest = r3(ytest), beta = beta)

## ---- BinomialExample: design of Kurnaz, Hoffmann and Filzmoser (2018) ----
## Correlated blocks (rho = 0.9 informative, 0.5 uninformative), latent model
## y = 1(1 + x'beta + e > 0); the first 10% of class 0 get informative
## predictors from N(20, 1) and keep label 0 (bad leverage points).
set.seed(1)
n <- 100; p <- 200; ntest <- 200; pa <- 20; pb <- p - pa
toep <- function(k, r) r^abs(outer(1:k, 1:k, "-"))
Ra <- chol(toep(pa, 0.9)); Rb <- chol(toep(pb, 0.5))
gen <- function(N) cbind(matrix(rnorm(N * pa), N, pa) %*% Ra, matrix(rnorm(N * pb), N, pb) %*% Rb)
beta <- c(rep(1, pa), rep(0, pb))
x <- gen(n); y <- as.integer(1 + x %*% beta + rnorm(n) > 0)
xtest <- gen(ntest); ytest <- as.integer(1 + xtest %*% beta + rnorm(ntest) > 0)
i0 <- which(y == 0); outliers <- i0[seq_len(floor(0.1 * length(i0)))]
x[outliers, 1:pa] <- matrix(rnorm(length(outliers) * pa, 20, 1), length(outliers), pa)
BinomialExample <- list(x = r3(x), y = y, outliers = outliers,
                        xtest = r3(xtest), ytest = ytest, beta = beta)

## ---- MultinomialExample: n = 150, p = 200, K = 3, 10% bad leverage points ----
set.seed(1)
n <- 150; p <- 200; K <- 3; ntest <- 300
B <- matrix(0, p, K); B[1:2, 1] <- 3; B[3:4, 2] <- 3; B[5:6, 3] <- 3
sim <- function(X) { lv <- X %*% B; pr <- exp(lv) / rowSums(exp(lv))
  apply(pr, 1, function(q) sample(seq_len(K), 1, prob = q)) }
x <- matrix(rnorm(n * p), n, p); y <- sim(x)
xtest <- matrix(rnorm(ntest * p), ntest, p); ytest <- sim(xtest)
outliers <- 1:15
x[outliers, 1:2] <- x[outliers, 1:2] + 5     # far out in the predictors of class 1 ...
y[outliers] <- 3                              # ... but labeled as class 3
MultinomialExample <- list(x = r3(x), y = factor(y), outliers = outliers,
                           xtest = r3(xtest), ytest = factor(ytest, levels = 1:3), beta = B)

save(GaussianExample, file = "data/GaussianExample.rda", compress = "xz")
save(BinomialExample, file = "data/BinomialExample.rda", compress = "xz")
save(MultinomialExample, file = "data/MultinomialExample.rda", compress = "xz")
