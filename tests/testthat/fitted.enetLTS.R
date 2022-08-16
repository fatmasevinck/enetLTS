test_that("fitted functions works", {
  ## for gaussian
  set.seed(86)
  n <- 100; p <- 25                             # number of observations and variables
  beta <- rep(0,p); beta[1:6] <- 1              # 10\% nonzero coefficients
  sigma <- 0.5                                  # controls signal-to-noise ratio
  x <- matrix(rnorm(n*p, sigma),nrow=n)
  e <- rnorm(n,0,1)                             # error terms
  eps <- 0.1                                    # contamination level
  m <- ceiling(eps*n)                           # observations to be contaminated
  eout <- e; eout[1:m] <- eout[1:m] + 10        # vertical outliers
  yout <- c(x %*% beta + sigma * eout)          # response
  xout <- x; xout[1:m,] <- xout[1:m,] + 10      # bad leverage points
  set.seed(86)
  fit1 <- enetLTS(xout,yout,crit.plot=FALSE)
  fitted.fit1 <- fitted(fit1)
  expect_equal(fitted.fit1,fit1$fitted.values)
  ## for binomial
  eps <-0.05                                     # \%10 contamination to only class 0
  m <- ceiling(eps*n)
  y <- sample(0:1,n,replace=TRUE)
  xout <- x
  xout[y==0,][1:m,] <- xout[1:m,] + 10;          # class 0
  yout <- y                                      # wrong classification for vertical outliers
  set.seed(86)
  fit2 <- enetLTS(xout,yout,family="binomial",crit.plot=FALSE,type.response="response")
  fitted.fit2 <- fitted(fit2)
  expect_equal(fitted.fit2,fit2$fitted.values)
})


