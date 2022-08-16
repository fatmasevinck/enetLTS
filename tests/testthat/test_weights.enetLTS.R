library(testthat)  
library(enetLTS)

test_that("weights functions works", {
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
  weights.fit11 <- weights.enetLTS(fit1)
  expect_equal( weights.fit11,fit1$wt)
  weights.fit12 <- weights.enetLTS(fit1,vers="raw")
  expect_equal( weights.fit12,fit1$raw.wt)
  ## for binomial
  eps <-0.05                                     # \%10 contamination to only class 0
  m <- ceiling(eps*n)
  y <- sample(0:1,n,replace=TRUE)
  xout <- x
  xout[y==0,][1:m,] <- xout[1:m,] + 10;          # class 0
  yout <- y                                      # wrong classification for vertical outliers
  set.seed(86)
  fit2 <- enetLTS(xout,yout,family="binomial",crit.plot=FALSE)
  weights.fit21 <- weights.enetLTS(fit2)
  expect_equal( weights.fit21,fit2$wt)
  weights.fit22 <- weights.enetLTS(fit2,vers="raw")
  expect_equal(weights.fit22,fit2$raw.wt)
  ## for multinomial
  n <- 120; p <- 15 
  NC <- 3                                        # number of groups
  X <- matrix(rnorm(n * p), n, p)               
  betas <- matrix(1:NC, ncol=NC, nrow=p, byrow=TRUE)
  betas[(p-5):p,]=0; betas <- rbind(rep(0,NC),betas)
  lv <- cbind(1,X)%*%betas                      
  probs <- exp(lv)/apply(exp(lv),1,sum)         
  y <- apply(probs,1,function(prob){sample(1:NC, 1, TRUE, prob)})
  xout <- X                                     
  eps <-0.05                                    # \%10 contamination to only class 0
  m <- ceiling(eps*n)                           
  xout[1:m,] <- xout[1:m,] + 10                 # bad leverage points
  yout <- y 
  fit3 <- enetLTS(xout,yout,family="multinomial",crit.plot=FALSE)
  weights.fit31 <- weights.enetLTS(fit3)
  expect_equal( weights.fit31,fit3$wt)
  weights.fit32 <- weights.enetLTS(fit3,vers="raw")
  expect_equal(weights.fit32,fit3$raw.wt)
})
