\name{print.enetLTS}
\alias{print.enetLTS}
%-------------------------------------------------
\title{
print from the \code{"enetLTS"} object
}
%-------------------------------------------------
\description{
Print a summary of the \code{enetLTS} object.
}
%-------------------------------------------------
\usage{
\method{print}{enetLTS}(x,vers=c("reweighted","raw"),...)
}
%------------------------------------
\arguments{
\item{x}{fitted \code{enetLTS} object}

\item{vers}{a character string specifying for which fit to make
  predictions.  Possible values are \code{"reweighted"} (the default) for
  predicting values from the reweighted fit, \code{"raw"} for predicting
  values from the raw fit.}


\item{\dots}{additional arguments from the \code{enetLTS} object if needed.}
}

%-------------------------------------------------
\details{
  The call that produced the \code{enetLTS} object is printed, followed by
  the coefficients, the number of nonzero coefficients and penalty parameters.
}
%-------------------------------------------------
\value{
  The produced object, the coefficients, the number of nonzero coefficients and penalty parameters are returned.
}
%-------------------------------------------------
\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fatmasevinckurnaz@gmail.com>;<fskurnaz@yildiz.edu.tr>}
%-------------------------------------------------
\seealso{
\code{enetLTS},
\code{predict.enetLTS},
\code{coef.enetLTS}
}
%-------------------------------------------------
\examples{
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
yout <- c(x \%*\% beta + sigma * eout)        # response
xout <- x; xout[1:m,] <- xout[1:m,] + 10      # bad leverage points

\donttest{
fit1 <- enetLTS(xout,yout)
print(fit1)
print(fit1,vers="raw")
}

## for binomial

eps <-0.05                                     # \%10 contamination to only class 0
m <- ceiling(eps*n)
y <- sample(0:1,n,replace=TRUE)
xout <- x
xout[y==0,][1:m,] <- xout[1:m,] + 10;          # class 0
yout <- y                                      # wrong classification for vertical outliers


\donttest{
fit2 <- enetLTS(xout,yout,family="binomial")
print(fit2)
print(fit2,vers="raw")
}

## for multinomial

n <- 120; p <- 15
NC <- 3
X <- matrix(rnorm(n * p), n, p)
betas <- matrix(1:NC, ncol=NC, nrow=p, byrow=TRUE)
betas[(p-5):p,]=0; betas <- rbind(rep(0,NC),betas)
lv <- cbind(1,X) \%*\% betas
probs <- exp(lv)/apply(exp(lv),1,sum)
y <- apply(probs,1,function(prob){sample(1:NC, 1, TRUE, prob)})
xout <- X
eps <-0.05                          # \%10 contamination to only class 0
m <- ceiling(eps*n)
xout[1:m,] <- xout[1:m,] + 10       # bad leverage points
yout <- y

\donttest{
fit3    <- enetLTS(xout,yout,family="multinomial")
print(fit3)
print(fit3,vers="raw")
}
}

%-------------------------------------------------
\keyword{models}
\keyword{regression}



