\name{predict.enetLTS}
\alias{predict.enetLTS}
%-------------------------------------------------
\title{
make predictions from the \code{"enetLTS"} object.
}
%-------------------------------------------------
\description{
Similar to other predict methods, this function predicts fitted values, logits,
coefficients and nonzero coefficients from a fitted \code{"enetLTS"} object.
}
%-------------------------------------------------
\usage{
\method{predict}{enetLTS}(object,newX,vers=c("reweighted","raw"),
    type=c("link","response","coefficients","nonzero","class"),...)
}
%------------------------------------
\arguments{
\item{object}{the model fit from which to make predictions.}

\item{newX}{new values for the predictor matrix \code{X}.
  Must be a matrix; can be sparse as in \code{Matrix} package.
  This argument is not used for \code{type=c("coefficients","nonzero")}.}

\item{vers}{a character string denoting which fit to use for the predictions.
  Possible values are \code{"reweighted"} (the default) for
  predicting values from the reweighted fit, \code{"raw"} for predicting
  values from the raw fit.}

\item{type}{type of prediction required. \code{type="link"} gives the link function.
  \code{type="response"} gives the
  fitted probabilities for \code{"binomial"} and gives the fitted values for
  \code{"gaussian"}. \code{type="coefficients"} computes the coefficients from the
  fitted model. \code{type="nonzero"} returns a list of the indices of the nonzero
  coefficients. \code{type="class"} is available only for \code{"binomial"} model,
  and produces the class label corresponding to the maximum probability.}

\item{\dots}{additional arguments from the \code{enetLTS} object if needed.}
}
%-------------------------------------------------
\details{
The \code{newdata} argument defaults to the matrix of predictors used to fit
the model such that the fitted values are computed.

\code{coef.enetLTS(...)} is equivalent to \code{predict.enetLTS(object,newX,type="coefficients",...)}, where newX argument is the matrix as in \code{enetLTS}.
}
%-------------------------------------------------
\value{
The requested predicted values are returned.
}
%-------------------------------------------------
\seealso{
\code{\link{enetLTS}},
\code{\link{coef.enetLTS}},
\code{\link{nonzeroCoef.enetLTS}}
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
predict(fit1,newX=xout)
predict(fit1,newX=xout,type="coefficients")
predict(fit1,newX=xout,type="nonzero",vers="raw")
# provide new X matrix
newX <- matrix(rnorm(n*p, sigma),nrow=n)
predict(fit1,newX=newX,type="response")
predict(fit1,newX=newX,type="coefficients")
predict(fit1,newX=newX,type="nonzero")}

## for binomial

eps <-0.05                                     # \%10 contamination to only class 0
m <- ceiling(eps*n)
y <- sample(0:1,n,replace=TRUE)
xout <- x
xout[y==0,][1:m,] <- xout[1:m,] + 10;          # class 0
yout <- y                                      # wrong classification for vertical outliers

\donttest{
fit2 <- enetLTS(xout,yout,family="binomial")
predict(fit2,newX=xout)
predict(fit2,newX=xout,type="coefficients")
predict(fit2,newX=xout,type="nonzero",vers="raw")
predict(fit2,newX=newX,type="class")
predict(fit2,newX=newX,type="coefficients",vers="raw")
predict(fit2,newX=newX,type="nonzero")
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
predict(fit3,newX=xout)
predict(fit3,newX=xout,type="coefficients")
predict(fit3,newX=xout,type="nonzero",vers="raw")
predict(fit3,newX=newX,type="class")
predict(fit3,newX=newX,type="coefficients",vers="raw")
predict(fit3,newX=newX,type="nonzero")
}
}

%-------------------------------------------------
\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fatmasevinckurnaz@gmail.com>;<fskurnaz@yildiz.edu.tr>}
%-------------------------------------------------
\keyword{regression}
\keyword{classification}


