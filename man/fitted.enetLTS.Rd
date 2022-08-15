\name{fitted.enetLTS}
\alias{fitted.enetLTS}
%-------------------------------------------------
\title{
 the fitted values from the \code{"enetLTS"} object.
}
%-------------------------------------------------
\description{
  A numeric vector which extract fitted values from the current model.
}
%-------------------------------------------------
\usage{
  \method{fitted}{enetLTS}(object,vers=c("reweighted","raw"),type=c("response","class"),...)
}
%-------------------------------------------------
\arguments{
\item{object}{the model fit from which to extract fitted values.}

\item{vers}{a character string specifying for which fit to make
  predictions. Possible values are \code{"reweighted"} (the default) for
  predicting values from the reweighted fit, \code{"raw"} for predicting
  values from the raw fit.}

\item{type}{type of prediction required. \code{type="response"} gives the
  fitted probabilities for \code{"multinomial"} and \code{"binomial"} and gives the fitted values
  for \code{"gaussian"}. \code{type="class"} is available only for \code{"multinomial"} and
  \code{"binomial"} model, and produces the class label corresponding to the maximum probability.}

\item{\dots}{additional arguments from the \code{enetLTS} object if needed.}
}
%-------------------------------------------------
\value{
A numeric vector containing the requested fitted values.
}
%-------------------------------------------------
\seealso{
\code{\link{enetLTS}},
\code{\link{predict.enetLTS}},
\code{\link{residuals.enetLTS}}
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
fit1 <- enetLTS(xout,yout,crit.plot=FALSE)
fitted(fit1)
fitted(fit1,vers="raw")
fitted(fit1,vers="reweighted",type="response")
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
fitted(fit2)
fitted(fit2,vers="raw")
fitted(fit2,vers="reweighted",type="class")
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
fit3 <- enetLTS(xout,yout,family="multinomial")
fitted(fit3)
fitted(fit3,vers="raw")
fitted(fit3,vers="reweighted",type="class")
}
}

\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fskurnaz@gmail.com>;<fskurnaz@yildiz.edu.tr>}

\keyword{regression}
\keyword{classification}




