\name{nonzeroCoef.enetLTS}
\alias{nonzeroCoef.enetLTS}

\title{
nonzero coefficients indices from the \code{"enetLTS"} object
}
\description{
A numeric vector which gives the indices of nonzero coefficients from the current model.
}
\usage{
nonzeroCoef.enetLTS(object,vers=c("reweighted","raw"))
}
\arguments{
\item{object}{the model fit from which to extract nonzero coefficients indices.}

\item{vers}{a character string denoting which model to use.
  Possible values are \code{"reweighted"} (the default) for
  plots from the reweighted fit, and \code{"raw"} for
  plots from the raw fit.}
}
\value{
A numeric vector (or a list object for family="multinomial") containing the request.
}
\seealso{
\code{\link{enetLTS}},
\code{\link{predict.enetLTS}},
\code{\link{coef.enetLTS}}
}

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
nonzeroCoef.enetLTS(fit1)
nonzeroCoef.enetLTS(fit1,vers="raw")
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
nonzeroCoef.enetLTS(fit2)
nonzeroCoef.enetLTS(fit2,vers="raw")
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
nonzeroCoef.enetLTS(fit3)
nonzeroCoef.enetLTS(fit3,vers="raw")

}
}


\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fatmasevinckurnaz@gmail.com>;<fskurnaz@yildiz.edu.tr>}

\keyword{regression}
\keyword{classification}



