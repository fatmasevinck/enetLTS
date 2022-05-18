\name{plotDiagnostic.enetLTS}
\alias{plotDiagnostic.enetLTS}
%-------------------------------------------------
\title{
diagnostics plots from the \code{"enetLTS"} object
}
%-------------------------------------------------
\description{
  Produce plots for the diagnostics of the current model.
}
%-------------------------------------------------
\usage{
plotDiagnostic.enetLTS(object,vers=c("reweighted","raw"),...)
}
%------------------------------------
\arguments{
 \item{object}{the model fit to be plotted.}

\item{vers}{a character string denoting which model to use for the plots.
  Possible values are \code{"reweighted"} (the default) for
  plots from the reweighted fit, and \code{"raw"} for
  plots from the raw fit.}

\item{\dots}{additional arguments from the \code{enetLTS} object if needed.}
}
%-------------------------------------------------
\value{
  An object of class \code{"ggplot"} (see \code{\link[ggplot2]{ggplot}}).
}
\note{
gives the plot of

        - First two components of estimated scores for multinomial logistic regression (for \code{family="multinomial"})

        - y vs fitted values/link function. (for for both \code{family="binomial"} and \code{family="gaussian"}).

}
%-------------------------------------------------
\seealso{
\code{\link[ggplot2]{ggplot}},
\code{\link{enetLTS}},
\code{\link{coef.enetLTS}},
\code{\link{predict.enetLTS}},
\code{\link{residuals.enetLTS}},
\code{\link{fitted.enetLTS}}
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
plotDiagnostic.enetLTS(fit1)
plotDiagnostic.enetLTS(fit1,vers="raw")}

## for binomial

eps <-0.05                                     # \%10 contamination to only class 0
m <- ceiling(eps*n)
y <- sample(0:1,n,replace=TRUE)
xout <- x
xout[y==0,][1:m,] <- xout[1:m,] + 10;          # class 0
yout <- y                                      # wrong classification for vertical outliers

\donttest{
fit2 <- enetLTS(xout,yout,family="binomial",crit.plot=FALSE)
plotDiagnostic.enetLTS(fit2)
plotDiagnostic.enetLTS(fit2,vers="raw")
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
fit3 <- enetLTS(xout,yout,family="multinomial",crit.plot=FALSE)
plotDiagnostic.enetLTS(fit3)
plotDiagnostic.enetLTS(fit3,vers="raw")
}
}
%-------------------------------------------------
\references{
Kurnaz, F.S., Hoffmann, I. and Filzmoser, P. (2017) Robust and sparse estimation methods
for high dimensional linear and logistic regression. \emph{Chemometrics and Intelligent Laboratory Systems}.
}
%-------------------------------------------------
\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fatmasevinckurnaz@gmail.com>; <fskurnaz@yildiz.edu.tr>}
%-------------------------------------------------
\keyword{regression}
\keyword{classification}


