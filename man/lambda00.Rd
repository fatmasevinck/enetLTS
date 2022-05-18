\name{lambda00}
\alias{lambda00}
%--------------------------
\title{Upper limit of the penalty parameter for \code{family="binomial"}}
%--------------------------
\description{
Use bivariate winsorization to estimate the smallest value of the upper limit for the penalty
parameter.
}
%--------------------------
\usage{
lambda00(x,y,normalize=TRUE,intercept=TRUE,const=2,prob=0.95,
      tol=.Machine$double.eps^0.5,eps=.Machine$double.eps,...)
}
%--------------------------
\arguments{
\item{x}{a numeric matrix containing the predictor variables.}

\item{y}{a numeric vector containing the response variable.}

\item{normalize}{a logical indicating whether the winsorized predictor
variables should be normalized or not (the
default is \code{TRUE}).}

\item{intercept}{a logical indicating whether a constant term should be
included in the model (the default is \code{TRUE}).}

\item{const}{numeric; tuning constant to be used in univariate
winsorization (the default is 2).}

\item{prob}{numeric; probability for the quantile of the
\eqn{\chi^{2}}{chi-squared} distribution to be used in bivariate
winsorization (the default is 0.95).}

\item{tol}{a small positive numeric value used to determine singularity
issues in the computation of correlation estimates for bivariate
winsorization.}

\item{eps}{a small positive numeric value used to determine whether the
robust scale estimate of a variable is too small (an effective zero).}

\item{\dots}{additional arguments if needed.}
}
%--------------------------
\value{
A robust estimate of the smallest value of the penalty parameter for
enetLTS regression (for \code{family="binomial"}).
}
%--------------------------
\details{
The estimation procedure is done with similar approach as in Alfons et al. (2013).
But the Pearson correlation between y and the jth predictor variable xj on winsorized data is
replaced to a robustified point-biserial correlation for logistic regression.
}
\note{
For linear regression, we take exactly same procedure as in Alfons et al., which is based on
the Pearson correlation between y and the jth predictor variable xj on winsorized
data. See Alfons et al. (2013).
}
%--------------------------
\examples{
set.seed(86)
n <- 100; p <- 25                             # number of observations and variables
beta <- rep(0,p); beta[1:6] <- 1              # 10\% nonzero coefficients
sigma <- 0.5                                  # controls signal-to-noise ratio
x <- matrix(rnorm(n*p, sigma),nrow=n)
e <- rnorm(n,0,1)                             # error terms
eps <-0.05                                    # \%10 contamination to only class 0
m <- ceiling(eps*n)
y <- sample(0:1,n,replace=TRUE)
xout <- x
xout[y==0,][1:m,] <- xout[1:m,] + 10;         # class 0
yout <- y                                     # wrong classification for vertical outliers

# compute smallest value of the upper limit for the penalty parameter
l00 <- lambda00(xout,yout)
}
%-----------------------------------
\references{
Kurnaz, F.S., Hoffmann, I. and Filzmoser, P. (2017) Robust and sparse estimation methods
for high dimensional linear and logistic regression.
\emph{Chemometrics and Intelligent Laboratory Systems.}

Alfons, A., Croux, C. and Gelper, S. (2013) Sparse least trimmed squares regression for
analyzing high-dimensional large data sets. \emph{The Annals of Applied Statistics}, 7(1), 226--248.
}
%-------------------------------------------------
\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fatmasevinckurnaz@gmail.com>;<fskurnaz@yildiz.edu.tr>}
%-----------------------------------
\seealso{
\code{\link{enetLTS}},
\code{\link[robustHD]{sparseLTS}},
\code{\link[robustHD]{lambda0}}
}
\keyword{robust}
