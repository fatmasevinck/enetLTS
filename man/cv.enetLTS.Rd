
\name{cv.enetLTS}

\alias{cv.enetLTS}

\title{Cross-validation for the \code{enetLTS} object}

\description{Does k-fold cross-validation for enetLTS, produces a plot,
  and returns optimal values for \code{alpha} and \code{lambda}. Combine the cross-validation
  functions internally used in the algorithm \code{enetLTS}.}

\usage{
cv.enetLTS(index=NULL,family,xx,yy,alphas,lambdas,nfold,repl,ncores,plot=TRUE)
}
\arguments{
  \item{index}{A user supplied index. The default is \code{NULL} in the algorithm \code{enetLTS.}}

  \item{family}{a description of the error distribution and link function to be used
  in the model. \code{"gaussian"} and \code{"binomial"} and \code{"multinomial"} options are available.}

  \item{xx}{matrix \code{xx} as in \code{enetLTS}.}

  \item{yy}{response \code{yy} as in \code{enetLTS}.}

  \item{alphas}{a user supplied alpha sequence for the elastic net penalty, which is
  the mixing proportion of the ridge and lasso penalties and takes value in [0,1]. Here
  \eqn{\alpha=1} is the lasso penalty, and \eqn{\alpha=0} the
  ridge penalty.}

  \item{lambdas}{a user supplied lambda sequence for the strength of the elastic net penalty.}

  \item{nfold}{a user supplied numeric value for fold number of k-fold cross-validation which
  used in varied functions of the algorithm. The default is 5-fold cross-validation.}

  \item{repl}{a user supplied posiitive number for more stable results, repeat the k-fold CV
  \code{repl} times and take the average of the corresponding evaluation measure. The default is 5.}

  \item{ncores}{a positive integer giving the number of processor cores to be used for parallel
  computing. The default is 4.}

  \item{plot}{a logical indicating if produces a plot for k-fold cross-validation
  based on alpha and lambda combinations. The default is TRUE.}
}

\note{This is an internal function. But, it is also available for direct usage to
obtain optimal values of alpha and lambda for user supplied index set.
}

\value{produces a plot,
  and returns optimal values for \code{alpha} and \code{lambda}
}

\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fskurnaz@gmail.com>;<fskurnaz@yildiz.edu.tr>
}

\keyword{models}
\keyword{regression}

