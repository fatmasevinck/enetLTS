\name{assess.enetLTS}
\alias{assess.enetLTS}
\alias{confusion.enetLTS}
\alias{roc.enetLTS}
\alias{print.confusion.enetLTS}
\title{Assess an \code{"enetLTS"} fit on test data}
\description{
Performance measures, confusion matrix and ROC curve of a fitted
\code{"enetLTS"} object on (new) data.
}
\usage{
assess.enetLTS(object, newX, newy, vers = c("both", "reweighted", "raw"))

confusion.enetLTS(object, newX, newy, vers = c("reweighted", "raw"))

roc.enetLTS(object, newX, newy, vers = c("reweighted", "raw"), plot = FALSE, ...)

\method{print}{confusion.enetLTS}(x, ...)
}
\arguments{
\item{object}{a fitted object of class \code{"enetLTS"}.}
\item{newX}{a numeric matrix with the predictors of the test data.}
\item{newy}{the response of the test data: numeric for \code{"gaussian"}, class
  labels (as in the training data) for \code{"binomial"} and \code{"multinomial"}.}
\item{vers}{which fit to assess: \code{"reweighted"} (the default of
  \code{confusion.enetLTS} and \code{roc.enetLTS}) or \code{"raw"}.
  \code{assess.enetLTS} also accepts \code{"both"} (its default), which reports
  the two fits side by side.}
\item{plot}{logical; if \code{TRUE} the ROC curve is drawn. The default is \code{FALSE}.}
\item{x}{an object of class \code{"confusion.enetLTS"}.}
\item{\dots}{not used.}
}
\details{
The test data should be free of outliers, since they stand for the data on which
the model is used.

\code{assess.enetLTS} computes the mean squared error \code{mse}, its root
\code{rmse} and the mean absolute error \code{mae} (gaussian); the deviance
\code{deviance}, the misclassification rate \code{class}, \code{auc}, \code{mse}
and \code{mae} of the predicted probabilities (binomial); and the same without
\code{auc} (multinomial; \code{mse} and \code{mae} are summed over the classes).
As in \code{glmnet::assess.glmnet}, \code{deviance} is twice the mean negative
log-likelihood. For \code{"binomial"} the second class of \code{object$classnames}
is the positive class, and an observation is classified as positive if its
probability is larger than 0.5.

\code{confusion.enetLTS} and \code{roc.enetLTS} are available for the
classification families (\code{roc.enetLTS} for \code{"binomial"} only).
}
\value{
\code{assess.enetLTS}: a list with one element per measure, in the layout of
\code{glmnet::assess.glmnet}. Every element is a named vector with one value per
fit (\code{reweighted}, \code{raw}) and has the attribute \code{"measure"} with the
name of the measure.

\code{confusion.enetLTS}: an object of class \code{"confusion.enetLTS"}, a list
with the contingency table \code{table} (predicted by true class) and the
proportion of correct classifications \code{correct}.

\code{roc.enetLTS}: a data frame with the columns \code{FPR} and \code{TPR}, with
the AUC as attribute \code{"auc"}, invisibly.
}
\references{
Hastie, T., Qian, J. and Tay, K. (2026). An Introduction to \code{glmnet}.
Vignette of the R package \pkg{glmnet}.
\url{https://glmnet.stanford.edu/articles/glmnet.html}
}
\seealso{
\code{\link{enetLTS}}, \code{\link{predict.enetLTS}}, \code{\link{outliers.enetLTS}}
}
\examples{
data("GaussianExample")
data("BinomialExample")
\donttest{
fit <- enetLTS(GaussianExample$x, GaussianExample$y)
assess.enetLTS(fit, GaussianExample$xtest, GaussianExample$ytest)

fitb <- enetLTS(BinomialExample$x, BinomialExample$y, family = "binomial")
assess.enetLTS(fitb, BinomialExample$xtest, BinomialExample$ytest)
confusion.enetLTS(fitb, BinomialExample$xtest, BinomialExample$ytest)
roc.enetLTS(fitb, BinomialExample$xtest, BinomialExample$ytest, plot = TRUE)
}
}
\author{
Fatma Sevinc Kurnaz, Irene Hoffmann, Peter Filzmoser
\cr Maintainer: Fatma Sevinc Kurnaz <fatmasevinckurnaz@gmail.com>
}
\keyword{models}
\keyword{classif}
