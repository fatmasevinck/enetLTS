addColnames <- function(x) {
   # 'x' needs to be a matrix
   if(is.null(colnames(x))) colnames(x) <- paste("x", seq_len(ncol(x)), sep="")
   x
}

## add intercept column to design matrix
addIntercept <- function(x, check = FALSE) {
   if(!check || is.na(match("(Intercept)", colnames(x)))) {
      cbind("(Intercept)"=rep.int(1, nrow(x)), x)
   } else x
}


uptrimMSE<- function(x,trim=0.1){
   # computes trim% upper trimmed mean
   return(mean(x[x<quantile(x,1-trim)]))
}

uptrimCrit <- function(x,trim=0.1) {
   # return(mean(x[x<quantile(x,1-trim,na.rm=TRUE)],na.rm=TRUE))
   return(mean(sort(x)[1:(length(x)*(1-trim))],na.rm=TRUE))
}

## use in lambda00
winsorize.default <- function(x, standardized = FALSE, centerFun = median,
                              scaleFun = mad, const = 2,
                              return = c("data", "weights"), ...) {
   ## initializations
   standardized <- isTRUE(standardized)
   if(standardized) return <- match.arg(return)
   else {
      # standardize data
      x <- robStandardize(x, centerFun=centerFun, scaleFun=scaleFun, ...)
      center <- attr(x, "center")
      scale <- attr(x, "scale")
   }
   ## winsorize standardized data
   #   ind <- abs(x) > const           # observations in 'x' that need to be shrunken
   #   x[ind] <- const * sign(x[ind])  # winsorize
   weights <- pmin(const / abs(x), 1)
   if(standardized && return == "weights") return(weights)
   x <- weights * x
   ## finalizations
   if(!standardized) {
      # transform back to original scale and remove attributes
      x <- c(x * scale + center)
   }
   x
}

utils::globalVariables(c("Comp.1", "Comp.2", "Wavelength", "residuals"))

