
##  F.S. Kurnaz

center.y <- function(y, index=NULL, centerFun = mean) {

   if (is.null(index)){
       center <- centerFun(y)  # compute center
   } else {
      center <- centerFun(y[index])
   }
   y <- y - center  # sweep out center
   # add attributes and return centered data
   attr(y, "center") <- center
   y
}

standardize.x <- function(x, index=NULL, centerFun = mean, scaleFun = sd) {

   # generic standardization with special cases mean/sd and median/MAD
   if (is.null(index)){
      if(identical(centerFun, mean)) {
         center <- colMeans(x)  # compute column means (faster than apply)
      } else center <- apply(x, 2, centerFun)  # compute column centers
      x <- sweep(x, 2, center, check.margin=FALSE)  # sweep out column centers

      if(identical(centerFun, mean) && identical(scaleFun, sd)) {
         # classical standardization with mean and standard deviation
         f <- function(v) sqrt(sum(v^2) / max(1, length(v)-1))
         scale <- apply(x, 2, f)
      } else if(identical(centerFun, median) && identical(scaleFun, mad)) {
         # robust standardization with median and MAD
         # compute column MADs with median already swept out
         scale <- apply(x, 2, mad, center=0)
      } else {
         scale <- apply(x, 2, scaleFun)}  # compute column scales
      x <- sweep(x, 2, scale, "/", check.margin=FALSE)  # sweep out column scales
   } else {
      if(identical(centerFun, mean)) {
         center <- colMeans(x[index,])  # compute column means (faster than apply)
      } else {
         center <- apply(x[index,], 2, centerFun)}  # compute column centers
      x <- sweep(x, 2, center, check.margin=FALSE)  # sweep out column centers

      if(identical(centerFun, mean) && identical(scaleFun, sd)) {
         # classical standardization with mean and standard deviation
         scale <- apply(x[index,], 2, sd)
      } else if(identical(centerFun, median) && identical(scaleFun, mad)) {
         # robust standardization with median and MAD
         # compute column MADs with median already swept out
         scale <- apply(x[index,], 2, mad, center=0)
      } else {
         scale <- apply(x[index,], 2, scaleFun)}  # compute column scales
      x <- sweep(x, 2, scale, "/", check.margin=FALSE)  # sweep out column scales
   }
   # add attributes and return standardized data
   attr(x, "center") <- center
   attr(x, "scale") <- scale
   x
}





