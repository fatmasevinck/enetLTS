
print.enetLTS <-
   function(x,vers=c("reweighted","raw"),zeros=FALSE,...){

      # require(enetLTS)
      # require(predict.enetLTS)
      # require(coef.enetLTS)

      vers <- match.arg(vers)
      cat("enetLTS estimator \n")

      cat("\nCall: ", deparse(x$call), "\n\n")

      coefficients <- coef.enetLTS(x,vers=vers,zeros=zeros)
      cat("\nCoefficients:\n")
      print(coefficients,...)

      nCoefficients <- sum(coefficients!=0)
      cat("\n number of the nonzero coefficients:\n")
      print(nCoefficients)

      cat(paste("\n alpha:",x$alpha))
      cat(paste("\n lambda:",x$lambda))
      cat(paste("\n lambdaw:",x$lambdaw))

}
