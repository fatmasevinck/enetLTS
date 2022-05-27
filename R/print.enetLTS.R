

print.enetLTS <-
   function(x,vers=c("reweighted","raw"),...){

      # require(enetLTS)
      # require(predict.enetLTS)
      # require(coef.enetLTS)

      vers <- match.arg(vers)
      cat("enetLTS estimator \n")

      cat("\nCall: ", deparse(x$call), "\n\n")
      
      if (x$inputs$family=="multinomial"){
         coefficients <- coef.enetLTS(x,vers=vers)
          cat("\nCoefficients:\n")
      print(coefficients)
      } else {
         coefficients <- drop(coef.enetLTS(x,vers=vers))
          cat("\nCoefficients:\n")
      print(unlist(coefficients))
      }
      
     

      nCoefficients <- sum(unlist(coefficients)!=0)
      cat("\n number of the nonzero coefficients:\n")
      print(nCoefficients)

      cat("\n alpha:",x$alpha)
      cat("\n lambda:",x$lambda)
      cat("\n lambdaw:",x$lambdaw)

}
