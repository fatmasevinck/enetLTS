

plot.enetLTS <- function(x,method=c("coefficients","resid","diagnostic"),
                         vers=c("reweighted","raw"),...){
   #  require(ggplot2)
   
   vers <- match.arg(vers)
   par(ask=T)
   readline(prompt = "Press <Enter> to continue...")
   ## call plot function
 
      if (missing("method")){
         plotCoef.enetLTS(x,vers=vers,...)
         plotResid.enetLTS(x,vers=vers,...)
         plotDiagnostic.enetLTS(x,vers=vers,...) }
      else if(method == "coefficients") plotCoef.enetLTS(x,vers=vers,...)
      else if(method == "resid") plotResid.enetLTS(x,vers=vers,...)
      else if(method == "diagnostic") plotDiagnostic.enetLTS(x,vers=vers,...)
   
}

