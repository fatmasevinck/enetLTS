plot.enetLTS <- function(x, method = c("coefficients", "resid", "diagnostic"),
                         vers = c("reweighted", "raw"), ...){

   vers <- match.arg(vers)
   if (missing(method)) {
      ## all three plots; in an interactive session ask before each new page
      if (interactive()) {
         old <- grDevices::devAskNewPage(TRUE)
         on.exit(grDevices::devAskNewPage(old))
      }
      plotCoef.enetLTS(x, vers = vers, ...)
      plotResid.enetLTS(x, vers = vers, ...)
      plotDiagnostic.enetLTS(x, vers = vers, ...)
   } else {
      method <- match.arg(method)
      switch(method,
             coefficients = plotCoef.enetLTS(x, vers = vers, ...),
             resid        = plotResid.enetLTS(x, vers = vers, ...),
             diagnostic   = plotDiagnostic.enetLTS(x, vers = vers, ...))
   }
   invisible(x)
}
