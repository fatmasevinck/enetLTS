#### gives raw and reweighted or both weights from the object  

weights.enetLTS <-
   function(object,vers=c("reweighted","raw","both"),index=FALSE,...){

      vers <- match.arg(vers)

      if (vers=="reweighted"){
         wt <- object$wt
         if (index==TRUE){
             names(wt) <- 1:length(wt)
         }
         wfit <- list(wt=wt)
      } else if (vers=="raw"){
         raw.wt <- object$raw.wt
         if (index==TRUE){
             names(raw.wt) <- 1:length(raw.wt)
         }
         wfit <- list(raw.wt=raw.wt)
      } else if (vers=="both"){
         wt <- object$wt
         raw.wt <- object$raw.wt
         if (index==TRUE){
            names(wt) <- 1:length(wt)
            names(raw.wt) <- 1:length(raw.wt)
         }
         wfit <- list(wt=wt,raw.wt=raw.wt)
      }
      return(wfit)
   }
