#### gives raw and reweighted

weights.enetLTS <-
   function(object,vers=c("reweighted","raw"),index=FALSE,...){

      vers <- match.arg(vers)

      if (vers=="reweighted"){
         wt <- object$wt
         if (index==TRUE){
             names(wt) <- 1:length(wt)
         }
         return(wt)
      } 
      if (vers=="raw"){
         raw.wt <- object$raw.wt
         if (index==TRUE){
             names(raw.wt) <- 1:length(raw.wt)
         }
         return(raw.wt)
      } 
   
   }
