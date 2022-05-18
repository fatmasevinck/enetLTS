
### weights for multinomial

  # weight function for Mahalanobis distances:
weight.multinomial <- function(md,c1=2.5,c2=5){
    w <- (1 - ((md - c1)/(c2 - c1))^2)^2
    w[md<c1] <- 1
    w[md>c2] <- 0
    return(w)
  } 

### weights for binomial
weight.binomial <- function(x, y, beta, intercept, del){
   if(intercept==TRUE){
      pi <- exp(cbind(1,x)%*%beta)/(1+exp(cbind(1,x)%*%beta))
      res <- (y - pi) / sqrt(pi*(1-pi))     ### pearson resiuals
   } else{
      pi <- exp(x%*%beta[-1])/(1+exp(x%*%beta[-1]))
      res <- (y - pi) / sqrt(pi*(1-pi))
   }
   we <- as.integer(abs(res) <= qnorm(1-del))
   return(we)
}

### weights for gaussian
weight.gaussian <- function(resi, ind, del){
   if(is.logical(ind)){
      h <- length(which(ind==TRUE))
   }else{
      h <- length(ind)
   }
   n <- length(resi)
   mu <- mean(resi[ind])
   rc <- (resi - mu)
   qn <- qnorm((h+n)/ (2*n))                         # required quantile
   cdelta <- 1 / sqrt(1 - (2*n)/(h/qn) * dnorm(qn))
   s <- sqrt(mean(rc[ind]^2)) * cdelta
   we <- as.integer(abs(rc/s) <= qnorm(1-del))
   out <- list(we=we,mu=mu,s=s)
   
   return(out)
}

