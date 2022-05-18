##  used for computation of the estimator of Bianco and Yohai (1996) in logistic regression

##  For more details we refer to
##     Croux, C., and Haesbroeck, G. (2003),
##     ``Implementing the Bianco and Yohai estimator for Logistic Regression'',
##     Computational Statistics and Data Analysis, 44, 273-295

dev1 <- function(s,y) log(1+exp(-abs(s))) + abs(s)*((y-0.5)*s<0)
dev2 <- function(s,y) log1p(exp(-abs(s))) + abs(s)*((y-0.5)*s<0)
dev3 <- function(s,y) -( y  * plogis(s, log.p=TRUE) +
                            (1-y)*plogis(s, lower.tail=FALSE, log.p=TRUE))


devBY <- dev1
rm(dev1, dev2, dev3)

GBY3Fs <- function(s,c3)
{
   e.f <- exp(0.25)*sqrt(pi)
   Fs <- exp(-devBY(s,1))
   resGinf <- e.f*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fs))))-1)
   resGinf <- (resGinf+(Fs*exp(-sqrt(-log(Fs)))))*as.numeric(s <= -log(exp(c3)-1))
   resGsup <- ((Fs*exp(-sqrt(c3)))+(e.f*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1))) *
      as.numeric(s > -log(exp(c3)-1))
   resGinf + resGsup
}

GBY3Fsm <- function(s,c3)
{
   e.f <- exp(0.25)*sqrt(pi)
   Fsm <- exp(-devBY(-s,1))
   resGinf <- e.f*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fsm))))-1)
   resGinf <- (resGinf+(Fsm*exp(-sqrt(-log(Fsm))))) * as.numeric(s >= log(exp(c3)-1))
   resGsup <- ((Fsm*exp(-sqrt(c3)))+(e.f*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1))) *
      as.numeric(s < log(exp(c3)-1))
   resGinf + resGsup
}

rhoBY3 <- function(t,c3)
{
   ec3 <- exp(-sqrt(c3))
   t*ec3* (t <= c3) +
      (ec3*(2+(2*sqrt(c3))+c3) - 2*exp(-sqrt(t))*(1+sqrt(t)))* (t > c3)
}


phiBY3 <- function(s,y,c3)
{
   s <- as.double(s)
   dev. <- devBY(s,y)
   rhoBY3(dev.,c3) + GBY3Fs(s,c3) + GBY3Fsm(s,c3)
}