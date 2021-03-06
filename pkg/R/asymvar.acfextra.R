## estimates long run variance by extrapolating the acf
## input:
# x: time series
# order.max: number of estimated autocorrelations which are used for extrapolation (if aic=TRUE maximal possible AR order)
# aic: should AR order be estimated based on aic criterium?
# obs: determines whether the asymptotical variance of the time series or of the ranks of the time series is estimated
## output:
# asy: estimated long run variance
# order: order of fitted AR model

asymvar.acfextra <- function(x, obs = c("untransformed", "ranks"), order.max = 2, aic = FALSE){
  n <- length(x)
  obs <- match.arg(obs)
  if (obs=="ranks") {
    x <- rank(x)/length(x)
	}
  armodel <- ar(x,order.max=order.max,aic=aic)
  if (length(armodel$ar)==0) acfest <- rep(0,n-1) else{
 	  acfest <- ARMAacf(ar=armodel$ar,lag.max=n-1)[-1]
  }
  w <- (n-1:(n-1))/n
  asy <- 1+2*sum(acfest*w)
  asy <- asy*var(x)
  erg <- list(lrv=asy,order=length(armodel$ar))	
  return(erg)
}
