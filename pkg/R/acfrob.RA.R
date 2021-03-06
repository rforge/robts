#####################
# calculates the acf using a RA-approach
# input
# x: time series without NA as vector
# lag.max: maximal lag of interest
# psi: which Psi-Function should be used? ("huber" and "bisquare" are available)
# locfn: a mean estimator
# scalefn: an estimator of scatter 
# ...: Robustness-Parameter for Psi-Functions
# output: acf
#####################

acfrob.RA <- function(x, lag.max, psi = c("huber", "bisquare"), k, locfn = median, scalefn = mad, biascorr = TRUE, ...) {
  psi <- match.arg(psi)
  n <- length(x)
  lags <- 1:lag.max
  
  # protective measures:
  if (!is.function(locfn)) {
  	warning("This is no suitable location estimator. Median is used instead")
  	locfn <- median
 	}
  if (!is.function(scalefn)) {
  	warning("This is no suitable scatter estimator. MAD is used instead")
  	scalefn <- mad
 	}
  
  # calculation of mean and variance:
  locval <- try(locfn(x), silent=TRUE)
  if (inherits(locval, "try-error")) {
  	warning("Something went wrong with the mean estimation")
  	return(NA)
 	}
  scaleval <- try(scalefn(x, ...), silent=TRUE)
  if (inherits(scaleval, "try-error")) {
  	warning("Something went wrong with the scale estimation")
  	return(NA)
 	}
  
  # transformation of values:
  x_transformed <- M_psi((x-locval)/scaleval, type=psi, k=k)
  # calculation of acf (biased!): 
  acfvalues_biased <- acf(x_transformed, demean=FALSE, plot=FALSE, lag.max=lag.max)$acf[-1]
  if(!biascorr) return(acfvalues_biased)
  
  # transformation for unbiasedness:
  if (psi=="huber") {
    if(!missing(k) && k!=1.37) warning("The automatic bias correction is only valid for k=1.37")
    load(system.file("extdata", "acfbiascorr_RAhuber", package = "robts")) # loading the simulated expected values
  }
  if (psi=="bisquare") {
    if(!missing(k) && k!=4.68) warning("The automatic bias correction is only valid for k=4.68")
    load(system.file("extdata", "acfbiascorr_RAbisquare", package = "robts")) # loading the simulated expected values
  }
  if (psi %in% c("huber", "bisquare")) {
    acfvalues <- sapply(acfvalues_biased, linearinterpol, a=get("expectations"), b=get("values"))
  }else{
    warning("No bias correction has been made. The estimation could be slightly biased")
 	} 	
 	
 	res <- list(
   acfvalues = acfvalues,
   are = NA
  )
  	
  return(res)
}