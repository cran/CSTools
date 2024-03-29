% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CST_RFSlope.R
\name{RFSlope}
\alias{RFSlope}
\title{RainFARM spectral slopes from an array (reduced version)}
\usage{
RFSlope(
  data,
  kmin = 1,
  time_dim = NULL,
  lon_dim = "lon",
  lat_dim = "lat",
  ncores = NULL
)
}
\arguments{
\item{data}{Array containing the spatial precipitation fields to downscale.
The input array is expected to have at least two dimensions named "lon" and 
"lat" by default (these default names can be changed with the \code{lon_dim} 
and \code{lat_dim} parameters) and one or more dimensions over which to 
average the slopes, which can be specified by parameter \code{time_dim}.}

\item{kmin}{First wavenumber for spectral slope (default \code{kmin=1}).}

\item{time_dim}{String or character array with name(s) of dimension(s)
(e.g. "ftime", "sdate", "member" ...) over which to compute spectral slopes.
If a character array of dimension names is provided, the spectral slopes
will be computed as an average over all elements belonging to those dimensions.
If omitted one of c("ftime", "sdate", "time")  is searched and the first one
with more than one element is chosen.}

\item{lon_dim}{Name of lon dimension ("lon" by default).}

\item{lat_dim}{Name of lat dimension ("lat" by default).}

\item{ncores}{is an integer that indicates the number of cores for parallel 
computations using multiApply function. The default value is one.}
}
\value{
RFSlope() returns spectral slopes using the RainFARM convention
(the logarithmic slope of k*|A(k)|^2 where A(k) are the spectral amplitudes).
The returned array has the same dimensions as the input array,
minus the dimensions specified by \code{lon_dim}, \code{lat_dim} and \code{time_dim}.
}
\description{
This function computes spatial spectral slopes from an array,
to be used for RainFARM stochastic precipitation downscaling method.
}
\examples{
# Example for the 'reduced' RFSlope function 
# Create a test array with dimension 8x8 and 20 timesteps, 
# 3 starting dates and 20 ensemble members.
pr <- 1:(4*3*8*8*20)
dim(pr) <- c(ensemble = 4, sdate = 3, lon = 8, lat = 8, ftime = 20)
# Compute the spectral slopes ignoring the wavenumber
# corresponding to the largest scale (the box)
slopes <- RFSlope(pr, kmin = 2, time_dim = 'ftime')
}
\author{
Jost von Hardenberg - ISAC-CNR, \email{j.vonhardenberg@isac.cnr.it}
}
