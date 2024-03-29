% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CST_RFWeights.R
\name{CST_RFWeights}
\alias{CST_RFWeights}
\title{Compute climatological weights for RainFARM stochastic precipitation downscaling}
\usage{
CST_RFWeights(
  climfile,
  nf,
  lon,
  lat,
  varname = NULL,
  fsmooth = TRUE,
  lonname = "lon",
  latname = "lat",
  ncores = NULL
)
}
\arguments{
\item{climfile}{Filename of a fine-scale precipitation climatology. The file 
is expected to be in NetCDF format and should contain at least one 
precipitation field. If several fields at different times are provided,
a climatology is derived by time averaging. Suitable climatology files could 
be for example a fine-scale precipitation climatology from a high-resolution 
regional climate model (see e.g. Terzago et al. 2018), a local 
high-resolution gridded climatology from observations, or a reconstruction 
such as those which can be downloaded from the WORLDCLIM 
(\url{https://www.worldclim.org}) or CHELSA (\url{https://chelsa-climate.org/}) 
websites. The latter data will need to be converted to NetCDF format before 
being used (see for example the GDAL tools (\url{https://gdal.org/}). It  
could also be an 's2dv_cube' object.}

\item{nf}{Refinement factor for downscaling (the output resolution is 
increased by this factor).}

\item{lon}{Vector of longitudes.}

\item{lat}{Vector of latitudes. The number of longitudes and latitudes is 
expected to be even and the same. If not the function will perform a 
subsetting to ensure this condition.}

\item{varname}{Name of the variable to be read from \code{climfile}.}

\item{fsmooth}{Logical to use smooth conservation (default) or large-scale 
box-average conservation.}

\item{lonname}{A character string indicating the name of the longitudinal 
dimension set as 'lon' by default.}

\item{latname}{A character string indicating the name of the latitudinal 
dimension set as 'lat' by default.}

\item{ncores}{An integer that indicates the number of cores for parallel 
computations using multiApply function. The default value is one.}
}
\value{
An object of class 's2dv_cube' containing in matrix \code{data} the 
weights with dimensions (lon, lat).
}
\description{
Compute climatological ("orographic") weights from a fine-scale 
precipitation climatology file.
}
\examples{
# Create weights to be used with the CST_RainFARM() or RainFARM() functions
# using an external random data in the form of 's2dv_cube'.
obs <- rnorm(2 * 3 * 4 * 8 * 8)
dim(obs) <- c(dataset = 1, member = 2, sdate = 3, ftime = 4, lat = 8, lon = 8)
lon <- seq(10, 13.5, 0.5)
lat <- seq(40, 43.5, 0.5)
coords <- list(lon = lon, lat = lat)
data <- list(data = obs, coords = coords)
class(data) <- "s2dv_cube"
res <- CST_RFWeights(climfile = data, nf = 3, lon, lat, lonname = 'lon', 
                    latname = 'lat', fsmooth = TRUE)
}
\references{
Terzago, S., Palazzi, E., & von Hardenberg, J. (2018).
Stochastic downscaling of precipitation in complex orography: 
A simple method to reproduce a realistic fine-scale climatology.
Natural Hazards and Earth System Sciences, 18(11),
2825-2840. \doi{10.5194/nhess-18-2825-2018}.
}
\author{
Jost von Hardenberg - ISAC-CNR, \email{j.vonhardenberg@isac.cnr.it}
}
