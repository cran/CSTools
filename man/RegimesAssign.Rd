% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CST_RegimesAssign.R
\name{RegimesAssign}
\alias{RegimesAssign}
\title{Function for matching a field of anomalies with 
a set of maps used as a reference (e.g. clusters obtained from the WeatherRegime function).}
\usage{
RegimesAssign(
  data,
  ref_maps,
  lat,
  method = "distance",
  composite = FALSE,
  memb = FALSE,
  ncores = NULL
)
}
\arguments{
\item{data}{an array containing anomalies with named dimensions: dataset, member, sdate, ftime, lat and lon.}

\item{ref_maps}{array with 3-dimensions ('lon', 'lat', 'cluster') containing the maps/clusters that will be used as a reference for the matching.}

\item{lat}{a vector of latitudes corresponding to the positions provided in data and ref_maps.}

\item{method}{whether the matching will be performed in terms of minimum distance (default = ’distance’) or 
the maximum spatial correlation (method=’ACC’) between the maps.}

\item{composite}{a logical parameter indicating if the composite maps are computed or not (default=FALSE).}

\item{memb}{a logical value indicating whether to compute composites for separate members (default FALSE) or as unique ensemble (TRUE).
This option is only available for when parameter 'composite' is set to TRUE and the data object has a dimension named 'member'.}

\item{ncores}{the number of multicore threads to use for parallel computation.}
}
\value{
A list with elements \code{$composite} (3-d array (lon, lat, k) containing the composites k=1,..,K for case (*1)
        \code{$pvalue} ( array with the same structure as \code{$composite} containing the pvalue of the composites obtained through a t-test 
        that accounts for the serial dependence of the data with the same structure as Composite.) (only if composite='TRUE'),
        \code{$cluster} (array with the same dimensions as data (except latitude and longitude which are removed) indicating the ref_maps to which each point is allocated.) ,
        \code{$frequency} (A vector of integers (from k = 1, ... k n reference maps) indicating the percentage of assignations corresponding to each map.),
}
\description{
This function performs the matching between a field of anomalies and a set
of maps which will be used as a reference. The anomalies will be assigned to the reference map 
for which the minimum Eucledian distance (method=’distance’) or highest spatial correlation 
(method=‘ACC’) is obtained.
}
\examples{
\dontrun{
regimes <- WeatherRegime(data = lonlat_data$obs$data, lat = lonlat_data$obs$lat,
                        EOFs = FALSE, ncenters = 4)$composite
res1 <- RegimesAssign(data = lonlat_data$exp$data, ref_maps = drop(regimes), 
                     lat = lonlat_data$exp$lat, composite = FALSE)
}	
}
\references{
Torralba, V. (2019) Seasonal climate prediction for the wind energy sector: methods and tools for the development of a climate service. Thesis. Available online: \url{https://eprints.ucm.es/56841/}
}
\author{
Verónica Torralba - BSC, \email{veronica.torralba@bsc.es}
}