% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CST_EnsClustering.R
\name{EnsClustering}
\alias{EnsClustering}
\title{Ensemble clustering}
\usage{
EnsClustering(
  data,
  lat,
  lon,
  time_moment = "mean",
  numclus = NULL,
  lon_lim = NULL,
  lat_lim = NULL,
  variance_explained = 80,
  numpcs = NULL,
  time_percentile = 90,
  time_dim = NULL,
  cluster_dim = "member",
  verbose = T
)
}
\arguments{
\item{data}{A matrix of dimensions 'dataset member sdate ftime lat lon' 
containing the variables to be analysed. Latitudinal dimension accepted 
names: 'lat', 'lats', 'latitude', 'y', 'j', 'nav_lat'. Longitudinal 
dimension accepted names: 'lon', 'lons','longitude', 'x', 'i', 'nav_lon'.}

\item{lat}{Vector of latitudes.}

\item{lon}{Vector of longitudes.}

\item{time_moment}{Decides the moment to be applied to the time dimension. Can 
be either 'mean' (time mean), 'sd' (standard deviation along time) or 'perc' 
(a selected percentile on time). If 'perc' the keyword 'time_percentile' is 
also used.}

\item{numclus}{Number of clusters (scenarios) to be calculated. If set to NULL
the number of ensemble members divided by 10 is used, with a minimum of 2 
and a maximum of 8.}

\item{lon_lim}{List with the two longitude margins in `c(-180,180)` format.}

\item{lat_lim}{List with the two latitude margins.}

\item{variance_explained}{variance (percentage) to be explained by the set of 
EOFs. Defaults to 80. Not used if numpcs is specified.}

\item{numpcs}{Number of EOFs retained in the analysis (optional).}

\item{time_percentile}{Set the percentile in time you want to analyse (used 
for `time_moment = "perc").}

\item{time_dim}{String or character array with name(s) of dimension(s) over 
which to compute statistics. If omitted c("ftime", "sdate", "time") are 
searched in this order.}

\item{cluster_dim}{Dimension along which to cluster. Typically "member" or 
"sdate". This can also be a list like c("member", "sdate").}

\item{verbose}{Logical for verbose output}
}
\value{
A list with elements \code{$cluster} (cluster assigned for each member),
\code{$freq} (relative frequency of each cluster), \code{$closest_member}
(representative member for each cluster), \code{$repr_field} (list of fields for
each representative member), \code{composites} (list of mean fields for each 
cluster), \code{$lon} (selected longitudes of output fields), \code{$lat} 
(selected longitudes of output fields).
}
\description{
This function performs a clustering on members/starting dates
and returns a number of scenarios, with representative members for each of 
them. The clustering is performed in a reduced EOF space.
}
\examples{
exp <- array(abs(rnorm(1152))*275, dim = c(dataset = 1, member = 4, 
                                          sdate = 6, ftime = 3, 
                                          lat = 4, lon = 4))
lon <- seq(0, 3)
lat <- seq(48, 45)
res <- EnsClustering(exp, lat = lat, lon = lon, numclus = 2,
                    cluster_dim = c("member", "dataset", "sdate"))

}
\author{
Federico Fabiano - ISAC-CNR, \email{f.fabiano@isac.cnr.it}

Ignazio Giuntoli - ISAC-CNR, \email{i.giuntoli@isac.cnr.it}

Danila Volpi - ISAC-CNR, \email{d.volpi@isac.cnr.it}

Paolo Davini - ISAC-CNR, \email{p.davini@isac.cnr.it}

Jost von Hardenberg - ISAC-CNR, \email{j.vonhardenberg@isac.cnr.it}
}
