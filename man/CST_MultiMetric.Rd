% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CST_MultiMetric.R
\name{CST_MultiMetric}
\alias{CST_MultiMetric}
\title{Multiple Metrics applied in Multiple Model Anomalies}
\usage{
CST_MultiMetric(
  exp,
  obs,
  metric = "correlation",
  multimodel = TRUE,
  time_dim = "ftime",
  memb_dim = "member",
  sdate_dim = "sdate"
)
}
\arguments{
\item{exp}{An object of class \code{s2dv_cube} as returned by 
\code{CST_Anomaly} function, containing the anomaly of the seasonal forecast 
experiments data in the element named \code{$data}.}

\item{obs}{An object of class \code{s2dv_cube} as returned by 
\code{CST_Anomaly} function, containing the anomaly of observed data in the 
element named \code{$data}.}

\item{metric}{A character string giving the metric for computing the maximum 
skill. This must be one of the strings 'correlation', 'rms', 'rmsss' and 
'rpss'. If 'rpss' is chossen the terciles probabilities are evaluated.}

\item{multimodel}{A logical value indicating whether a Multi-Model Mean should 
be computed.}

\item{time_dim}{Name of the temporal dimension where a mean will be applied. 
It can be NULL, the default value is 'ftime'.}

\item{memb_dim}{Name of the member dimension. It can be NULL, the default 
value is 'member'.}

\item{sdate_dim}{Name of the start date dimension or a dimension name 
identifiying the different forecast. It can be NULL, the default value is 
'sdate'.}
}
\value{
An object of class \code{s2dv_cube} containing the statistics of the 
selected metric in the element \code{$data} which is a list of arrays: for the
metric requested and others for statistics about its signeificance. The arrays 
have two dataset dimensions equal to the 'dataset' dimension in the 
\code{exp$data} and \code{obs$data} inputs. If \code{multimodel} is TRUE, the 
first position in the first 'nexp' dimension correspons to the Multi-Model Mean.
}
\description{
This function calculates correlation (Anomaly Correlation 
Coefficient; ACC), root mean square error (RMS) and the root mean square error 
skill score (RMSSS) of individual anomaly models and multi-models mean (if 
desired) with the observations.
}
\examples{
mod <- rnorm(2*2*4*5*2*2)
dim(mod) <- c(dataset = 2, member = 2, sdate = 4, ftime = 5, lat = 2, lon = 2)
obs <- rnorm(1*1*4*5*2*2)
dim(obs) <- c(dataset = 1, member = 1, sdate = 4, ftime = 5, lat = 2, lon = 2)
lon <- seq(0, 30, 5)
lat <- seq(0, 25, 5)
coords <- list(lat = lat, lon = lon)
exp <- list(data = mod, coords = coords)
obs <- list(data = obs, coords = coords)
attr(exp, 'class') <- 's2dv_cube'
attr(obs, 'class') <- 's2dv_cube'
a <- CST_MultiMetric(exp = exp, obs = obs)
}
\references{
Mishra, N., Prodhomme, C., & Guemas, V. (n.d.). Multi-Model Skill 
Assessment of Seasonal Temperature and Precipitation Forecasts over Europe, 
29-31. \doi{10.1007/s00382-018-4404-z}
}
\seealso{
\code{\link[s2dv]{Corr}}, \code{\link[s2dv]{RMS}}, 
\code{\link[s2dv]{RMSSS}} and \code{\link{CST_Load}}
}
\author{
Mishra Niti, \email{niti.mishra@bsc.es}

Perez-Zanon Nuria, \email{nuria.perez@bsc.es}
}
