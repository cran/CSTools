% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CST_WeatherRegimes.R
\name{CST_WeatherRegimes}
\alias{CST_WeatherRegimes}
\title{Function for Calculating the Cluster analysis}
\usage{
CST_WeatherRegimes(
  data,
  ncenters = NULL,
  EOFs = TRUE,
  neofs = 30,
  varThreshold = NULL,
  method = "kmeans",
  iter.max = 100,
  nstart = 30,
  ncores = NULL
)
}
\arguments{
\item{data}{An 's2dv_cube' object.}

\item{ncenters}{Number of clusters to be calculated with the clustering 
function.}

\item{EOFs}{Whether to compute the EOFs (default = 'TRUE') or not (FALSE) to 
filter the data.}

\item{neofs}{Number of modes to be kept (default = 30).}

\item{varThreshold}{Value with the percentage of variance to be explained by 
the PCs. Only sufficient PCs to explain this much variance will be used in 
the clustering.}

\item{method}{Different options to estimate the clusters. The most traditional 
approach is the k-means analysis (default=’kmeans’) but the function also 
support the different methods included in the hclust . These methods are:
"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" 
(= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). For more details 
about these methods see the hclust function documentation included in the 
stats package.}

\item{iter.max}{Parameter to select the maximum number of iterations allowed 
(Only if method='kmeans' is selected).}

\item{nstart}{Parameter for the cluster analysis determining how many random 
sets to choose (Only if method='kmeans' is selected).}

\item{ncores}{The number of multicore threads to use for parallel computation.}
}
\value{
A list with two elements \code{$data} (a 's2dv_cube' object containing 
the composites cluster = 1,..,K for case (*1) or only k = 1 for any specific 
cluster, i.e., case (*2)) and \code{$statistics} that includes \code{$pvalue} 
(array with the same structure as \code{$data} containing the pvalue of the 
composites obtained through a t-test that accounts for the serial dependence.),
\code{cluster} (A matrix or vector with integers (from 1:k) indicating the 
cluster to which each time step is allocated.), \code{persistence} (Percentage 
of days in a month/season before a cluster is replaced for a new one (only if 
method=’kmeans’ has been selected.)), \code{frequency} (Percentage of days in 
a month/season belonging to each cluster (only if method=’kmeans’ has been 
selected).),
}
\description{
This function computes the weather regimes from a cluster 
analysis. It is applied on the array \code{data} in a 's2dv_cube' object. The 
dimensionality of this object can be also reduced by using PCs obtained from 
the application of the #'EOFs analysis to filter the dataset. The cluster 
analysis can be performed with the traditional k-means or those methods
included in the hclust (stats package).
}
\examples{
data <- array(abs(rnorm(1280, 283.7, 6)), dim = c(dataset = 2, member = 2, 
                                                 sdate = 3, ftime = 3, 
                                                 lat = 4, lon = 4))
coords <- list(lon = seq(0, 3), lat = seq(47, 44))
obs <- list(data = data, coords = coords)
class(obs) <- 's2dv_cube'

res1 <- CST_WeatherRegimes(data = obs, EOFs = FALSE, ncenters = 4)
res2 <- CST_WeatherRegimes(data = obs, EOFs = TRUE, ncenters = 3)

}
\references{
Cortesi, N., V., Torralba, N., González-Reviriego, A., Soret, and 
F.J., Doblas-Reyes (2019). Characterization of European wind speed variability 
using weather regimes. Climate Dynamics,53, 4961–4976, 
\doi{10.1007/s00382-019-04839-5}.

Torralba, V. (2019) Seasonal climate prediction for the wind 
energy sector: methods and tools for the development of a climate service. 
Thesis. Available online: \url{https://eprints.ucm.es/56841/}.
}
\author{
Verónica Torralba - BSC, \email{veronica.torralba@bsc.es}
}
