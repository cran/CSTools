% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BEI_PDFBest.R
\name{PDFIndexHind}
\alias{PDFIndexHind}
\title{Computing the Index PDFs for a dataset of SFSs for a hindcats period.}
\usage{
PDFIndexHind(
  index_hind,
  index_obs,
  method = "ME",
  time_dim_name = "time",
  na.rm = FALSE
)
}
\arguments{
\item{index_hind}{Index (e.g. NAO index) array from SFSs
with at least two dimensions (time , member) or (time, statistic). 
The temporal dimension, by default 'time', must be greater than 2. 
The dimension 'member' must be greater than 1. 
The dimension 'statistic' must be equal to 2, for containing the two 
paramenters of a normal distribution (mean and sd) representing the ensemble 
of a SFS. It is not possible to have the dimension 'member' and  'statistic' 
together.}

\item{index_obs}{Index (e.g. NAO index) array from an observational database
or reanalysis with at least a temporal dimension (by default 'time'), 
which must be greater than 2.}

\item{method}{A character string indicating which methodology is applied
to compute the PDFs. One of "ME" (default) or "LMEV".}

\item{time_dim_name}{A character string indicating the name of the temporal 
dimension, by default 'time'.}

\item{na.rm}{Logical (default = FALSE). Should missing values be removed?}
}
\value{
An array with at least two dimensions (time, statistic = 4). The firt 
statistic is the parameter 'mean' of the PDF with not bias corrected.
The second statistic is the parameter 'standard deviation' of the PDF with not 
bias corrected. The third statistic is the parameter 'mean' of the PDF with 
bias corrected. The fourth statistic is the parameter 'standard deviation' of 
the PDF with bias corrected.
}
\description{
This function implements the computation to obtain the index PDFs
(e.g. NAO index) to improve the index estimate from SFSs for a hindcast period.
}
\examples{
# Example for the PDFIndexHind function
# Example 1 
index_obs <- 1 : (5 * 3 ) 
dim(index_obs) <- c(time = 5, season = 3)
index_hind <- 1 : (5 * 4 * 3)
dim(index_hind) <- c(time = 5, member = 4, season = 3)
res <- PDFIndexHind(index_hind, index_obs)
dim(res)
# time statistic  season
#   5         4        3
# Example 2
index_obs <- 1 : (5 * 3) 
dim(index_obs) <- c(time = 5, season = 3)
index_hind <- 1 : (5 * 2 * 3)
dim(index_hind) <- c(time = 5, statistic = 2, season = 3)
res <- PDFIndexHind(index_hind, index_obs)
}
\references{
Regionally improved seasonal forecast of precipitation through Best
estimation of winter NAO, Sanchez-Garcia, E. et al.,
Adv. Sci. Res., 16, 165174, 2019, \doi{10.5194/asr-16-165-2019}
}
\author{
Eroteida Sanchez-Garcia - AEMET, \email{esanchezg@aemet.es}
}
