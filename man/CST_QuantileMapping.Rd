% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CST_QuantileMapping.R
\name{CST_QuantileMapping}
\alias{CST_QuantileMapping}
\title{Quantile Mapping for seasonal or decadal forecast data}
\usage{
CST_QuantileMapping(
  exp,
  obs,
  exp_cor = NULL,
  sdate_dim = "sdate",
  memb_dim = "member",
  window_dim = NULL,
  method = "QUANT",
  na.rm = FALSE,
  ncores = NULL,
  ...
)
}
\arguments{
\item{exp}{An object of class \code{s2dv_cube}.}

\item{obs}{An object of class \code{s2dv_cube}.}

\item{exp_cor}{An object of class \code{s2dv_cube} in which the quantile 
mapping correction should be applied. If it is not specified, the correction
is applied in object 'exp'.}

\item{sdate_dim}{A character string indicating the dimension name in which 
cross-validation would be applied when exp_cor is not provided. 'sdate' by 
default.}

\item{memb_dim}{A character string indicating the dimension name where
ensemble members are stored in the experimental arrays. It can be NULL if 
there is no ensemble member dimension. It is set as 'member' by default.}

\item{window_dim}{A character string indicating the dimension name where 
samples have been stored. It can be NULL (default) in case all samples are 
used.}

\item{method}{A character string indicating the method to be used:'PTF', 
'DIST', 'RQUANT', 'QUANT', 'SSPLIN'. By default, the empirical quantile 
mapping 'QUANT' is used.}

\item{na.rm}{A logical value indicating if missing values should be removed
(FALSE by default).}

\item{ncores}{An integer indicating the number of cores for parallel 
computation using multiApply function. The default value is NULL (1).}

\item{...}{Additional parameters to be used by the method choosen. See qmap 
package for details.}
}
\value{
An object of class \code{s2dv_cube} containing the experimental data
after applying the quantile mapping correction.
}
\description{
This function is a wrapper of fitQmap and doQmap from package
'qmap' to be applied on the object of class 's2dv_cube'. The quantile mapping
adjustment between an experiment, typically a hindcast, and observation is 
applied to the experiment itself or to a provided forecast.
}
\examples{
# Use synthetic data
exp <- NULL
exp$data <- 1 : c(1 * 3 * 5 * 4 * 3 * 2)
dim(exp$data) <- c(dataset = 1, member = 3, sdate = 5, ftime = 4,
                  lat = 3, lon = 2)
class(exp) <- 's2dv_cube'
obs <- NULL
obs$data <- 101 : c(100 + 1 * 1 * 5 * 4 * 3 * 2)
dim(obs$data) <- c(dataset = 1, member = 1, sdate = 5, ftime = 4,
                  lat = 3, lon = 2)
class(obs) <- 's2dv_cube'
res <- CST_QuantileMapping(exp, obs)

}
\seealso{
\code{\link[qmap]{fitQmap}} and \code{\link[qmap]{doQmap}}
}
\author{
Nuria Perez-Zanon, \email{nuria.perez@bsc.es}
}
