% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CST_BiasCorrection.R
\name{CST_BiasCorrection}
\alias{CST_BiasCorrection}
\title{Bias Correction based on the mean and standard deviation adjustment}
\usage{
CST_BiasCorrection(
  exp,
  obs,
  exp_cor = NULL,
  na.rm = FALSE,
  memb_dim = "member",
  sdate_dim = "sdate",
  dat_dim = NULL,
  ncores = NULL
)
}
\arguments{
\item{exp}{An object of class \code{s2dv_cube} as returned by \code{CST_Start} 
function, containing the seasonal forecast experiment data in the element 
named \code{$data} with at least time and member dimensions.}

\item{obs}{An object of class \code{s2dv_cube} as returned by \code{CST_Start} 
function, containing the observed data in the element named \code{$data} 
with at least time dimension.}

\item{exp_cor}{An object of class \code{s2dv_cube} as returned by 
\code{CST_Start} function, containing the seasonal forecast experiment to be 
corrected with at least time dimension. If it is NULL, the 'exp' forecast 
will be corrected. If there is only one corrected dataset, it should not  
have dataset dimension. If there is a corresponding corrected dataset for  
each 'exp' forecast, the dataset dimension must have the same length as in 
'exp'. The default value is NULL.}

\item{na.rm}{A logical value indicating whether missing values should be 
stripped before the computation proceeds, by default it is set to FALSE.}

\item{memb_dim}{A character string indicating the name of the member 
dimension. By default, it is set to 'member'.}

\item{sdate_dim}{A character string indicating the name of the start date 
dimension. By default, it is set to 'sdate'.}

\item{dat_dim}{A character string indicating the name of dataset dimension. 
The length of this dimension can be different between 'exp' and 'obs'. 
The default value is NULL.}

\item{ncores}{An integer that indicates the number of cores for parallel 
computations using multiApply function. The default value is NULL.}
}
\value{
An object of class \code{s2dv_cube} containing the bias corrected 
forecasts with the dimensions nexp, nobs and same dimensions as in the 'exp' 
object. nexp is the number of experiment (i.e., 'dat_dim' in exp), and nobs is 
the number of observation (i.e., 'dat_dim' in obs). If dat_dim is NULL, nexp 
and nobs are omitted. If 'exp_cor' is provided the returned array will be with 
the same dimensions as 'exp_cor'.
}
\description{
This function applies the simple bias adjustment technique 
described in Torralba et al. (2017). The adjusted forecasts have an equivalent 
standard deviation and mean to that of the reference dataset.
}
\examples{
mod1 <- 1 : (1 * 3 * 4 * 5 * 6 * 7)
dim(mod1) <- c(dataset = 1, member = 3, sdate = 4, time = 5, lat = 6, lon = 7)
obs1 <- 1 : (1 * 1 * 4 * 5 * 6 * 7)
dim(obs1) <- c(dataset = 1, member = 1, sdate = 4, time = 5, lat = 6, lon = 7)
lon <- seq(0, 30, 5)
lat <- seq(0, 25, 5)
coords <- list(lat = lat, lon = lon)
exp <- list(data = mod1, coords = coords)
obs <- list(data = obs1, coords = coords)
attr(exp, 'class') <- 's2dv_cube'
attr(obs, 'class') <- 's2dv_cube'
a <- CST_BiasCorrection(exp = exp, obs = obs)
}
\references{
Torralba, V., F.J. Doblas-Reyes, D. MacLeod, I. Christel and M. 
Davis (2017). Seasonal climate prediction: a new source of information for 
the management of wind energy resources. Journal of Applied Meteorology and 
Climatology, 56, 1231-1247, \doi{10.1175/JAMC-D-16-0204.1}. (CLIM4ENERGY, 
EUPORIAS, NEWA, RESILIENCE, SPECS)
}
\author{
Verónica Torralba, \email{veronica.torralba@bsc.es}
}
