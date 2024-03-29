% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CST_DynBiasCorrection.R
\name{DynBiasCorrection}
\alias{DynBiasCorrection}
\title{Performing a Bias Correction conditioned by the dynamical
properties of the data.}
\usage{
DynBiasCorrection(
  exp,
  obs,
  method = "QUANT",
  wetday = FALSE,
  proxy = "dim",
  quanti,
  ncores = NULL
)
}
\arguments{
\item{exp}{A multidimensional array with named dimensions with the 
experiment data.}

\item{obs}{A multidimensional array with named dimensions with the 
observation data.}

\item{method}{A character string indicating the method to apply bias 
correction among these ones:
"PTF", "RQUANT", "QUANT", "SSPLIN".}

\item{wetday}{Logical indicating whether to perform wet day correction 
or not OR a numeric threshold below which all values are set to zero (by 
default is set to 'FALSE').}

\item{proxy}{A character string indicating the proxy for local dimension 
'dim' or inverse of persistence 'theta' to apply the dynamical conditioned 
bias correction method.}

\item{quanti}{A number lower than 1 indicating the quantile to perform the 
computation of local dimension and theta.}

\item{ncores}{The number of cores to use in parallel computation.}
}
\value{
A multidimensional array with named dimensions with a bias correction 
performed conditioned by local dimension 'dim' or inverse of persistence 'theta'.
}
\description{
This function perform a bias correction conditioned by the 
dynamical properties of the dataset. This function used the functions 
'CST_Predictability' to divide in terciles the two dynamical proxies 
computed with 'CST_ProxiesAttractor'. A bias correction
between the model and the observations is performed using the division into
terciles of the local dimension 'dim' and inverse of the persistence 'theta'.
For instance, model values with lower 'dim' will be corrected with observed 
values with lower 'dim', and the same for theta. The function gives two options
of bias correction: one for 'dim' and/or one for 'theta'
}
\examples{
expL <- rnorm(1:2000)
dim (expL) <- c(time =100,lat = 4, lon = 5)
obsL <- c(rnorm(1:1980),expL[1,,]*1.2)
dim (obsL) <- c(time = 100,lat = 4, lon = 5)
dynbias <- DynBiasCorrection(exp = expL, obs = obsL, method='QUANT',
                            proxy= "dim", quanti = 0.6)
}
\references{
Faranda, D., Alvarez-Castro, M.C., Messori, G., Rodriguez, D., 
and Yiou, P. (2019). The hammam effect or how a warm ocean enhances large 
scale atmospheric predictability.Nature Communications, 10(1), 1316. 
\doi{10.1038/s41467-019-09305-8}"

Faranda, D., Gabriele Messori and Pascal Yiou. (2017).
Dynamical proxies of North Atlantic predictability and extremes. 
Scientific Reports, 7-41278, 2017.
}
\author{
Carmen Alvarez-Castro, \email{carmen.alvarez-castro@cmcc.it}

Maria M. Chaves-Montero, \email{mdm.chaves-montero@cmcc.it}

Veronica Torralba, \email{veronica.torralba@cmcc.it}

Davide Faranda, \email{davide.faranda@lsce.ipsl.fr}
}
