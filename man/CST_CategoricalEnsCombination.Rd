% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CST_CategoricalEnsCombination.R
\name{CST_CategoricalEnsCombination}
\alias{CST_CategoricalEnsCombination}
\title{Make categorical forecast based on a multi-model forecast with potential for calibrate}
\usage{
CST_CategoricalEnsCombination(
  exp,
  obs,
  cat.method = "pool",
  eval.method = "leave-one-out",
  amt.cat = 3,
  ...
)
}
\arguments{
\item{exp}{an object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the seasonal forecast experiment data in the element named \code{$data}. The amount of forecasting models is equal to the size of the \code{dataset} dimension of the data array. The amount of members per model may be different. The  size of the \code{member} dimension of the data array is equal to the maximum of the ensemble members among the models. Models with smaller ensemble sizes have residual indices of \code{member} dimension in the data array filled with NA values.}

\item{obs}{an object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the observed data in the element named \code{$data}.}

\item{cat.method}{method used to produce the categorical forecast, can be either \code{pool}, \code{comb}, \code{mmw} or \code{obs}. The method pool assumes equal weight for all ensemble members while the method comb assumes equal weight for each model. The weighting method is descirbed in Rajagopalan et al. (2002), Robertson et al. (2004) and Van Schaeybroeck and Vannitsem (2019). Finally, the \code{obs} method classifies the observations into the different categories and therefore contains only 0 and 1 values.}

\item{eval.method}{is the sampling method used, can be either \code{"in-sample"} or \code{"leave-one-out"}. Default value is the \code{"leave-one-out"} cross validation.}

\item{amt.cat}{is the amount of categories. Equally-sized quantiles will be calculated based on the amount of categories.}

\item{...}{other parameters to be passed on to the calibration procedure.}
}
\value{
an object of class \code{s2dv_cube} containing the categorical forecasts in the element called \code{$data}. The first two dimensions of the returned object are named dataset and member and are both of size one. An additional dimension named category is introduced and is of size amt.cat.
}
\description{
This function converts a multi-model ensemble forecast 
into a categorical forecast by giving the probability
for each category. Different methods are available to combine 
the different ensemble forecasting models into 
probabilistic categorical forecasts. 

Motivation:
Beyond the short range, the unpredictable component of weather 
predictions becomes substantial due to the chaotic nature of the earth 
system. Therefore, predictions can mostly be skillful when used in a probabilistic sense. 
In practice this is done using ensemble forecasts. It is then common to
convert the ensemble forecasts to occurence probabilities for different categories. 
These categories typically are taken as terciles from climatolgical distributions.
For instance for temperature, there is a cold, normal and warm class. 
Commonly multiple ensemble forecasting systems
are available but some models may be more competitive than others 
for the variable, region and user need under consideration. Therefore, 
when calculating the category probabilities, the ensemble members of 
the different forecasting system may be differently weighted.
Such weighting is typically done by comparison of the ensemble forecasts
with observations. 

Description of the tool:
The tool considers all forecasts (all members from all forecasting systems) 
and converts them into occurrence probabilities of different categories.
The amount of categories can be changed and are taken as the 
climatological quantiles (e.g. terciles), extracted 
from the observational data. 
The methods that are available to combine the ensemble forecasting models into 
probabilistic categorical forecasts are: 1) ensemble pooling where 
all ensemble members of all ensemble systems are weighted equally, 
2) model combination where each model system is weighted equally, and,
3) model weighting. 
The model weighting method is described in Rajagopalan et al. (2002),
Robertson et al. 2004 and Van Schaeybroeck and Vannitsem (2019). 
More specifically, this method uses different weights for the
occurence probability predicted by the available models and by a climatological model
and optimizes the weights by minimizing the ignorance score. 
Finally, the function can also be used to categorize the observations 
in the categorical quantiles.
}
\examples{

mod1 <- 1 : (2 * 3 * 4 * 5 * 6 * 7)
dim(mod1) <- c(dataset = 2, member = 3, sdate = 4, ftime = 5, lat = 6, lon = 7)
mod1[ 2, 3, , , , ] <- NA
dimnames(mod1)[[1]] <- c("MF", "UKMO")
obs1 <- 1 : (1 * 1 * 4 * 5 * 6 * 7)
dim(obs1) <- c(dataset = 1, member = 1, sdate = 4, ftime = 5, lat = 6, lon = 7)
lon <- seq(0, 30, 5)
lat <- seq(0, 25, 5)
exp <- list(data = mod1, lat = lat, lon = lon)
obs <- list(data = obs1, lat = lat, lon = lon)
attr(exp, 'class') <- 's2dv_cube'
attr(obs, 'class') <- 's2dv_cube'
\donttest{
a <- CST_CategoricalEnsCombination(exp = exp, obs = obs, amt.cat = 3, cat.method = "mmw")
}
}
\references{
Rajagopalan, B., Lall, U., & Zebiak, S. E. (2002). Categorical climate forecasts through regularization and optimal combination of multiple GCM ensembles. Monthly Weather Review, 130(7), 1792-1811.

Robertson, A. W., Lall, U., Zebiak, S. E., & Goddard, L. (2004). Improved combination of multiple atmospheric GCM ensembles for seasonal prediction. Monthly Weather Review, 132(12), 2732-2744.

Van Schaeybroeck, B., & Vannitsem, S. (2019). Postprocessing of Long-Range Forecasts. In Statistical Postprocessing of Ensemble Forecasts (pp. 267-290).
}
\author{
Bert Van Schaeybroeck, \email{bertvs@meteo.be}
}
