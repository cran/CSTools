% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CST_CategoricalEnsCombination.R
\name{CategoricalEnsCombination}
\alias{CategoricalEnsCombination}
\title{Make categorical forecast based on a multi-model forecast with potential for 
calibrate}
\usage{
CategoricalEnsCombination(fc, obs, cat.method, eval.method, amt.cat, ...)
}
\arguments{
\item{fc}{A multi-dimensional array with named dimensions containing the 
seasonal forecast experiment data in the element named \code{$data}. The 
amount of forecasting models is equal to the size of the \code{dataset} 
dimension of the data array. The amount of members per model may be 
different. The  size of the \code{member} dimension of the data array is 
equal to the maximum of the ensemble members among the models. Models with 
smaller ensemble sizes have residual indices of \code{member} dimension in 
the data array filled with NA values.}

\item{obs}{A multidimensional array with named dimensions containing the 
observed data in the element named \code{$data}.}

\item{cat.method}{Method used to produce the categorical forecast, can be 
either \code{pool}, \code{comb}, \code{mmw} or \code{obs}. The method pool 
assumes equal weight for all ensemble members while the method comb assumes 
equal weight for each model. The weighting method is descirbed in 
Rajagopalan et al. (2002), Robertson et al. (2004) and Van Schaeybroeck and 
Vannitsem (2019). Finally, the \code{obs} method classifies the observations
into the different categories and therefore contains only 0 and 1 values.}

\item{eval.method}{Is the sampling method used, can be either 
\code{"in-sample"} or \code{"leave-one-out"}. Default value is the 
\code{"leave-one-out"} cross validation.}

\item{amt.cat}{Is the amount of categories. Equally-sized quantiles will be 
calculated based on the amount of categories.}

\item{...}{Other parameters to be passed on to the calibration procedure.}
}
\value{
An array containing the categorical forecasts in the element called 
\code{$data}. The first two dimensions of the returned object are named 
dataset and member and are both of size one. An additional dimension named 
category is introduced and is of size amt.cat.
}
\description{
This function converts a multi-model ensemble forecast into a 
categorical forecast by giving the probability for each category. Different 
methods are available to combine the different ensemble forecasting models 
into probabilistic categorical forecasts. 

See details in ?CST_CategoricalEnsCombination
}
\references{
Rajagopalan, B., Lall, U., & Zebiak, S. E. (2002). Categorical 
climate forecasts through regularization and optimal combination of multiple 
GCM ensembles. Monthly Weather Review, 130(7), 1792-1811.

Robertson, A. W., Lall, U., Zebiak, S. E., & Goddard, L. (2004). 
Improved combination of multiple atmospheric GCM ensembles for seasonal 
prediction. Monthly Weather Review, 132(12), 2732-2744.

Van Schaeybroeck, B., & Vannitsem, S. (2019). Postprocessing of 
Long-Range Forecasts. In Statistical Postprocessing of Ensemble Forecasts (pp. 267-290).
}
\author{
Bert Van Schaeybroeck, \email{bertvs@meteo.be}
}
