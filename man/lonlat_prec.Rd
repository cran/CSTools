% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\docType{data}
\name{lonlat_prec}
\alias{lonlat_prec}
\title{Sample Of Experimental Precipitation Data In Function Of Longitudes And Latitudes}
\description{
This sample data set contains a small cutout of gridded seasonal precipitation forecast data from the Copernicus Climate Change ECMWF-System 5 forecast system, to be used to demonstrate downscaling. Specifically, for the 'pr' (precipitation) variable, for the first 6 forecast ensemble members, daily values, for all 31 days in March following the forecast starting dates in November of years 2010 to 2012, for a small 4x4 pixel cutout in a region in the North-Western Italian Alps (44N-47N, 6E-9E). The data resolution is 1 degree.
}
\details{
The `CST_Load` call used to generate the data set in the infrastructure of the Marconi machine at CINECA is shown next, working on files which were extracted from forecast data available in the MEDSCOPE internal archive.

\preformatted{
library(CSTools)

infile <- list(path = '../medscope/nwalps/data/$VAR_NAME$_$START_DATE$_nwalps.nc')
lonlat_prec <- CST_Load('prlr', exp = list(infile), obs = NULL,
                        sdates = c('20101101', '20111101', '20121101'),
                        leadtimemin = 121, leadtimemax = 151,
                        latmin = 44, latmax = 47,
                        lonmin = 5, lonmax = 9,
                        nmember = 25,
                        storefreq = "daily", sampleperiod = 1,
                        output = "lonlat"
                       )$exp
}
}
\author{
Jost von Hardenberg \email{j.vonhardenberg@isac.cnr.it}
}
\keyword{data}
