% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CST_SaveExp.R
\name{CST_SaveExp}
\alias{CST_SaveExp}
\title{Save CSTools objects of class 's2dv_cube' containing experiments or observed 
data in NetCDF format}
\usage{
CST_SaveExp(data, destination = "./CST_Data")
}
\arguments{
\item{data}{an object of class \code{s2dv_cube}.}

\item{destination}{a character string containing the directory name in which 
to save the data. NetCDF file for each starting date are saved into the 
folder tree: destination/experiment/variable/. By default the function 
creates and saves the data into the folder "CST_Data" in the working 
directory.}
}
\description{
This function allows to divide and save a object of class 
's2dv_cube' into a NetCDF file, allowing to reload the saved data using 
\code{CST_Load} function.
}
\examples{
\dontrun{
library(CSTools)
data <- lonlat_data$exp
destination <- "./path2/"
CST_SaveExp(data = data, destination = destination)
}

}
\seealso{
\code{\link{CST_Load}}, \code{\link{as.s2dv_cube}} and \code{\link{s2dv_cube}}
}
\author{
Perez-Zanon Nuria, \email{nuria.perez@bsc.es}
}
