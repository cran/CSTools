% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotWeeklyClim.R
\name{PlotWeeklyClim}
\alias{PlotWeeklyClim}
\title{Plots the observed weekly means and climatology of a timeseries data}
\usage{
PlotWeeklyClim(
  data,
  first_date,
  ref_period,
  last_date = NULL,
  data_years = NULL,
  time_dim = "time",
  sdate_dim = "sdate",
  ylim = NULL,
  title = NULL,
  subtitle = NULL,
  ytitle = NULL,
  legend = TRUE,
  palette = "Blues",
  fileout = NULL,
  device = NULL,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
}
\arguments{
\item{data}{A multidimensional array with named dimensions with at least sdate 
and time dimensions containing observed daily data. It can also be a 
dataframe with computed percentiles as input for ggplot. If it's a 
dataframe, it must contain the following column names: 'week', 'clim', 
'p10', 'p90', 'p33', 'p66', 'week_mean', 'day' and 'data'.}

\item{first_date}{The first date of the observed values of timeseries. It can 
be of class 'Date', 'POSIXct' or a character string in the format 
'yyyy-mm-dd'. If parameter 'data_years' is not provided, it must be a date
included in the reference period.}

\item{ref_period}{A vector of numeric values indicating the years of the 
reference period. If parameter 'data_years' is not specified, it must 
be of the same length of dimension 'sdate_dim' of parameter 'data'.}

\item{last_date}{Optional parameter indicating the last date of the target 
period of the daily timeseries. It can be of class 'Date', 'POSIXct' or a 
character string in the format 'yyyy-mm-dd'. If it is NULL, the last date of 
the daily timeseries will be set as the last date of 'data'. As the data is 
plotted by weeks, only full groups of 7 days will be plotted. If the last 
date of the timeseries is not a multiple of 7 days, the last week will 
not be plotted.}

\item{data_years}{A vector of numeric values indicating the years of the 
data. It must be of the same length of dimension 'sdate_dim' of parameter 
'data'. It is optional, if not specified, all the years will be used as the 
target period.}

\item{time_dim}{A character string indicating the daily time dimension name. 
The default value is 'time'.}

\item{sdate_dim}{A character string indicating the start year dimension name. 
The default value is 'sdate'.}

\item{ylim}{A numeric vector of length two providing limits of the scale. 
Use NA to refer to the existing minimum or maximum. For more information, 
see 'ggplot2' documentation of 'scale_y_continuous' parameter.}

\item{title}{The text for the top title of the plot. It is NULL by default.}

\item{subtitle}{The text for the subtitle of the plot. It is NULL bu default.}

\item{ytitle}{Character string to be drawn as y-axis title. It is NULL by 
default.}

\item{legend}{A logical value indicating whether a legend should be included 
in the plot. If it is TRUE or NA, the legend will be included. If it is 
FALSE, the legend will not be included. It is TRUE by default.}

\item{palette}{A palette name from the R Color Brewer’s package. The default 
value is 'Blues'.}

\item{fileout}{A character string indicating the file name where to save the 
plot. If not specified (default) a graphics device will pop up.}

\item{device}{A character string indicating the device to use. Can either be 
a device function (e.g. png), or one of "eps", "ps", "tex" (pictex), 
"pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only).}

\item{width}{A numeric value of the plot width in units ("in", "cm", "mm", or 
"px"). It is set to 8 by default.}

\item{height}{A numeric value of the plot height in units ("in", "cm", "mm", 
or "px"). It is set to 6 by default.}

\item{units}{Units of the size of the device (file or window) to plot in. 
Inches (’in’) by default.}

\item{dpi}{A numeric value of the plot resolution. It is set to 300 by 
default.}
}
\value{
A ggplot object containing the plot.
}
\description{
This function plots the observed weekly means and climatology of 
a timeseries data using ggplot package. It compares the weekly climatology in 
a specified period (reference period) to the observed conditions during the 
target period analyzed in the case study.
}
\examples{
data <- array(rnorm(49*20*3, 274), dim = c(time = 49, sdate = 20, member = 3))
PlotWeeklyClim(data = data, first_date = '2002-08-09', 
              last_date = '2002-09-15', ref_period = 2010:2019, 
              data_years = 2000:2019, time_dim = 'time', sdate_dim = 'sdate',
              title = "Observed weekly means and climatology", 
              subtitle = "Target years: 2010 to 2019", 
              ytitle = paste0('tas', " (", "deg.C", ")"))

}
