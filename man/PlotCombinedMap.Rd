% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotCombinedMap.R
\name{PlotCombinedMap}
\alias{PlotCombinedMap}
\title{Plot Multiple Lon-Lat Variables In a Single Map According to a Decision Function}
\usage{
PlotCombinedMap(
  maps,
  lon,
  lat,
  map_select_fun,
  display_range,
  map_dim = "map",
  brks = NULL,
  cols = NULL,
  bar_limits = NULL,
  triangle_ends = c(F, F),
  col_inf = NULL,
  col_sup = NULL,
  col_unknown_map = "white",
  mask = NULL,
  col_mask = "grey",
  dots = NULL,
  bar_titles = NULL,
  legend_scale = 1,
  cex_bar_titles = 1.5,
  plot_margin = NULL,
  bar_extra_margin = c(2, 0, 2, 0),
  fileout = NULL,
  width = 8,
  height = 5,
  size_units = "in",
  res = 100,
  drawleg = T,
  return_leg = FALSE,
  ...
)
}
\arguments{
\item{maps}{List of matrices to plot, each with (longitude, latitude) 
dimensions, or 3-dimensional array with the dimensions (longitude, latitude, 
map). Dimension names are required.}

\item{lon}{Vector of longitudes. Must match the length of the corresponding 
dimension in 'maps'.}

\item{lat}{Vector of latitudes. Must match the length of the corresponding 
dimension in 'maps'.}

\item{map_select_fun}{Function that selects, for each grid point, which value 
to take among all the provided maps. This function receives as input a 
vector of values for a same grid point for all the provided maps, and must 
return a single selected value (not its index!) or NA. For example, the 
\code{min} and \code{max} functions are accepted.}

\item{display_range}{Range of values to be displayed for all the maps. This 
must be a numeric vector c(range min, range max). The values in the 
parameter 'maps' can go beyond the limits specified in this range. If the 
selected value for a given grid point (according to 'map_select_fun') falls 
outside the range, it will be coloured with 'col_unknown_map'.}

\item{map_dim}{Optional name for the dimension of 'maps' along which the 
multiple maps are arranged. Only applies when 'maps' is provided as a 
3-dimensional array. Takes the value 'map' by default.}

\item{brks}{Colour levels to be sent to PlotEquiMap. This parameter is 
optional and adjusted automatically by the function.}

\item{cols}{List of vectors of colours to be sent to PlotEquiMap for the 
colour bar of each map. This parameter is optional and adjusted 
automatically by the function (up to 5 maps). The colours provided for each 
colour bar will be automatically interpolated to match the number of breaks. 
Each item in this list can be named, and the name will be used as title for 
the corresponding colour bar (equivalent to the parameter 'bar_titles').}

\item{bar_limits}{Parameter from s2dv::ColorBar. Vector of two numeric values 
with the extremes of the range of values represented in the colour bar. If 
'var_limits' go beyond this interval, the drawing of triangle extremes is 
triggered at the corresponding sides, painted in 'col_inf' and 'col_sup'. 
Either of them can be set as NA and will then take as value the 
corresponding extreme in 'var_limits' (hence a triangle end won't be 
triggered for these sides). Takes as default the extremes of 'brks' if 
available, else the same values as 'var_limits'.}

\item{triangle_ends}{Parameter from s2dv::ColorBar. Vector of two logical 
elements, indicating whether to force the drawing of triangle ends at each 
of the extremes of the colour bar. This choice is automatically made from 
the provided 'brks', 'bar_limits', 'var_limits', 'col_inf' and 'col_sup', 
but the behaviour can be manually forced to draw or not to draw the triangle 
ends with this parameter. If 'cols' is provided, 'col_inf' and 'col_sup' 
will take priority over 'triangle_ends' when deciding whether to draw the 
triangle ends or not.}

\item{col_inf}{Parameter from s2dv::ColorBar. Colour to fill the inferior 
triangle end with. Useful if specifying colours manually with parameter 
'cols', to specify the colour and to trigger the drawing of the lower 
extreme triangle, or if 'cols' is not specified, to replace the colour 
automatically generated by ColorBar().}

\item{col_sup}{Parameter from s2dv::ColorBar. Colour to fill the superior 
triangle end with. Useful if specifying colours manually with parameter 
'cols', to specify the colour and to trigger the drawing of the upper 
extreme triangle, or if 'cols' is not specified, to replace the colour 
automatically generated by ColorBar().}

\item{col_unknown_map}{Colour to use to paint the grid cells for which a map 
is not possible to be chosen according to 'map_select_fun' or for those 
values that go beyond 'display_range'. Takes the value 'white' by default.}

\item{mask}{Optional numeric array with dimensions (latitude, longitude), with 
values in the range [0, 1], indicating the opacity of the mask over each 
grid point. Cells with a 0 will result in no mask, whereas cells with a 1 
will result in a totally opaque superimposed pixel coloured in 'col_mask'.}

\item{col_mask}{Colour to be used for the superimposed mask (if specified in 
'mask'). Takes the value 'grey' by default.}

\item{dots}{Array of same dimensions as 'var' or with dimensions 
c(n, dim(var)), where n is the number of dot/symbol layers to add to the 
plot. A value of TRUE at a grid cell will draw a dot/symbol on the 
corresponding square of the plot. By default all layers provided in 'dots' 
are plotted with dots, but a symbol can be specified for each of the 
layers via the parameter 'dot_symbol'.}

\item{bar_titles}{Optional vector of character strings providing the titles to 
be shown on top of each of the colour bars.}

\item{legend_scale}{Scale factor for the size of the colour bar labels. Takes 
1 by default.}

\item{cex_bar_titles}{Scale factor for the sizes of the bar titles. Takes 1.5 
by default.}

\item{plot_margin}{Numeric vector of length 4 for the margin sizes in the 
following order: bottom, left, top, and right. If not specified, use the 
default of par("mar"), c(5.1, 4.1, 4.1, 2.1). Used as 'margin_scale' in 
s2dv::PlotEquiMap.}

\item{bar_extra_margin}{Parameter from s2dv::ColorBar. Extra margins to be 
added around the colour bar, in the format c(y1, x1, y2, x2). The units are 
margin lines. Takes rep(0, 4) by default.}

\item{fileout}{File where to save the plot. If not specified (default) a 
graphics device will pop up. Extensions allowed: eps/ps, jpeg, png, pdf, bmp 
and tiff}

\item{width}{File width, in the units specified in the parameter size_units 
(inches by default). Takes 8 by default.}

\item{height}{File height, in the units specified in the parameter size_units 
(inches by default). Takes 5 by default.}

\item{size_units}{Units of the size of the device (file or window) to plot in. 
Inches ('in') by default. See ?Devices and the creator function of the 
corresponding device.}

\item{res}{Resolution of the device (file or window) to plot in. See ?Devices 
and the creator function of the corresponding device.}

\item{drawleg}{Where to draw the common colour bar. Can take values TRUE, 
FALSE or:\cr
'up', 'u', 'U', 'top', 't', 'T', 'north', 'n', 'N'\cr
'down', 'd', 'D', 'bottom', 'b', 'B', 'south', 's', 'S' (default)\cr
'right', 'r', 'R', 'east', 'e', 'E'\cr
'left', 'l', 'L', 'west', 'w', 'W'}

\item{return_leg}{A logical value indicating if the color bars information 
should be returned by the function. If TRUE, the function doesn't plot the
color bars but still creates the layout with color bar areas, and the 
arguments for GradientCatsColorBar() or ColorBar() will be returned. It is
convenient for users to adjust the color bars manually. The default is 
FALSE, the color bars will be plotted directly.}

\item{...}{Additional parameters to be passed on to \code{PlotEquiMap}.}
}
\description{
Plot a number a two dimensional matrices with (longitude, 
latitude) dimensions on a single map with the cylindrical equidistant 
latitude and longitude projection.
}
\examples{
# Simple example
x <- array(1:(20 * 10), dim = c(lat = 10, lon = 20)) / 200
a <- x * 0.6
b <- (1 - x) * 0.6
c <- 1 - (a + b)
lons <- seq(0, 359.5, length = 20)
lats <- seq(-89.5, 89.5, length = 10)
\dontrun{
PlotCombinedMap(list(a, b, c), lons, lats, 
               toptitle = 'Maximum map',
               map_select_fun = max,
               display_range = c(0, 1),
               bar_titles = paste('\% of belonging to', c('a', 'b', 'c')), 
               brks = 20, width = 12, height = 10)
}

Lon <- c(0:40, 350:359)
Lat <- 51:26
data <- rnorm(51 * 26 * 3)
dim(data) <- c(map = 3, lon = 51, lat = 26)
mask <-  sample(c(0,1), replace = TRUE, size = 51 * 26)
dim(mask) <- c(lat = 26, lon = 51)
\dontrun{
PlotCombinedMap(data, lon = Lon, lat = Lat, map_select_fun = max,
               display_range = range(data), mask = mask,
               width = 14, height = 10) 
}

}
\seealso{
\code{PlotCombinedMap} and \code{PlotEquiMap}
}
\author{
Nicolau Manubens, \email{nicolau.manubens@bsc.es}

Veronica Torralba, \email{veronica.torralba@bsc.es}
}
