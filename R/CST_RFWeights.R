#' Compute climatological weights for RainFARM stochastic precipitation downscaling
#'
#' @author Jost von Hardenberg - ISAC-CNR, \email{j.vonhardenberg@isac.cnr.it}
#'
#' @description Compute climatological ("orographic") weights from a fine-scale precipitation climatology file.
#' @references Terzago, S., Palazzi, E., & von Hardenberg, J. (2018).
#' Stochastic downscaling of precipitation in complex orography: 
#' A simple method to reproduce a realistic fine-scale climatology.
#' Natural Hazards and Earth System Sciences, 18(11),
#' 2825-2840. http://doi.org/10.5194/nhess-18-2825-2018 .
#' @param climfile Filename of a fine-scale precipitation climatology.
#' The file is expected to be in NetCDF format and should contain
#' at least one precipitation field. If several fields at different times are provided,
#' a climatology is derived by time averaging.
#' Suitable climatology files could be for example a fine-scale precipitation climatology
#' from a high-resolution regional climate model (see e.g. Terzago et al. 2018), a local
#' high-resolution gridded climatology from observations, or a reconstruction such as those which 
#' can be downloaded from the WORLDCLIM (http://www.worldclim.org) or CHELSA (http://chelsa-climate.org)
#' websites. The latter data will need to be converted to NetCDF format before being used (see for example the GDAL tools (https://www.gdal.org).
#' @param nf Refinement factor for downscaling (the output resolution is increased by this factor).
#' @param lon Vector of longitudes.
#' @param lat Vector of latitudes.
#' The number of longitudes and latitudes is expected to be even and the same. If not
#' the function will perform a subsetting to ensure this condition.
#' @param varname Name of the variable to be read from \code{climfile}.
#' @param fsmooth Logical to use smooth conservation (default) or large-scale box-average conservation.
#' @return A matrix containing the weights with dimensions (lon, lat).
#' @import ncdf4
#' @import rainfarmr
#' @examples
#' # Create weights to be used with the CST_RainFARM() or RainFARM() functions
#' # using an external fine-scale climatology file.
#'
#' \dontrun{
#' # Specify lon and lat of the input
#' lon <- seq(10,13.5,0.5)
#' lat <- seq(40,43.5,0.5)
#' nf <- 8
#' ww <- CST_RFWeights("./worldclim.nc", nf, lon, lat, fsmooth = TRUE)
#' }
#' @export

CST_RFWeights <- function(climfile, nf, lon, lat, varname = "", fsmooth=TRUE) {

  # Ensure input  grid is square and with even dimensions
  if ((length(lat) != length(lon)) | (length(lon) %% 2 == 1)) {
    warning("Input data are expected to be on a square grid",
            " with an even number of pixels per side.")
    nmin <- min(length(lon), length(lat))
    nmin <- floor(nmin / 2) * 2
    lon <- lon[1:nmin]
    lat <- lat[1:nmin]
    warning("The input data have been cut to the range.")
    warning(paste0("lon: [", lon[1], ", ", lon[length(lon)], "] ",
                 " lat: [", lat[1], ", ", lat[length(lat)], "]"))
  }
  
  ncin <- nc_open(climfile) 
  latin <- ncvar_get(ncin, grep("lat", attributes(ncin$dim)$names, value = TRUE))
  lonin <- ncvar_get(ncin, grep("lon", attributes(ncin$dim)$names, value = TRUE))
  if ( varname == "") {
    varname <- grep("bnds", attributes(ncin$var)$names,
                    invert = TRUE, value = TRUE)[1]
  }
  zclim <- ncvar_get(ncin, varname)

# Check if lon and lat need to be reversed
  if ( lat[1] > lat[2] ) {
    lat <- rev(lat) 
    frev = TRUE
  } else {
    frev = FALSE
  }
  if ( latin[1] > latin[2] ) {
    latin <- rev(latin)
    zclim = zclim[,seq(dim(zclim)[2],1)]
  }
  ww <- rfweights(zclim, lonin, latin, lon, lat, nf, fsmooth = fsmooth)
  if ( frev ) {
    ww = ww[,seq(dim(ww)[2],1)]
  }
  attributes(dim(ww))$names <- c("lon","lat")
  return(ww)
}
