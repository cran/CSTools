#'CSTools data retrieval function using Start
#' 
#'This function aggregates, subsets and retrieves sub-seasonal, seasonal, 
#'decadal or climate projection data from NetCDF files in a local file system 
#'and arranges it for easy application of the CSTools functions. It calls the 
#'function \code{Start} from startR, which is an R package started at BSC with
#'the aim to develop a tool that allows the user to automatically process large 
#'multidimensional distributed data sets. Then, the output is transformed into 
#'`s2dv_cube` object.
#'
#'It receives any number of parameters (`...`) that are automatically forwarded 
#'to the `startR::Start` function. See details in `?startR::Start`. The 
#'auxiliary startR functions (e.g. indices(), values(), Sort(), CircularSort(), 
#'CDORemapper(), ...) can be used to define dimensions.
#' 
#'CST_Start() uses as.s2dv_cube() to transform the output into an s2dv_cube
#'object. The as.s2dv_cube() function is designed to be used with data that 
#'has been retrieved into memory. To avoid errors, please ensure that 
#'CST_Start(..., retrieve = TRUE) is specified. 
#'
#'@param ... Parameters that are automatically forwarded to the `startR::Start` 
#'  function. See details in `?startR::Start`.
#'@return An object of class \code{'s2dv_cube'} containing the subsetted and aggregated 
#'  data retrieved from the NetCDF files using \code{startR}. The data structure 
#'  is compatible with CSTools functions.
#'@examples 
#'\dontrun{
#'  sdates <- c('20101101', '20111101', '20121101')
#'  latmin <- 44
#'  latmax <- 47
#'  lonmin <- 6
#'  lonmax <- 9
#'  data <- CST_Start(dat = path,
#'                    var = 'prlr',
#'                    ensemble = indices(1:6),
#'                    sdate = sdates,
#'                    time = 121:151,
#'                    latitude = values(list(latmin, latmax)),
#'                    longitude = values(list(lonmin, lonmax)),
#'                    synonims = list(longitude = c('lon', 'longitude'),
#'                                    latitude = c('lat', 'latitude')),
#'                    return_vars = list(time = 'sdate',
#'                                       longitude = NULL, latitude = NULL),
#'                    retrieve = TRUE)
#'} 
#'\dontshow{
#' exp <- CSTools::lonlat_temp_st$exp
#' obs <- CSTools::lonlat_temp_st$obs
#' data <- CSTools::lonlat_prec
#'}
#'@import startR
#'@export
CST_Start <- function(...) {
  inputs <- as.list(substitute(list(...)))[-1]
  for (i in 1:length(inputs)) {
    inputs[[i]] <- eval(inputs[[i]])
  }
  res <- do.call(Start, inputs) 
  res <- as.s2dv_cube(res)
  return(res)
}

