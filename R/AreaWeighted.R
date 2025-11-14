#'Calculate the spatial area-weighted average of multidimensional arrays.
#'
#'This function computes a spatial area-weighted average of n-dimensional arrays 
#'given a spatial mask containing the area (e.g. in m^2) for each grid box.
#'
#'@param data An object of class \code{s2dv_cube} with at least lon_dim 
#'  and lat_dim containing climate data.
#'@param area A multidimensional array with named dimensions (at least lon_dim 
#'  and lat_dim) containing area of each grid box. The resolution, length and
#'   order of the lon_dim and lat_dim dimensions should be identical as for
#'   \code{data}. Missing values are allowed and treated as 0s.
#'@param lon_dim A character string indicating the name of the longitudinal
#'  dimension. The default value is 'lon'.
#'@param lat_dim A character string indicating the name of the latitudinal
#'  dimension. The default value is 'lat'.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'@param extra_info A named list with additional metadata to add to the 
#'  s2dv_cube. There can be one entry for each dimension in 'area' that is
#'  not also present in 'data'.
#'  
#'@return An object of class \code{s2dv_cube} with same dimensions as in object 
#'\code{data}, except with \code{lon_dim} and \code{lat_dim} removed. Any extra 
#'dimensions in \code{area} not present in \code{data} are preserved in the output (e.g., 
#'region dimension).
#' 
#'@examples
#'data <- array(data = 1:10, dim = c(year = 10, lon = 5, lat = 3))
#'coords <- list(lat = 1:3, lon = 1:5)
#'data <- list(data = data, coords = coords)
#'attr(data, 'class') <- 's2dv_cube'
#'area <- array(data = 1:10, dim = c(lon = 5, lat = 3))
#'CST_AreaWeighted(data, area)
#'# With extra region dimension
#'area <- array(data = 1:10, dim = c(lon = 5, lat = 3, region = 4))
#'CST_AreaWeighted(data, area)
#'
#'@import multiApply
#'@export
CST_AreaWeighted <- function(data, area, lon_dim = 'lon', 
                             lat_dim = 'lat', ncores = NULL,
                             extra_info = NULL) {
  
  # Check 's2dv_cube'
  if (!inherits(data, "s2dv_cube")) {
    stop("Parameter 'data' must be of the class 's2dv_cube'.")
  }
  # Check 'extra_info'
  if (!is.null(extra_info)) {
    if (!is.list(extra_info)) {
      stop("Parameter 'extra_info' should be a list with named elements.")
    } else {
      if (!all(names(extra_info) %in% names(dim(area)))) {
        warning("'extra_info' contains some names not found in the dimensions",
                "of 'area'. These will not be used.")
      }
    }
  }
  
  res <- AreaWeighted(data$data, area, lon_dim = lon_dim,
                      lat_dim = lat_dim, ncores = ncores)
  
  data$data <- res
  # dims
  data$dims <- dim(res)
  new_dims <- names(data$dims)[!names(data$dims) %in% names(data$coords)]
  # coords
  data$coords[c(lat_dim, lon_dim)] <- NULL
  for (dimname in new_dims) {
    if (!is.null(extra_info)) {
      data$coords[dimname] <- extra_info[dimname]
      attr(data$coords[[dimname]], "values") <- TRUE
      attr(data$coords[[dimname]], "indices") <- FALSE
    } else {
      data$coords[dimname] <- 1:data$dims[[dimname]]
      attr(data$coords[[dimname]], "values") <- FALSE
    }
  }

  # metadata
  data$attrs$Variable$metadata[c(lat_dim, lon_dim)] <- NULL
  
  return(data)
  
} 
#'Calculate the spatial area-weighted average of multidimensional arrays.
#'
#'This function computes a spatial area-weighted average of n-dimensional arrays 
#'given a spatial mask containing the area (e.g. in m^2) for each grid box.
#'
#'@param data A multidimensional array with named dimensions (at least lon_dim 
#'  and lat_dim) containing climate data.
#'@param area A multidimensional array with named dimensions (at least lon_dim 
#'  and lat_dim) containing area of each grid box. The resolution, length and
#'   order of the lon_dim and lat_dim dimensions should be identical as for
#'   \code{data}. Missing values are allowed and treated as 0s.
#'@param lon_dim A character string indicating the name of the longitudinal
#'  dimension. The default value is 'lon'.
#'@param lat_dim A character string indicating the name of the latitudinal
#'  dimension. The default value is 'lat'.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'  
#'@return An array with the same dimensions as \code{data}, except with 
#'\code{lon_dim} and \code{lat_dim} removed. Any extra dimensions in 
#'\code{area} not present in \code{data} are preserved in the output (e.g., 
#'region dimension).
#' 
#' @examples
#'data <- array(data = 1:10, dim = c(year = 10, lon = 5, lat = 3))
#'area <- array(data = 1:10, dim = c(lon = 5, lat = 3))
#'AreaWeighted(data, area)
#'# With extra region dimension
#'area <- array(data = 1:10, dim = c(lon = 5, lat = 3, region = 4))
#'AreaWeighted(data, area)
#'
#'@import multiApply
#'@export
AreaWeighted <- function(data, area, lon_dim = 'lon', 
                         lat_dim = 'lat', ncores = NULL) {
  
  # Check inputs
  # data
  if (is.null(data)) {
    stop("Parameter 'data' cannot be NULL.")
  }
  if (!is.array(data)) {
    stop("Parameter 'data' must be a numeric array.")
  }
  dim_names <- names(dim(data))
  if (is.null(dim_names)) {
    stop("Parameter 'data' must have dimension names.")
  }
  # area
  if (is.null(area)) {
    stop("Parameter 'area' cannot be NULL.")
  }
  if (!is.array(area)) {
    stop("Parameter 'area' must be a numeric array.")
  }
  dim_names <- names(dim(area))
  if (is.null(dim_names)) {
    stop("Parameter 'area' must have dimension names.")
  }
  # lon_dim
  if (!is.character(lon_dim) | length(lon_dim) > 1) {
    stop("Parameter 'lon_dim' must be a character string.")
  }
  if (!lon_dim %in% names(dim(data))) {
    stop("Parameter 'lon_dim' is not found in 'data'.")
  }
  if (!lon_dim %in% names(dim(area))) {
    stop("Parameter 'lon_dim' is not found in 'area'.")
  }
  # lat_dim
  if (!is.character(lat_dim) | length(lat_dim) > 1) {
    stop("Parameter 'lat_dim' must be a character string.")
  }
  if (!lat_dim %in% names(dim(data))) {
    stop("Parameter 'lat_dim' is not found in 'data'.")
  }
  if (!lat_dim %in% names(dim(area))) {
    stop("Parameter 'lat_dim' is not found in 'area'.")
  }
  # ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) | length(ncores) > 1) {
      stop("Parameter 'ncores' must be either NULL or a positive integer.")
    } else if (ncores %% 1 != 0 | ncores <= 0) {
      stop("Parameter 'ncores' must be a positive integer.")
    }
  }
  
  # Computation
  res <- multiApply::Apply(data = list(data = data, area = area),
                           target_dims = c(lon_dim, lat_dim),
                           fun = .AreaWeighted,
                           ncores = ncores)$output1
  return(res)
}

.AreaWeighted <- function(data, area){
  
  x <- data * area
  if (all(is.na(x))) {
    return(NA_real_)
  } else {
    denom <- sum(area, na.rm = TRUE)
    if (denom == 0) return(NA_real_)
    return(sum(x, na.rm = TRUE) / denom)
  }
  
}
