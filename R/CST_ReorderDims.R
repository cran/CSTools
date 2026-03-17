#'Reorder the dimensions of an s2dv_cube
#'
#'Reorder the dimensions of an s2dv_cube. The order can be provided
#'either as indices or the dimension names. If the order is dimension names, 
#'the function looks for names(dim(x)). If it doesn't exist, the function checks
#' if attributes "dimensions" exists; this attribute is in the objects generated
#' by Load(). 
#'
#'@param x An array of which the dimensions to be reordered.
#'@param order A vector of indices or character strings indicating the new 
#'  order of the dimensions.
#'
#'@return An s2dv_cube which has the data and information as parameter 'x',
#'  but with different dimension order.
#'
#'@examples
#'# Example with sample data:
#'# Check original dimensions and coordinates
#'lonlat_temp$exp$dims
#'names(lonlat_temp$exp$coords)
#'# Reorder the s2dv_cube
#'new_order <- c("dataset", "member", "sdate", "ftime", "lon", "lat")
#'res <- CST_ReorderDims(x = lonlat_temp$exp, order = new_order)
#'# Check new dimensions and coordinates
#'dim(res$data)
#'res$dims
#'names(res$coords)
#' 
#'@seealso \link[s2dv]{Reorder}
#'
#'@importFrom s2dv Reorder
#'@export

CST_ReorderDims <- function(x, order) {

  # Check that x is s2dv_cube
  if (!inherits(x, 's2dv_cube')) {
    stop("Parameter 'x' must be of the class 's2dv_cube'.")
  }
  
  # Check that order is  A vector of indices or character strings
  if (!(is.numeric(order) || is.character(order))) {
    stop("'order' should be a vector of indices or character strings ",
    "indicating the new order of the dimensions")
  }

  # Check that dimension order is consistent
  if (!identical(names(dim(x$data)), names(x$dims))) {
    stop("The dimension names in 'x$data' and 'x$dims' are different or",
         " in a different order. Please check that the data and dims are",
         " consistent.")
  }

  # Reorder data
  x$data <- s2dv::Reorder(data = x$data, order = order)

  # Reorder dims
  x$dims <- dim(x$data)

  # Reorder coordinates
  x$coords <- x$coords[order]

  # Attributes
  time_dims <- names(dim(x$attrs$Dates))
  if (is.character(order)) {
    new_time_dims <- order[which(order %in% time_dims)] 
  } else {
    new_time_dims <- names(x$dims)[names(x$dims) %in% time_dims]
  }
  x$attrs$Dates <- s2dv::Reorder(data = x$attrs$Dates, order = new_time_dims)
  if (!is.null(x$attrs$time_bounds)) {
    x$attrs$time_bounds[[1]] <- s2dv::Reorder(data = x$attrs$time_bounds[[1]],
                                              order = new_time_dims)
    x$attrs$time_bounds[[2]] <- s2dv::Reorder(data = x$attrs$time_bounds[[2]],
                                              order = new_time_dims)
  }
  if (!is.null(x$attrs$Variable$metadata[['time']])) {
    x$attrs$Variable$metadata[['time']] <- s2dv::Reorder(data = x$attrs$Variable$metadata[['time']],
                                                         order = new_time_dims)
  }
  return(x)
}
