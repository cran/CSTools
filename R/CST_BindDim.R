#'Bind two objects of class s2dv_cube
#'
#'This function combines the data inside two or more objects of class
#'\code{s2dv_cube} along the specified \code{along} dimension, and modifies the
#'dimensions, coordinates and attributes accordingly, producing a result that
#'contains the complete metadata for all variables, time steps and spatial
#'coordinates that are bound in the process. It ensures that the information
#'inside the s2dv_cube remains coherent with the data it contains.\cr\cr
#'If the dimension specified in \code{along} is among the time dimensions in
#'\code{attrs$Dates}, the dates arrays are also bound along this dimension. 
#'The \code{load_parameters} and \code{when} attributes of the first cube
#'are preserved. The \code{source_files} attribute is bound along the
#'\code{var_dim} and \code{dat_dim} dimensions.
#' 
#'@author Agudetse Roures Victoria, \email{victoria.agudetse@bsc.es}
#'
#'@param x Two or more objects of class \code{s2dv_cube} to be bound together.
#'@param along A character string indicating the name of the binding dimension.
#'@param dat_dim A character string indicating the name of dataset dimension.
#'  The default value is NULL. Specifying this dimension ensures the dataset
#'  metadata is correctly preserved.
#'@param var_dim A character string indicating the name of the variable
#'  dimension. The default value is NULL. Specifying this dimension ensures
#   the variable metadata is correctly preserved.
#'
#'@return An object of class \code{s2dv_cube} with the combined data, 
#'  dimensions, coordinates and attributes of the elements in \code{x}.
#'
#'@examples
#'  # Example with sample data:
#'  # Check original dimensions and coordinates
#'  lonlat_temp$exp$dims
#'  lonlat_temp$obs$dims
#'  # Bind both datasets along the member dimension
#'  res <- CST_BindDim(x = list(lonlat_temp$exp, lonlat_temp$obs),
#'                     along = "member",
#'                     var_dim = NULL,
#'                     dat_dim = "dat")
#'  # Check new dimensions and coordinates
#'  res$dims
#'  names(res$coords)
#' 
#'@seealso \link[abind]{abind}
#'
#'@importFrom abind abind
#'@export

CST_BindDim <- function(x, along, var_dim = NULL, dat_dim = NULL) {
  # Check that at least two objects have been passed
  if (length(x) < 2) {
    stop("'x' must contain at least two 's2dv_cube' objects.")
  }
  # Check list
  if (!is.list(x)) {
    x <- as.list(x)
  }
  # Check that x elements are of s2dv_cube class
  if (!all(vapply(x, function (x) {inherits(x, "s2dv_cube")}, logical(1)))) {
    stop("All elements of 'x' must be of the class 's2dv_cube'.")
  }
  # Check 'along' dimension
  inner_dimensions <- lapply(x,
                             function(x) {
                               x$dims[names(x$dims) %in% along]
                             })
  if (any(lengths(inner_dimensions) == 0)) {
    stop("The dimension specified in 'along' must be present in all the ",
         "s2dv_cubes in 'x'.")
  }
  # Check that all s2dv_cubes have the same dimensions
  outer_dimensions <- lapply(x,
                             function(x) {
                               x$dims[!names(x$dims) %in% along]
                             })
  if (length(unique(outer_dimensions)) > 1) {
    stop("All elements of 'x' must have the same dimension names and length, ",
         "except for the dimension specified in 'along'.")
  }
  # Check var_dim
  if (!is.null(var_dim)) {
    if ((!is.character(var_dim)) || (length(var_dim) > 1)) {
      stop("Parameter 'var_dim' must be a character string.")
    }
  } else {
    warning("Parameter 'var_dim' not specified, variable metadata ",
            "of the output might be inconsistent with the input.")
  }
  # Check dat_dim
  if (!is.null(dat_dim)) {
    if ((!is.character(dat_dim)) || (length(dat_dim) > 1)) {
      stop("Parameter 'dat_dim' must be a character string.")
    }
  } else {
    warning("Parameter 'dat_dim' not specified, dataset metadata ",
            "of the output might be inconsistent with the input.")
  }

  # Bind data
  res <- list()
  res$data <- BindDim(x = lapply(x, "[[", "data"),
                   along = along)
  # Adjust dimensions
  res$dims <- dim(res$data)
  # Adjust coordinates and coordinate attributes
  ## TODO: Can probably be improved
  coords_attrs <- x[[1]]$coords
  res$coords <- coords_attrs
  # If indices or unsure, simply create new sequence of indices
  if (is.null(attr(coords_attrs[[along]], "indices")) ||
      attr(coords_attrs[[along]], "indices")) {
    res$coords[[along]] <- seq(1:res$dims[[along]])
    attr(res$coords[[along]], "indices") <- TRUE
  } else {
    res$coords[[along]] <- sapply(x,
                                  function(x, along) {
                                    return(x$coords[[along]])
                                  },
                                  along = along)
    attributes(res$coords[[along]]) <- attributes(coords_attrs[[along]])
  }
  # Variable
  res$attrs <- x$attrs
  attrs <- lapply(x, "[[", "attrs")
  if (along %in% names(dim(x[[1]]$attrs$Dates))) {
    original_timezone <- attr(x[[1]]$attrs$Dates[1], "tzone")
    res$attrs$Dates <- BindDim(x = lapply(attrs, "[[", "Dates"),
                               along = along)
    # Transform dates back to POSIXct
    res$attrs$Dates <- as.POSIXct(res$attrs$Dates,
                                  origin = "1970-01-01",
                                  tz = original_timezone)

  } else if (!is.null(var_dim) && along == var_dim) {
    var_names <- sapply(x,
                        function(x) {
                          return(x$attrs$Variable$varName)
                        })
    res$attrs$Variable$varName <- var_names
    for (i in 1:length(x)) {
      var_name <- var_names[i]
      res$attrs$Variable$metadata[var_name] <- x[i]$attrs$Variable$metadata[var_name]
    }
  } else if (!is.null(dat_dim) && along == dat_dim) {
    res$attrs$Datasets <- sapply(x,
                                 function(x) {
                                   return(x$attrs$Datasets)
                                 })
  }
  # Source files
  source_files <- lapply(attrs, "[[", "source_files")
  res$attrs$source_files <- as.vector(sapply(list(source_files), c))
  # When
  res$attrs$when <- Sys.time()
  # Class
  class(res) <- "s2dv_cube"
  return(res)
}

#'Bind two arrays by a specified named dimension
#'
#'This function combines the data inside two or more arrays with named
#'dimensions along the specified \code{along} dimension. It is a wrapper of the
#'abind() function from the abind package.
#'
#'@author Agudetse Roures Victoria, \email{victoria.agudetse@bsc.es}
#'
#'@param x A list of two or more arrays with named dimensions to be bound 
#'  together.
#'@param along A character string indicating the name of the binding dimension.
#'
#'@return A single array combining the arrays in \code{x} along the specified
#'  dimension, with dimension names. The order of the dimensions will be the
#'  same as in the first array provided in the list.
#'
#'@examples
#'array1 <- array(1:100, dim = c(time = 1, lon = 10, lat = 10))
#'array2 <- array(101:200, dim = c(lon = 10, lat = 10, time = 1))
#'# Bind arrays
#'array3 <- BindDim(x = list(array1, array2),
#'                  along = "time")
#'# Check new dimensions
#'dim(array3)
#' 
#'@seealso \link[abind]{abind}
#'
#'@importFrom abind abind
#'@importFrom s2dv Reorder
#'@export

BindDim <- function(x, along) {
  original_dims <- names(dim(x[[1]]))
  # All array dimensions must be in the same order
  x <- lapply(x, s2dv::Reorder, order = original_dims)
  # Bind arrays and restore dimension names
  res <- abind(x, along = which(original_dims == along))
  names(dim(res)) <- original_dims
  return(res)
}
