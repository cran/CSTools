#' Bind two objects of class s2dv_cube
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
#'@param sort_dates A logical indicating whether or not to sort the data
#    and metadata in chronological order when binding objects along a time
#    dimension. The default value is FALSE.
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

CST_BindDim <- function(x, along, var_dim = NULL, dat_dim = NULL,
                        sort_dates = FALSE) {
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

  # Check time_bounds.
  time_bounds <- lapply(x, 
                        function(x) {
                          x$attrs$time_bounds
                        })
  if (any(time_bounds %in% list(NULL)) && !all(time_bounds %in% list(NULL))) {
    stop("Some of the s2dv_cubes have '$attrs$time_bounds', but others don't.",
         " The structure of all the s2dv_cubes must be consistent.")
  }

  # Check sort_dates
  if ((!is.logical(sort_dates))) {
    stop("Parameter 'sort_dates' must be a logical value (TRUE or FALSE).")
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
  res$attrs <- x[[1]]$attrs
  attrs <- lapply(x, "[[", "attrs")
  
  if (along %in% names(dim(x[[1]]$attrs$Dates))) {
    original_timezone <- attr(x[[1]]$attrs$Dates[1], "tzone")
    res$attrs$Dates <- BindDim(x = lapply(attrs, "[[", "Dates"),
                               along = along)
    # Transform dates back to POSIXct
    res$attrs$Dates <- as.POSIXct(res$attrs$Dates,
                                  origin = "1970-01-01",
                                  tz = original_timezone)
    # Time bounds
    if (!is.null(x[[1]]$attrs$time_bounds) && along %in% names(dim(x[[1]]$attrs$time_bounds[[1]]))) {
      for (bound in c("start", "end")) {
        res$attrs$time_bounds[[bound]] <- BindDim(x = lapply(
                                                        attrs,
                                                        function(x, bound) {
                                                          return(x$time_bounds[[bound]])
                                                        },
                                                        bound = bound),
                                                  along = along)
        res$attrs$time_bounds[[bound]] <- as.POSIXct(res$attrs$time_bounds[[bound]],
                                                     origin = "1970-01-01",
                                                     tz = original_timezone)
      }
}

    if (sort_dates) {

      dates_dims <- dim(res$attrs$Dates)
      dates_dimnames <- names(dates_dims)
      time_pos <- which(dates_dimnames == along)

      # Get all dates as a vector with their original indices
      dates_vector <- as.vector(res$attrs$Dates)
      n_total <- length(dates_vector)

      # Create indices for the original positions in multi-dimensional array
      indices_df <- expand.grid(lapply(dates_dims, seq_len))
      names(indices_df) <- dates_dimnames
      indices_df$date_value <- dates_vector
      indices_df$original_idx <- seq_len(n_total)

      # Sort by date value
      indices_df <- indices_df[order(indices_df$date_value), ]

      # Get the sorted dates
      sorted_dates <- indices_df$date_value

      # Reshape sorted dates back to original dimensions, preserving POSIXct class
      res$attrs$Dates <- array(sorted_dates, dim = dates_dims, dimnames = dimnames(res$attrs$Dates))
      class(res$attrs$Dates) <- c("POSIXct", "POSIXt")
      attr(res$attrs$Dates, "tzone") <- original_timezone

      data_vector <- as.vector(res$data)
      data_dims <- dim(res$data)
      data_dimnames <- names(data_dims)

      # Now, remap the dimensions of data according to the new order of dates
      common_dims <- intersect(dates_dimnames, data_dimnames)

      data_indices_df <- expand.grid(lapply(data_dims, seq_len))
      names(data_indices_df) <- data_dimnames
      data_indices_df$original_idx <- seq_len(nrow(data_indices_df))

      date_reorder_map <- indices_df[, c(dates_dimnames, "original_idx")]
      date_reorder_map$new_idx <- seq_len(nrow(date_reorder_map))

      # Merge with data indices to create full reordering
      merge_cols <- common_dims
      data_with_map <- merge(data_indices_df,
        date_reorder_map[, c(merge_cols, "original_idx", "new_idx")],
        by = merge_cols,
        suffixes = c("_data", "_date")
      )

      data_with_map <- data_with_map[order(data_with_map$new_idx, data_with_map$original_idx_data), ]

      # Extract the reordering for data
      data_reorder <- data_with_map$original_idx_data

      # Apply reordering to data
      sorted_data <- data_vector[data_reorder]
      res$data <- array(sorted_data, dim = data_dims, dimnames = dimnames(res$data))
    }
  } else if (!is.null(var_dim) && along == var_dim) {
    var_names <- sapply(x,
      function(x) {
        return(x$attrs$Variable$varName)
                        })
    res$attrs$Variable$varName <- var_names
    for (i in 1:length(x)) {
      var_name <- var_names[i]
      res$attrs$Variable$metadata[var_name] <- x[[i]]$attrs$Variable$metadata[var_name]
    }
    unique_vars <- unique(c(var_names, names(res$attrs$Variable$metadata)))
    res$attrs$Variable$metadata <- res$attrs$Variable$metadata[unique_vars]
  } else if (!is.null(dat_dim) && along == dat_dim) {
    res$attrs$Datasets <- sapply(x,
      function(x) {
        return(x$attrs$Datasets)
                                 })
  }
  # Bind medatada
  for (variable in 1:length(names(x[[1]]$attrs$Variable$metadata))) {
    if (any(along %in% names(dim(x[[1]]$attrs$Variable$metadata[[variable]])))) {
      metadata <- lapply(x,
                         function(x) {
                           return(x$attrs$Variable$metadata[[variable]])
                         })
      bound_metadata <- as.POSIXct(BindDim(metadata, along = along),
                                   origin = "1970-01-01",
                                   original_timezone)
      res$attrs$Variable$metadata[[along]] <- bound_metadata
    }
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
