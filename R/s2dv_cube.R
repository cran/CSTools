#'Creation of a 's2dv_cube' object
#'
#'@description This function allows to create a 's2dv_cube' object by passing information through its parameters. This function will be needed if the data hasn't been loaded using CST_Load or has been transformed with other methods. A 's2dv_cube' object has many different components including metadata. This function will allow to create 's2dv_cube' objects even if not all elements are defined and for each expected missed parameter a warning message will be returned.
#'
#'@author Perez-Zanon Nuria, \email{nuria.perez@bsc.es}
#'
#'@param data an array with any number of named dimensions, typically an object output from CST_Load, with the following dimensions: dataset, member, sdate, ftime, lat and lon.
#'@param lon an array with one dimension containing the longitudes and attributes: dim, cdo_grid_name, data_across_gw, array_across_gw, first_lon, last_lon and projection.
#'@param lat an array with one dimension containing the latitudes and attributes: dim, cdo_grid_name, first_lat, last_lat and projection.
#'@param Variable a list of two elements: \code{varName} a character string indicating the abbreviation of a variable name and \code{level} a character string indicating the level (e.g., "2m"), if it is not required it could be set as NULL.
#'@param Datasets a named list with the dataset model with two elements: \code{InitiatlizationDates}, containing a list of the start dates for each member named with the names of each member, and \code{Members} containing a vector with the member names (e.g., "Member_1")
#'@param Dates a named list of two elements: \code{start}, an array of dimensions (sdate, time) with the POSIX initial date of each forecast time of each starting date, and \code{end}, an array of dimensions (sdate, time) with the POSIX final date of each forecast time of each starting date.
#'@param when a time stamp of the date issued by the Load() call to obtain the data.
#'@param source_files a vector of character strings with complete paths to all the found files involved in the Load() call.
#'
#'@return The function returns an object of class 's2dv_cube'.
#'
#'@seealso \code{\link[s2dverification]{Load}} and \code{\link{CST_Load}}
#'@examples
#'exp_original <- 1:100
#'dim(exp_original) <- c(lat = 2, time = 10, lon = 5)
#'exp1 <- s2dv_cube(data = exp_original)
#'class(exp1)
#'exp2 <- s2dv_cube(data = exp_original, lon = seq(-10, 10, 5), lat = c(45, 50)) 
#'class(exp2)
#'exp3 <- s2dv_cube(data = exp_original, lon = seq(-10, 10, 5), lat = c(45, 50),
#'                   Variable = list(varName = 'tas', level = '2m')) 
#'class(exp3)
#'exp4 <- s2dv_cube(data = exp_original, lon = seq(-10, 10, 5), lat = c(45, 50),
#'                  Variable = list(varName = 'tas', level = '2m'),
#'                  Dates = list(start = paste0(rep("01", 10), rep("01", 10), 1990:1999),
#'                               end = paste0(rep("31", 10), rep("01", 10), 1990:1999)))  
#'class(exp4)
#'exp5 <- s2dv_cube(data = exp_original, lon = seq(-10, 10, 5), lat = c(45, 50),
#'                  Variable = list(varName = 'tas', level = '2m'),
#'                  Dates = list(start = paste0(rep("01", 10), rep("01", 10), 1990:1999),
#'                               end = paste0(rep("31", 10), rep("01", 10), 1990:1999)),
#'                  when = "2019-10-23 19:15:29 CET")  
#'class(exp5)
#'exp6 <- s2dv_cube(data = exp_original, lon = seq(-10, 10, 5), lat = c(45, 50),
#'                  Variable = list(varName = 'tas', level = '2m'),
#'                  Dates = list(start = paste0(rep("01", 10), rep("01", 10), 1990:1999),
#'                               end = paste0(rep("31", 10), rep("01", 10), 1990:1999)),
#'                  when = "2019-10-23 19:15:29 CET", 
#'                  source_files = c("/path/to/file1.nc", "/path/to/file2.nc"))  
#'class(exp6)
#'exp7 <- s2dv_cube(data = exp_original, lon = seq(-10, 10, 5), lat = c(45, 50),
#'                  Variable = list(varName = 'tas', level = '2m'),
#'                  Dates = list(start = paste0(rep("01", 10), rep("01", 10), 1990:1999),
#'                               end = paste0(rep("31", 10), rep("01", 10), 1990:1999)),
#'                  when = "2019-10-23 19:15:29 CET", 
#'                  source_files = c("/path/to/file1.nc", "/path/to/file2.nc"),
#'                  Datasets = list(
#'                    exp1 = list(InitializationsDates = list(Member_1 = "01011990", 
#'                                                            Members = "Member_1"))))  
#'class(exp7)
#'dim(exp_original) <- c(dataset = 1, member = 1, sdate = 2, ftime = 5, lat = 2, lon = 5)
#'exp8 <- s2dv_cube(data = exp_original, lon = seq(-10, 10, 5), lat = c(45, 50),
#'                  Variable = list(varName = 'tas', level = '2m'),
#'                  Dates = list(start = paste0(rep("01", 10), rep("01", 10), 1990:1999),
#'                               end = paste0(rep("31", 10), rep("01", 10), 1990:1999)))  
#'class(exp8)
#'@export
s2dv_cube <- function(data, lon = NULL, lat = NULL, Variable = NULL, Datasets = NULL,
             Dates = NULL, when = NULL, source_files = NULL) {
  
  if (is.null(data) | !is.array(data) | is.null(names(dim(data)))) {
    stop("Parameter 'data' must be an array with named dimensions.")
  }
  dims <- dim(data)
  if (is.null(lon)) {
     if (!any(c('lon', 'longitude') %in% names(dims))) {
        warning("Parameter 'lon' is not provided but data contains a ",
                "longitudinal dimension.")
     } else {
        warning("Parameter 'lon' is not provided so the data is from an ",
                "unknown location.")
     }
  }    
  if (is.null(lat)) {
     if (!any(c('lat', 'latitude') %in% names(dims))) {
        warning("Parameter 'lat' is not provided but data contains a ",
                "latitudinal dimension.")
     } else {
        warning("Parameter 'lat' is not provided so the data is from an ",
                "unknown location.")
     }
  }
  if (is.null(Variable)) {
    warning("Parameter 'Variable' is not provided so the metadata ",
            "of 's2dv_cube' object will be incomplete.")
  }
  if (is.null(Datasets)) {
    warning("Parameter 'Datasets' is not provided so the metadata ",
            "of 's2dv_cube' object will be incomplete.")
  }
  if (is.null(Dates)) {
     if (!any(c('time', 'ftime', 'sdate') %in% names(dims))) {
        warning("Parameter 'Dates' is not provided but data contains a ",
                "temporal dimension.")
     } else {
        warning("Parameter 'Dates' is not provided so the data is from an ",
                "unknown time period.")
     }
  }
  if (is.null(when)) {
      warning("Parameter 'when' is not provided so the metadata ",
            "of 's2dv_cube' object will be incomplete.")
  }
  if (is.null(source_files)) {
      warning("Parameter 'source_files' is not provided so the metadata ",
            "of 's2dv_cube' object will be incomplete.")
  }
  if (!is.null(Variable)) {
    if (!is.list(Variable)) {
       Variable <- list(Variable)
    }
    if (names(Variable)[1] != 'varName' | names(Variable)[2] != 'level') {
       warning("The name of the first elment of parameter 'Variable' is ",
               "expected to be 'varName' and the second 'level'.")
    }
    if (!is.character(Variable[[1]])) {
       warning("The element 'Varname' of parameter 'Variable' must be ",
               "a character.")
    }
  }
  # Dimensions comparison
  if (!is.null(lon)) {
     if (any(names(dims) %in% c('lon', 'longitude'))) {
        if (dims[(names(dims) %in% c('lon', 'longitude'))] != length(lon) &
           dims[(names(dims) %in% c('lon', 'longitude'))] != 1) {
           stop("Length of parameter 'lon' doesn't match the length of ",
                "longitudinal dimension in parameter 'data'.")
         }
      }
  }
  if (!is.null(lat)) { 
     if (any(names(dims) %in% c('lat', 'latitude'))) {
        if (dims[(names(dims) %in% c('lat', 'latitude'))] != length(lat) &
           dims[(names(dims) %in% c('lat', 'latitude'))] != 1) {
           stop("Length of parameter 'lat' doesn't match the length of ",
                "latitudinal dimension in parameter 'data'.")
         }
      }
  }
  if (!is.null(Dates)) {
     if (!is.list(Dates)) {
        stop("Parameter 'Dates' must be a list.")
     } else {
        if (length(Dates) > 2) {
          warning("Parameter 'Dates' is a list with more than 2 ",
                  "elements and only the first two will be used.")
          Dates <- Dates[1 : 2]
        } 
        if (names(Dates)[1] != 'start' | names(Dates)[2] != 'end') {
          warning("The name of the first element of parameter 'Dates' ",
                  "is expected to be 'start' and the second 'end'.")
        }
        if (length(Dates[[1]]) != length(Dates[[2]]) & 
            length(Dates) == 2) {
          stop("The length of the elements in parameter 'Dates' must ",
               "be equal.")
        }
        time_dims <- dims[names(dims) %in% c('sdate', 'time', 'ftime')]
        if (prod(time_dims) != length(Dates[[1]])) {
          stop("The length of the temporal dimension doesn't match ", 
               " with the length of elements in parameter 'Dates'.")
        }
     }
  }
  object <- list(data = data, lon = lon, lat = lat, Variable = Variable, 
                 Datasets = Datasets, Dates = Dates, when = when, 
                 source_files = source_files)
  class(object) <- 's2dv_cube'
  return(object)
}
