#'Save objects of class 's2dv_cube' to data in NetCDF format
#'
#'@author Perez-Zanon Nuria, \email{nuria.perez@bsc.es}
#'
#'@description This function allows to divide and save a object of class 
#''s2dv_cube' into a NetCDF file, allowing to reload the saved data using 
#'\code{Start} function from StartR package. If the original 's2dv_cube' object 
#'has been created from \code{CST_Load()}, then it can be reloaded with 
#'\code{Load()}.
#'
#'@param data An object of class \code{s2dv_cube}.
#'@param destination A character string containing the directory name in which 
#'  to save the data. NetCDF file for each starting date are saved into the 
#'  folder tree: \cr
#'  destination/Dataset/variable/. By default the function 
#'  creates and saves the data into the working directory.
#'@param sdate_dim A character string indicating the name of the start date 
#'  dimension. By default, it is set to 'sdate'. It can be NULL if there is no
#'  start date dimension.
#'@param ftime_dim A character string indicating the name of the forecast time  
#'  dimension. By default, it is set to 'time'. It can be NULL if there is no 
#'  forecast time dimension.
#'@param dat_dim A character string indicating the name of dataset dimension. 
#'  By default, it is set to 'dataset'. It can be NULL if there is no dataset  
#'  dimension.
#'@param var_dim A character string indicating the name of variable dimension. 
#'  By default, it is set to 'var'. It can be NULL if there is no variable  
#'  dimension.
#'@param memb_dim A character string indicating the name of the member dimension.
#'  By default, it is set to 'member'. It can be NULL if there is no member 
#'  dimension.  
#'@param single_file A logical value indicating if all object is saved in a 
#'  single file (TRUE) or in multiple files (FALSE). When it is FALSE, 
#'  the array is separated for Datasets, variable and start date. It is FALSE  
#'  by default.
#'@param extra_string A character string to be include as part of the file name, 
#'  for instance, to identify member or realization. It would be added to the 
#'  file name between underscore characters.
#'
#'@return Multiple or single NetCDF files containing the data array.\cr
#'\item{\code{single_file = TRUE}}{
#'  All data is saved in a single file located in the specified destination  
#'  path with the following name: 
#'  <variable_name>_<extra_string>_<first_sdate>_<last_sdate>.nc. Multiple 
#'  variables are saved separately in the same file. The forecast time units 
#'  is extracted from the frequency of the time steps (hours, days, months). 
#'  The first value of forecast time is 1. If no frequency is found, the units 
#'  will be 'hours since' each start date and the time steps are assumed to be 
#'  equally spaced.
#'}
#'\item{\code{single_file = FALSE}}{
#'  The data array is subset and stored into multiple files. Each file 
#'  contains the data subset for each start date, variable and dataset. Files 
#'  with different variables and Datasets are stored in separated directories 
#'  within the following directory tree: destination/Dataset/variable/. 
#'  The name of each file will be: 
#'  <variable_name>_<extra_string>_<sdate>.nc.
#'}
#' 
#'@seealso \code{\link[startR]{Start}}, \code{\link{as.s2dv_cube}} and 
#'\code{\link{s2dv_cube}}
#'
#'@examples
#'\dontrun{
#'data <- lonlat_temp$exp
#'destination <- "./"
#'CST_SaveExp(data = data, destination = destination, ftime_dim = 'ftime', 
#'            var_dim = NULL, ftime_dim = 'ftime', var_dim = NULL)
#'}
#'
#'@import ncdf4
#'@importFrom s2dv Reorder
#'@importFrom ClimProjDiags Subset
#'@import multiApply
#'@export
CST_SaveExp <- function(data, destination = "./", sdate_dim = 'sdate',  
                        ftime_dim = 'time', dat_dim = 'dataset',
                        var_dim = 'var', memb_dim = 'member', 
                        single_file = FALSE, extra_string = NULL) {
  # Check 's2dv_cube'
  if (!inherits(data, 's2dv_cube')) {
    stop("Parameter 'data' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load.")
  }
  # Check object structure
  if (!all(c('data', 'attrs') %in% names(data))) {
    stop("Parameter 'data' must have at least 'data' and 'attrs' elements ",
         "within the 's2dv_cube' structure.")
  }
  if (!inherits(data$attrs, 'list')) {
    stop("Level 'attrs' must be a list with at least 'Dates' element.")
  }
  if (!all(c('coords') %in% names(data))) {
    warning("Element 'coords' not found. No coordinates will be used.")
  }
  # metadata
  if (is.null(data$attrs$Variable$metadata)) {
    warning("No metadata found in element Variable from attrs.")
  } else {
    if (!inherits(data$attrs$Variable$metadata, 'list')) {
      stop("Element metadata from Variable element in attrs must be a list.")
    }
    if (!any(names(data$attrs$Variable$metadata) %in% names(data$coords))) {
      warning("Metadata is not found for any coordinate.")
    } else if (!any(names(data$attrs$Variable$metadata) %in% 
                    data$attrs$Variable$varName)) {
      warning("Metadata is not found for any variable.")
    }
  }
  # Dates
  if (is.null(data$attrs$Dates)) {
    stop("Element 'Dates' from 'attrs' level cannot be NULL.")
  }
  if (is.null(dim(data$attrs$Dates))) {
    stop("Element 'Dates' from 'attrs' level must have time dimensions.")
  }
  # sdate_dim
  if (!is.null(sdate_dim)) {
    if (!is.character(sdate_dim)) {
      stop("Parameter 'sdate_dim' must be a character string.")
    }
    if (length(sdate_dim) > 1) {
      warning("Parameter 'sdate_dim' has length greater than 1 and ",
              "only the first element will be used.")
      sdate_dim <- sdate_dim[1]
    }
  } else {
    if (length(dim(data$attrs$Dates)) == 1) {
      sdate_dim <- 'sdate'
      dim(data$data) <- c(sdate = 1, dim(data$data))
      data$dims <- dim(data$data)
      dim(data$attrs$Dates) <- c(sdate = 1, dim(data$attrs$Dates))
      data$coords[[sdate_dim]] <- data$attrs$Dates[1]
    }
  }

  SaveExp(data = data$data,
          destination = destination, 
          Dates = data$attrs$Dates, 
          coords = data$coords,
          varname = data$attrs$Variable$varName,
          metadata = data$attrs$Variable$metadata,
          Datasets = data$attrs$Datasets, 
          startdates = data$coords[[sdate_dim]],
          dat_dim = dat_dim, sdate_dim = sdate_dim, 
          ftime_dim = ftime_dim, var_dim = var_dim, 
          memb_dim = memb_dim,
          extra_string = extra_string, 
          single_file = single_file)
}
#'Save a multidimensional array with metadata to data in NetCDF format
#'@description This function allows to save a data array with metadata into a 
#'NetCDF file, allowing to reload the saved data using \code{Start} function 
#'from StartR package. If the original 's2dv_cube' object has been created from 
#'\code{CST_Load()}, then it can be reloaded with \code{Load()}.
#'
#'@author Perez-Zanon Nuria, \email{nuria.perez@bsc.es}
#'
#'@param data A multi-dimensional array with named dimensions.
#'@param destination A character string indicating the path where to store the 
#'  NetCDF files.
#'@param Dates A named array of dates with the corresponding sdate and forecast 
#'  time dimension.
#'@param coords A named list with elements of the coordinates corresponding to 
#'  the dimensions of the data parameter. The names and length of each element 
#'  must correspond to the names of the dimensions. If any coordinate is not 
#'  provided, it is set as an index vector with the values from 1 to the length 
#'  of the corresponding dimension.
#'@param varname A character string indicating the name of the variable to be 
#'  saved.
#'@param metadata A named list where each element is a variable containing the
#'  corresponding information. The information must be contained in a list of 
#'  lists for each variable.
#'@param Datasets A vector of character string indicating the names of the 
#'  datasets.
#'@param startdates A vector of dates indicating the initialization date of each 
#'  simulations.
#'@param sdate_dim A character string indicating the name of the start date 
#'  dimension. By default, it is set to 'sdate'. It can be NULL if there is no
#'  start date dimension.
#'@param ftime_dim A character string indicating the name of the forecast time  
#'  dimension. By default, it is set to 'time'. It can be NULL if there is no 
#'  forecast time dimension.
#'@param dat_dim A character string indicating the name of dataset dimension. 
#'  By default, it is set to 'dataset'. It can be NULL if there is no dataset  
#'  dimension.
#'@param var_dim A character string indicating the name of variable dimension. 
#'  By default, it is set to 'var'. It can be NULL if there is no variable  
#'  dimension.
#'@param memb_dim A character string indicating the name of the member dimension.
#'  By default, it is set to 'member'. It can be NULL if there is no member 
#'  dimension.  
#'@param single_file A logical value indicating if all object is saved in a 
#'  unique file (TRUE) or in separated directories (FALSE). When it is FALSE, 
#'  the array is separated for Datasets, variable and start date. It is FALSE  
#'  by default.
#'@param extra_string A character string to be include as part of the file name, 
#'  for instance, to identify member or realization. It would be added to the 
#'  file name between underscore characters.
#'
#'@return Multiple or single NetCDF files containing the data array.\cr
#'\item{\code{single_file = TRUE}}{
#'  All data is saved in a single file located in the specified destination  
#'  path with the following name: 
#'  <variable_name>_<extra_string>_<first_sdate>_<last_sdate>.nc. Multiple 
#'  variables are saved separately in the same file. The forecast time units 
#'  is extracted from the frequency of the time steps (hours, days, months). 
#'  The first value of forecast time is 1. If no frequency is found, the units 
#'  will be 'hours since' each start date and the time steps are assumed to be 
#'  equally spaced.
#'}
#'\item{\code{single_file = FALSE}}{
#'  The data array is subset and stored into multiple files. Each file 
#'  contains the data subset for each start date, variable and dataset. Files 
#'  with different variables and Datasets are stored in separated directories 
#'  within the following directory tree: destination/Dataset/variable/. 
#'  The name of each file will be: 
#'  <variable_name>_<extra_string>_<sdate>.nc.
#'}
#' 
#'@examples
#'\dontrun{
#'data <- lonlat_temp$exp$data
#'lon <- lonlat_temp$exp$coords$lon
#'lat <- lonlat_temp$exp$coords$lat
#'coords <- list(lon = lon, lat = lat)
#'Datasets <- lonlat_temp$exp$attrs$Datasets
#'varname <- 'tas'
#'Dates <- lonlat_temp$exp$attrs$Dates
#'destination = './'
#'metadata <- lonlat_temp$exp$attrs$Variable$metadata
#'SaveExp(data = data, destination = destination, coords = coords, 
#'        Datasets = Datasets, varname = varname, Dates = Dates, 
#'        metadata = metadata, single_file = TRUE, ftime_dim = 'ftime', 
#'        var_dim = NULL)
#'}
#' 
#'@import ncdf4
#'@importFrom s2dv Reorder
#'@import multiApply
#'@importFrom ClimProjDiags Subset
#'@export
SaveExp <- function(data, destination = "./", Dates = NULL, coords = NULL, 
                    varname = NULL, metadata = NULL, Datasets = NULL, 
                    startdates = NULL, dat_dim = 'dataset', sdate_dim = 'sdate', 
                    ftime_dim = 'time', var_dim = 'var', memb_dim = 'member',
                    single_file = FALSE, extra_string = NULL) {
  ## Initial checks
  # data
  if (is.null(data)) {
    stop("Parameter 'data' cannot be NULL.")
  }
  dimnames <- names(dim(data))
  if (is.null(dimnames)) {
    stop("Parameter 'data' must be an array with named dimensions.")
  }
  # destination
  if (!is.character(destination) | length(destination) > 1) {
    stop("Parameter 'destination' must be a character string of one element ",
         "indicating the name of the file (including the folder if needed) ",
         "where the data will be saved.")
  }
  # Dates
  if (!is.null(Dates)) {
    if (!inherits(Dates, "POSIXct") & !inherits(Dates, "Date")) {
      stop("Parameter 'Dates' must be of 'POSIXct' or 'Dates' class.")
    }
    if (is.null(dim(Dates))) {
      stop("Parameter 'Dates' must have dimension names.")
    }
  }
  # coords
  if (!is.null(coords)) {
    if (!all(names(coords) %in% dimnames)) {
      coords <- coords[-which(!names(coords) %in% dimnames)]
    }
    for (i_coord in dimnames) {
      if (i_coord %in% names(coords)) {
        if (length(coords[[i_coord]]) != dim(data)[i_coord]) {
          warning(paste0("Coordinate '", i_coord, "' has different lenght as ",
                         "its dimension and it will not be used."))
          coords[[i_coord]] <- 1:dim(data)[i_coord]
        }
      } else {
        warning(paste0("Coordinate '", i_coord, "' is not provided ",
                       "and it will be set as index in element coords."))
        coords[[i_coord]] <- 1:dim(data)[i_coord]
      }
    }
  } else {
    coords <- sapply(dimnames, function(x) 1:dim(data)[x])
  }
  # varname
  if (is.null(varname)) {
    warning("Parameter 'varname' is NULL. It will be assigned to 'X'.")
    varname <- 'X'
  } else if (length(varname) > 1) {
    multiple_vars <- TRUE
  } else {
    multiple_vars <- FALSE
  }
  if (!all(sapply(varname, is.character))) {
    stop("Parameter 'varname' must be a character string with the ",
         "variable names.")
  }
  # metadata
  if (is.null(metadata)) {
    warning("Parameter 'metadata' is not provided so the metadata saved ",
            "will be incomplete.")
  }
  # single_file
  if (!inherits(single_file, 'logical')) {
    warning("Parameter 'single_file' must be a logical value. It will be ", 
            "set as FALSE.")
    single_file <- FALSE
  }
  # extra_string
  if (!is.null(extra_string)) {
    if (!is.character(extra_string)) {
      stop("Parameter 'extra_string' must be a character string.")
    }
  }

  ## Dimensions checks
  # Spatial coordinates
  if (!any(dimnames %in% .KnownLonNames()) | 
      !any(dimnames %in% .KnownLatNames())) {
    warning("Spatial coordinate names do not match any of the names accepted by ",
            "the package.")
    lon_dim <- NULL
    lat_dim <- NULL
  } else {
    lon_dim <- dimnames[which(dimnames %in% .KnownLonNames())]
    lat_dim <- dimnames[which(dimnames %in% .KnownLatNames())]
    if (length(lon_dim) > 1) {
      warning("Found more than one longitudinal dimension. Only the first one ", 
              "will be used.")
      lon_dim <- lon_dim[1]
    }
    if (length(lat_dim) > 1) {
      warning("Found more than one latitudinal dimension. Only the first one ", 
              "will be used.")
      lat_dim <- lat_dim[1]
    }
  }
  # ftime_dim
  if (!is.null(ftime_dim)) {
    if (!is.character(ftime_dim)) {
      stop("Parameter 'ftime_dim' must be a character string.")
    }
    if (!all(ftime_dim %in% dimnames)) {
      stop("Parameter 'ftime_dim' is not found in 'data' dimension.")
    }
    if (length(ftime_dim) > 1) {
      warning("Parameter 'ftime_dim' has length greater than 1 and ",
              "only the first element will be used.")
      ftime_dim <- ftime_dim[1]
    }
  }
  # sdate_dim
  if (!is.null(sdate_dim)) {
    if (!is.character(sdate_dim)) {
      stop("Parameter 'sdate_dim' must be a character string.")
    }
    if (length(sdate_dim) > 1) {
      warning("Parameter 'sdate_dim' has length greater than 1 and ",
              "only the first element will be used.")
      sdate_dim <- sdate_dim[1]
    }
    if (!all(sdate_dim %in% dimnames)) {
      stop("Parameter 'sdate_dim' is not found in 'data' dimension.")
    }
  }
  # memb_dim
  if (!is.null(memb_dim)) {
    if (!is.character(memb_dim)) {
      stop("Parameter 'memb_dim' must be a character string.")
    }
    if (!all(memb_dim %in% dimnames)) {
      stop("Parameter 'memb_dim' is not found in 'data' dimension. Set it ", 
           "as NULL if there is no member dimension.")
    }
  }
  # dat_dim
  if (!is.null(dat_dim)) {
    if (!is.character(dat_dim)) {
      stop("Parameter 'dat_dim' must be a character string.")
    }
    if (!all(dat_dim %in% dimnames)) {
      stop("Parameter 'dat_dim' is not found in 'data' dimension. Set it ", 
           "as NULL if there is no Datasets dimension.")
    }
    if (length(dat_dim) > 1) {
      warning("Parameter 'dat_dim' has length greater than 1 and ",
              "only the first element will be used.")
      dat_dim <- dat_dim[1]
    }
    n_datasets <- dim(data)[dat_dim]
  } else {
    n_datasets <- 1
  }
  # var_dim
  if (!is.null(var_dim)) {
    if (!is.character(var_dim)) {
      stop("Parameter 'var_dim' must be a character string.")
    }
    if (!all(var_dim %in% dimnames)) {
      stop("Parameter 'var_dim' is not found in 'data' dimension. Set it ", 
           "as NULL if there is no variable dimension.")
    }
    if (length(var_dim) > 1) {
      warning("Parameter 'var_dim' has length greater than 1 and ",
              "only the first element will be used.")
      var_dim <- var_dim[1]
    }
    n_vars <- dim(data)[var_dim]
  } else {
    n_vars <- 1
  }
  # minimum dimensions
  if (all(dimnames %in% c(var_dim, dat_dim))) {
    if (!single_file) {
      warning("Parameter data has only ", 
              paste(c(var_dim, dat_dim), collapse = ' and '), " dimensions ", 
              "and it cannot be splitted in multiple files. All data will ", 
              "be saved in a single file.")
      single_file <- TRUE
    }
  }
  # Dates dimension check
  if (!is.null(Dates)) {
    if (all(names(dim(Dates)) == c(ftime_dim, sdate_dim)) | 
        all(names(dim(Dates)) == c(sdate_dim, ftime_dim))) {
      if (is.null(startdates)) {
        startdates <- Subset(Dates, along = ftime_dim, 1, drop = 'selected')
      } else if ((!inherits(startdates, "POSIXct") & !inherits(startdates, "Date")) &&
                 (!is.character(startdates) | (all(nchar(startdates) != 10) &
                  all(nchar(startdates) != 8) & all(nchar(startdates) != 6)))) {
        warning("Parameter 'startdates' should be a character string containing ", 
                "the start dates in the format 'yyyy-mm-dd', 'yyyymmdd', 'yyyymm', ", 
                "'POSIXct' or 'Dates' class.")
        startdates <- Subset(Dates, along = ftime_dim, 1, drop = 'selected')
      }
    } else {
      stop("Parameter 'Dates' must have start date dimension and ", 
          "forecast time dimension.")
    }
  }
  # startdates
  if (is.null(startdates)) {
    if (is.null(sdate_dim)) {
      startdates <- 'XXX'
    } else {
      startdates <- rep('XXX', dim(data)[sdate_dim])
    }
  } else {
    if (is.null(sdate_dim)) {
      if (length(startdates) != 1) {
        warning("Parameter 'startdates' has length more than 1. Only first ", 
                "value will be used.")
        startdates <- startdates[[1]]
      }
    }
  }
  # Datasets
  if (is.null(Datasets)) {
    if (!single_file) {
      warning("Parameter 'Datasets' is NULL. Files will be saved with a ", 
              "directory name of 'XXX'.")
    }
    Datasets <- rep('XXX', n_datasets )
  }
  if (inherits(Datasets, 'list')) {
    Datasets <- names(Datasets)
  }
  if (n_datasets > length(Datasets)) {
    warning("Dimension 'Datasets' in 'data' is greater than those listed in ",
            "element 'Datasets' and the first element will be reused.")
    Datasets <- c(Datasets, rep(Datasets[1], n_datasets - length(Datasets)))
  } else if (n_datasets < length(Datasets)) {
    warning("Dimension 'Datasets' in 'data' is smaller than those listed in ",
            "element 'Datasets' and only the firsts elements will be used.")
    Datasets <- Datasets[1:n_datasets]
  }

  ## Unknown dimensions check
  alldims <- c(dat_dim, var_dim, sdate_dim, lon_dim, lat_dim, memb_dim, ftime_dim)
  if (!all(dimnames %in% alldims)) {
    unknown_dims <- dimnames[which(!dimnames %in% alldims)]
    warning("Detected unknown dimension: ", paste(unknown_dims, collapse = ', '))
    memb_dim <- c(memb_dim, unknown_dims)
    alldims <- c(dat_dim, var_dim, sdate_dim, lon_dim, lat_dim, memb_dim, ftime_dim)
  }
  # Reorder
  if (any(dimnames != alldims)) {
    data <- Reorder(data, alldims)
    dimnames <- names(dim(data))
    if (!is.null(attr(data, 'dimensions'))) {
      attr(data, 'dimensions') <- dimnames
    }
  }

  ## NetCDF dimensions definition
  defined_dims <- NULL
  extra_info_dim <- NULL
  if (is.null(Dates)) {
    filedims <- dimnames[which(!dimnames %in% c(dat_dim, var_dim))]
  } else {
    filedims <- dimnames[which(!dimnames %in% c(dat_dim, var_dim, sdate_dim, ftime_dim))]
  }
  for (i_coord in filedims) {
    dim_info <- list()
    # vals
    if (i_coord %in% names(coords)) {
      if (is.numeric(coords[[i_coord]])) {
        dim_info[['vals']] <- as.vector(coords[[i_coord]])
      } else {
        dim_info[['vals']] <- 1:dim(data)[i_coord]
      }
    } else {
      dim_info[['vals']] <- 1:dim(data)[i_coord]
    }
    # name
    dim_info[['name']] <- i_coord
    # len
    dim_info[['len']] <- as.numeric(dim(data)[i_coord])
    # unlim
    dim_info[['unlim']] <- FALSE
    # create_dimvar
    dim_info[['create_dimvar']] <- TRUE
    ## metadata
    if (i_coord %in% names(metadata)) {
      if ('variables' %in% names(attributes(metadata[[i_coord]]))) {
        # from Start: 'lon' or 'lat'
        attrs <- attributes(metadata[[i_coord]])[['variables']][[i_coord]]
        i_coord_info <- attrs[!sapply(attrs, inherits, 'list')]
      } else if (inherits(metadata[[i_coord]], 'list')) {
        # from Start and Load: main var
        i_coord_info <- metadata[[i_coord]]
      } else if (!is.null(attributes(metadata[[i_coord]]))) {
        # from Load
        i_coord_info <- attributes(metadata[[i_coord]])
      } else {
        stop("Metadata is not correct.")
      }
      # len
      if ('size' %in% names(i_coord_info)) {
        if (i_coord_info[['size']] != dim(data)[i_coord]) {
          dim_info[['original_len']] <- i_coord_info[['size']]
          i_coord_info[['size']] <- NULL
        }
      }
      # units
      if (!('units' %in% names(i_coord_info))) {
        dim_info[['units']] <- ''
      } else {
        dim_info[['units']] <- i_coord_info[['units']]
        i_coord_info[['units']] <- NULL
      }
      # calendar
      if (!('calendar' %in% names(i_coord_info))) {
        dim_info[['calendar']] <- NA
      } else {
        dim_info[['calendar']] <- i_coord_info[['calendar']]
        i_coord_info[['calendar']] <- NULL
      }
      # longname
      if ('long_name' %in% names(i_coord_info)) {
        dim_info[['longname']] <- i_coord_info[['long_name']]
        i_coord_info[['long_name']] <- NULL
      } else if ('longname' %in% names(i_coord_info)) {
        dim_info[['longname']] <- i_coord_info[['longname']]
        i_coord_info[['longname']] <- NULL
      } else {
        if (i_coord %in% .KnownLonNames()) {
          dim_info[['longname']] <- 'longitude'
        } else if (i_coord %in% .KnownLatNames()) {
          dim_info[['longname']] <- 'latitude'
        }
      }
      # extra information
      if (!is.null(names(i_coord_info))) {
        extra_info_dim[[i_coord]] <- i_coord_info
      }
    } else {
      # units
      dim_info[['units']] <- "adim"
      # longname
      dim_info[['longname']] <- i_coord
      # calendar
      dim_info[['calendar']] <- NA
    }
    new_dim <- list(ncdim_def(name = dim_info[['name']], units = dim_info[['units']], 
                              vals = dim_info[['vals']], unlim = dim_info[['unlim']], 
                              create_dimvar = dim_info[['create_dimvar']], 
                              calendar = dim_info[['calendar']], 
                              longname = dim_info[['longname']]))
    names(new_dim) <- i_coord
    defined_dims <- c(defined_dims, new_dim)
  }

  defined_vars <- list()
  if (!single_file) {
    for (i in 1:n_datasets) {
      path <- file.path(destination, Datasets[i], varname)
      for (j in 1:n_vars) {
        dir.create(path[j], recursive = TRUE)
        startdates <- gsub("-", "", startdates)
        dim(startdates) <- c(length(startdates))
        names(dim(startdates)) <- sdate_dim
        if (is.null(dat_dim) & is.null(var_dim)) {
          data_subset <- data
        } else if (is.null(dat_dim)) {
          data_subset <- Subset(data, c(var_dim), list(j), drop = 'selected')
        } else if (is.null(var_dim)) {
          data_subset <- Subset(data, along = c(dat_dim), list(i), drop = 'selected')
        } else {
          data_subset <- Subset(data, c(dat_dim, var_dim), list(i, j), drop = 'selected')
        }
        if (is.null(Dates)) {
          input_data <- list(data_subset, startdates)
          target_dims <- list(c(lon_dim, lat_dim, memb_dim, ftime_dim), NULL)
        } else {
          input_data <- list(data_subset, startdates, Dates)
          target_dims = list(c(lon_dim, lat_dim, memb_dim, ftime_dim), NULL, ftime_dim)
        }
        Apply(data = input_data,
              target_dims = target_dims,
              fun = .saveExp, 
              destination = path[j],
              defined_dims = defined_dims, 
              ftime_dim = ftime_dim, 
              varname = varname[j], 
              metadata_var = metadata[[varname[j]]], 
              extra_info_dim = extra_info_dim, 
              extra_string = extra_string)
      }
    }
  } else {
    # Datasets definition
    # From here
    if (!is.null(dat_dim)) {
      new_dim <- list(ncdim_def(name = dat_dim, units = "adim",
                                vals = 1 : dim(data)[dat_dim],
                                longname = 'Datasets', create_dimvar = TRUE))
      names(new_dim) <- dat_dim
      defined_dims <- c(new_dim, defined_dims)
      extra_info_dim[[dat_dim]] <- list(Datasets = paste(Datasets, collapse = ', '))
    }
    first_sdate <- last_sdate <- NULL
    if (!is.null(Dates)) {
      # sdate definition
      sdates <- Subset(Dates, along = ftime_dim, 1, drop = 'selected')
      differ <- as.numeric((sdates - sdates[1])/3600)
      new_dim <- list(ncdim_def(name = sdate_dim, units = paste('hours since', sdates[1]),
                                vals = differ,
                                longname = sdate_dim, create_dimvar = TRUE))
      names(new_dim) <- sdate_dim
      defined_dims <- c(defined_dims, new_dim)
      first_sdate <- sdates[1]
      last_sdate <- sdates[length(sdates)]
      # ftime definition
      Dates <- Reorder(Dates, c(ftime_dim, sdate_dim))
      differ_ftime <- apply(Dates, 2, function(x){as.numeric((x - x[1])/3600)})
      dim(differ_ftime) <- dim(Dates)
      differ_ftime_subset <- Subset(differ_ftime, along = sdate_dim, 1, drop = 'selected')
      if (all(apply(differ_ftime, 1, function(x){length(unique(x)) == 1}))) {
        if (all(diff(differ_ftime_subset/24) == 1)) {
          # daily values
          dim_time <- list(ncdim_def(name = ftime_dim, units = 'days',
                                     vals = round(differ_ftime_subset/24) + 1, 
                                     calendar = 'proleptic_gregorian',
                                     longname = ftime_dim, unlim = TRUE))
          names(dim_time) <- ftime_dim
          defined_dims <- c(defined_dims, dim_time)     
        } else if (all(diff(differ_ftime_subset/24) %in% c(28, 29, 30, 31))) {
          # monthly values
          dim_time <- list(ncdim_def(name = ftime_dim, units = 'months',
                                     vals = round(differ_ftime_subset/730) + 1, 
                                     calendar = 'proleptic_gregorian',
                                     longname = ftime_dim, unlim = TRUE))
          names(dim_time) <- ftime_dim
          defined_dims <- c(defined_dims, dim_time)   
        } else {
          # other frequency
          dim_time <- list(ncdim_def(name = ftime_dim, units = 'hours',
                                     vals = differ_ftime_subset + 1, 
                                     calendar = 'proleptic_gregorian',
                                     longname = ftime_dim, unlim = TRUE))
          names(dim_time) <- ftime_dim
          defined_dims <- c(defined_dims, dim_time)  
        }
      } else {
        warning("Time steps are not equal for all start dates. Only ", 
                "forecast time values for the first start date will be saved ", 
                "correctly.")
        dim_time <- list(ncdim_def(name = ftime_dim, 
                                   units = paste('hours since', 
                                           paste(sdates, collapse = ', ')),
                                   vals = differ_ftime_subset, 
                                   calendar = 'proleptic_gregorian',
                                   longname = ftime_dim, unlim = TRUE))
        names(dim_time) <- ftime_dim
        defined_dims <- c(defined_dims, dim_time)
      }
    }

    # var definition
    defined_vars <- list()
    extra_info_var <- NULL
    for (j in 1:n_vars) {
      var_info <- list()
      i_var_info <- metadata[[varname[j]]][!sapply(metadata[[varname[j]]], inherits, 'list')]
      ## Define metadata
      # name
      var_info[['name']] <- varname[j]
      # units
      if ('units' %in% names(i_var_info)) {
        var_info[['units']] <- i_var_info[['units']]
        i_var_info[['units']] <- NULL
      } else {
        var_info[['units']] <- ''
      }
      # dim
      var_info[['dim']] <- defined_dims
      # missval
      if ('missval' %in% names(i_var_info)) {
        var_info[['missval']] <- i_var_info[['missval']]
        i_var_info[['missval']] <- NULL
      } else {
        var_info[['missval']] <- NULL
      }
      # longname
      if (any(c('longname', 'long_name') %in% names(i_var_info))) {
        longname <- names(i_var_info)[which(names(i_var_info) %in% c('longname', 'long_name'))]
        var_info[['longname']] <- i_var_info[[longname]]
        i_var_info[[longname]] <- NULL
      } else {
        var_info[['longname']] <- varname[j]
      }
      # prec
      if ('prec' %in% names(i_var_info)) {
        var_info[['prec']] <- i_var_info[['prec']]
        i_var_info[['prec']] <- NULL
      } else {
        prec <- typeof(data)
        if (prec == 'character') {
          var_info[['prec']] <- 'char'
        }
        if (any(prec %in% c('short', 'float', 'double', 'integer', 'char', 'byte'))) {
          var_info[['prec']] <- prec
        } else {
          var_info[['prec']] <- 'double'
        }
      }
      # extra information
      if (!is.null(names(i_var_info))) {
        extra_info_var[[varname[j]]] <- i_var_info
      }
      new_var <- list(ncvar_def(name = var_info[['name']],
                                units = var_info[['units']],
                                dim = var_info[['dim']], 
                                missval = var_info[['missval']],
                                longname = var_info[['longname']], 
                                prec = var_info[['prec']]))

      names(new_var) <- varname[j]
      defined_vars <- c(defined_vars, new_var)
    }
    if (is.null(extra_string)) {
      gsub("-", "", first_sdate)
      file_name <- paste0(paste(c(varname, 
                                  gsub("-", "", first_sdate), 
                                  gsub("-", "", last_sdate)), 
                                  collapse = '_'), ".nc")
    } else {
      file_name <- paste0(paste(c(varname, extra_string, 
                                  gsub("-", "", first_sdate), 
                                  gsub("-", "", last_sdate)), 
                                  collapse = '_'), ".nc")
    }
    full_filename <- file.path(destination, file_name)
    file_nc <- nc_create(full_filename, defined_vars)
    if (is.null(var_dim)) {
      ncvar_put(file_nc, varname, vals = data)
    } else {
      for (j in 1:n_vars) {
        ncvar_put(file_nc, defined_vars[[j]]$name, 
                  vals = Subset(data, var_dim, j, drop = 'selected'))
      }
    }
    # Additional dimension attributes
    for (dim in names(defined_dims)) {
      if (dim %in% names(extra_info_dim)) {
        for (info_dim in names(extra_info_dim[[dim]])) {
          ncatt_put(file_nc, dim, info_dim, as.character(extra_info_dim[[dim]][[info_dim]]))
        }
      }
    }
    # Additional dimension attributes
    for (var in names(defined_vars)) {
      if (var %in% names(extra_info_var)) {
        for (info_var in names(extra_info_var[[var]])) {
          ncatt_put(file_nc, var, info_var, as.character(extra_info_var[[var]][[info_var]]))
        }
      }
    }
    nc_close(file_nc)
  }
}

.saveExp <- function(data, startdates = NULL, dates = NULL, destination = "./", 
                     defined_dims, ftime_dim = 'time',  varname = 'var', 
                     metadata_var = NULL, extra_info_dim = NULL, 
                     extra_string = NULL) {
  # ftime_dim
  if (!is.null(dates)) {
    differ <- as.numeric((dates - dates[1])/3600)
    dim_time <- list(ncdim_def(name = ftime_dim, units = paste('hours since', dates[1]),
                               vals = differ, calendar = 'proleptic_gregorian',
                               longname = ftime_dim, unlim = TRUE))
    names(dim_time) <- ftime_dim
    defined_dims <- c(defined_dims, dim_time)
  }

  ## Define var metadata
  var_info <- NULL
  extra_info_var <- NULL
  i_var_info <- metadata_var[!sapply(metadata_var, inherits, 'list')]

  # name
  var_info[['name']] <- varname
  # units
  if ('units' %in% names(i_var_info)) {
    var_info[['units']] <- i_var_info[['units']]
    i_var_info[['units']] <- NULL
  } else {
    var_info[['units']] <- ''
  }
  # dim
  var_info[['dim']] <- defined_dims
  # missval
  if ('missval' %in% names(i_var_info)) {
    var_info[['missval']] <- i_var_info[['missval']]
    i_var_info[['missval']] <- NULL
  } else {
    var_info[['missval']] <- NULL
  }
  # longname
  if (any(c('longname', 'long_name') %in% names(i_var_info))) {
    longname <- names(i_var_info)[which(names(i_var_info) %in% c('longname', 'long_name'))]
    var_info[['longname']] <- i_var_info[[longname]]
    i_var_info[[longname]] <- NULL
  } else {
    var_info[['longname']] <- varname
  }
  # prec
  if ('prec' %in% names(i_var_info)) {
    var_info[['prec']] <- i_var_info[['prec']]
    i_var_info[['prec']] <- NULL
  } else {
    prec <- typeof(data)
    if (prec == 'character') {
      var_info[['prec']] <- 'char'
    }
    if (any(prec %in% c('short', 'float', 'double', 'integer', 'char', 'byte'))) {
      var_info[['prec']] <- prec
    } else {
      var_info[['prec']] <- 'double'
    }
  }
  # extra information
  if (!is.null(names(i_var_info))) {
    extra_info_var <- i_var_info
  }

  datanc <- ncvar_def(name = var_info[['name']],
                      units = var_info[['units']],
                      dim = var_info[['dim']], 
                      missval = var_info[['missval']],
                      longname = var_info[['longname']], 
                      prec = var_info[['prec']])

  if (is.null(extra_string)) {
    file_name <- paste0(varname, "_", startdates, ".nc")
  } else {
    file_name <- paste0(varname, "_", extra_string, "_", startdates, ".nc")
  }
  full_filename <- file.path(destination, file_name)
  file_nc <- nc_create(full_filename, datanc)
  ncvar_put(file_nc, datanc, data)

  # Additional attributes
  for (dim in names(defined_dims)) {
    if (dim %in% names(extra_info_dim)) {
      for (info_dim in names(extra_info_dim[[dim]])) {
        ncatt_put(file_nc, dim, info_dim, as.character(extra_info_dim[[dim]][[info_dim]]))
      }
    }
  }
  # Additional dimension attributes
  if (!is.null(extra_info_var)) {
    for (info_var in names(extra_info_var)) {
      ncatt_put(file_nc, varname, info_var, as.character(extra_info_var[[info_var]]))
    }
  }

  nc_close(file_nc)
}
