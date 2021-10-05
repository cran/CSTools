#'Save CSTools objects of class 's2dv_cube' containing experiments or observed 
#'data in NetCDF format
#'
#'@author Perez-Zanon Nuria, \email{nuria.perez@bsc.es}
#'
#'@description This function allows to divide and save a object of class 
#''s2dv_cube' into a NetCDF file, allowing to reload the saved data using 
#'\code{CST_Load} function.
#'
#'@param data an object of class \code{s2dv_cube}.
#'@param destination a character string containing the directory name in which 
#'to save the data. NetCDF file for each starting date are saved into the 
#'folder tree: destination/experiment/variable/. By default the function 
#'creates and saves the data into the folder "CST_Data" in the working 
#'directory.
#'@param extra_string a character string to be include as part of the file name, for instance, to identify member or realization. It would be added to the file name between underscore characters.
#'
#'@seealso \code{\link{CST_Load}}, \code{\link{as.s2dv_cube}} and \code{\link{s2dv_cube}}
#'
#'@import ncdf4
#'@importFrom s2dv Reorder InsertDim
#'@import multiApply
#'
#'@examples
#'\dontrun{
#'library(CSTools)
#'data <- lonlat_data$exp
#'destination <- "./path2/"
#'CST_SaveExp(data = data, destination = destination)
#'}
#'
#'@export
CST_SaveExp <- function(data, destination = "./CST_Data", extra_string = NULL) {
  if (!is.character(destination) & length(destination) > 1) {
    stop("Parameter 'destination' must be a character string of one element ",
         "indicating the name of the file (including the folder if needed) ",
         "where the data will be saved.")
  }
  if (!inherits(data, 's2dv_cube')) {
    stop("Parameter 'data' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load.")
  }
  sdates <- lapply(1:length(data$Datasets), function(x) {
                          unique(data$Datasets[[x]]$InitializationDates)})[[1]]
  if (!is.character(attributes(data$Variable)$units)) {
      units <- attributes(data$Variable)$variable$units
  } else {
      units <- attributes(data$Variable)$units
  }
  cdo_grid_name = attr(data$lon, 'cdo_grid_name')
  projection = attr(data$lon, 'projection')
  var_name <- data$Variable$varName
  time_values <- data$Dates$start
  dim(time_values) <- c(time = length(time_values) / length(sdates), 
                        sdate = length(sdates)) 
  SaveExp(data = data$data, lon = data$lon, lat = data$lat,
          Dataset = names(data$Datasets), var_name = var_name,
          units = units, cdo_grid_name = cdo_grid_name, projection = projection,
          startdates = sdates, Dates = time_values, destination, 
          extra_string = extra_string) 
} 
#'Save an experiment in a format compatible with CST_Load
#'@description This function is created for compatibility with CST_Load/Load for saving post-processed datasets such as those calibrated of downscaled with CSTools functions
#'
#'@author Perez-Zanon Nuria, \email{nuria.perez@bsc.es}
#'
#'@param data an multi-dimensional array with named dimensions (longitude, latitude, time, member, sdate)
#'@param lon vector of logitud corresponding to the longitudinal dimension in data
#'@param lat vector of latitud corresponding to the latitudinal dimension in data
#'@param Dataset a vector of character string indicating the names of the datasets
#'@param var_name a character string indicating the name of the variable to be saved
#'@param units a character string indicating the units of the variable
#'@param startdates a vector of dates indicating the initialization date of each simulations
#'@param Dates a matrix of dates with two dimension 'time' and 'sdate'.
#'@param cdo_grid_name a character string indicating the name of the grid e.g.: 'r360x181'
#'@param projection a character string indicating the projection name
#'@param destination a character string indicating the path where to store the NetCDF files
#'@param extra_string a character string to be include as part of the file name, for instance, to identify member or realization.
#'
#'@return the function creates as many files as sdates per dataset. Each file could contain multiple members. It would be added to the file name between underscore characters.
#' The path will be created with the name of the variable and each Datasets.
#' 
#'@import ncdf4
#'@importFrom s2dv Reorder InsertDim
#'@import multiApply
#'
#'@examples
#'\dontrun{
#'data <- lonlat_data$exp$data
#'lon <- lonlat_data$exp$lon
#'lat <- lonlat_data$exp$lat
#'Dataset <- 'XXX'
#'var_name <- 'tas'
#'units <- 'k'
#'startdates <- lapply(1:length(lonlat_data$exp$Datasets),
#'                     function(x) {
#'                         lonlat_data$exp$Datasets[[x]]$InitializationDates[[1]]})[[1]]
#'Dates <- lonlat_data$exp$Dates$start
#'dim(Dates) <- c(time = length(Dates)/length(startdates), sdate = length(startdates))
#'cdo_grid_name = attr(lonlat_data$exp$lon, 'cdo_grid_name')
#'projection = attr(lonlat_data$exp$lon, 'projection')
#'destination = './path/'
#'SaveExp(data, lon, lat, Dataset, var_name, units, startdates, Dates,
#'                    cdo_grid_name, projection, destination)
#'}
#'@export
SaveExp <- function(data, lon, lat, Dataset, var_name, units, startdates, Dates,
                    cdo_grid_name, projection, destination, 
                    extra_string = NULL) {  
  dimname <- names(dim(data))
  if (any(dimname == "ftime")) {
    dimname[which(dimname == "ftime")] <- "time"
    names(dim(data))[which(dimname == "ftime")] <- "time"
  }
  if (any(dimname == "memb")) {
    dimname[which(dimname == "memb")] <- "member"
    names(dim(data))[which(dimname == "memb")] <- "member"
  }
  if (any(dimname == "ensemble")) {
    dimname[which(dimname == "ensemble")] <- "member"
    names(dim(data))[which(dimname == "ensemble")] <- "member"
  }
  if (any(dimname == "lon")) {
    dimname[which(dimname == "lon")] <- "longitude"
    names(dim(data))[which(dimname == "lon")] <- "longitude"
  }
  if (any(dimname == "lat")) {
    dimname[which(dimname == "lat")] <- "latitude"
    names(dim(data))[which(dimname == "lat")] <- "latitude"
  }
  names(dim(data)) <- dimname
  if (is.null(dimname)) {
    stop("Element 'data' in parameter 'data' must have named dimensions.")
  }
  sdate_pos <- which(dimname == "sdate")

  if (length(sdate_pos) == 0) {
    stop("Element 'data' in parameter 'data' hasn't 'sdate' dimension.")
  } else if (length(sdate_pos) > 1) {
    stop("Element 'data' in parameter 'data' has more than one 'sdate'",
         " dimension.")
  }
  if (!is.null(extra_string)) {
    if (!is.character(extra_string)) {
      stop("Parameter 'extra_string' must be a character string.")
    }
  }
  dataset_pos <- which(dimname == "dataset" | dimname == "dat")
  dims <- dim(data)
  if (length(dataset_pos) == 0) {
    warning("Element 'data' in parameter 'data' hasn't 'dataset' dimension. ",
            "All data is stored in the same 'dataset' folder.")
    data$data <- InsertDim(data, posdim = 1, lendim = 1)
    names(dim(data))[1] <- "dataset"
    dimname <- c("dataset", dimname)
    dataset_pos = 1
  } else if (length(dataset_pos) > 1) {
    stop("Element 'data' in parameter 'data' has more than one 'dataset'",
         " dimension.")
  }
  n_datasets <- dim(data)[dataset_pos] # number of folder by dataset
  # dataset names:
  datasets <- Dataset
  if (n_datasets > length(datasets)) {
    warning("Dimension 'dataset' in element 'data' from parameter 'data' ",
            "is greater than those listed in element 'Datasets' and the ",
            "first element is reused.")
    datasets <- c(datasets, rep(datasets[1], n_datasets - length(datasets)))
  } else if (n_datasets < length(datasets)) {
    warning("Dimension 'dataset' in element 'data' from parameter 'data', ",
            "is smaller than those listed in element 'Datasets' and only the",
            " first element will be used.")
    datasets <- datasets[1 : n_datasets]
  }
  # var names:
  if ('var' %in% dimname) {
    var_pos <- which(dimname == 'var')
      if (dims[var_pos] == 1) {
          data <- adrop(data, drop = var_pos)
          dimname <- names(dim(data))
      }
  }
  if (length(var_name) != 1) {
    stop("One variable name must be included in element 'Variable$varName' ",
         "of parameter 'data'.") 
  }
  if (!is.character(var_name)) {
    stop("Element 'Variable$varName' of parameter 'data' ",
         "must be a character string.")
  }

  known_dim_names <- c("var", "lat", "latitude", "lon", "longitude", "time", 
                       "ftime", "sdate", "dataset", "dat", "nlevel", "levels")
  dims_var <- NULL
  list_pos <- 1

  if (any(dimname == 'longitude') | any(dimname == 'lon')) {
    dim_lon <- ncdim_def(name = 'lon', units = 'degrees',
                         vals = as.vector(lon), longname = 'longitude')
    dims_var[[list_pos]] <- dim_lon
    list_pos <- list_pos + 1
  }
    if (any(dimname == 'latitude') | any(dimname == 'lat')) {
    dim_lat <- ncdim_def(name = 'lat', units = 'degrees_north',
                         vals = as.vector(lat), longname = 'latitude')
    dims_var[[list_pos]] <- dim_lat
    list_pos <- list_pos + 1
  }
  if (any(!(dimname %in% known_dim_names))) {
    dims_member <- dimname[!(dimname %in% known_dim_names)]
    if (length(dims_member) > 1) {
      stop("Ask for saving realizations or further dimensions to the mantainer.")
    } else {
      dim_memb <- ncdim_def(name = 'ensemble', units = "adim",
                          vals = 1 : dim(data)[which(dimname == 'member')],
                          longname = 'ensemble', create_dimvar = TRUE)
    dims_var[[list_pos]] <- dim_memb
    list_pos <- list_pos + 1
    }
  }


  if (any(dimname == 'level')) {
    stop("Ask for saving 3Dim fields to the mantainer.")
  }

  for (i in 1 : n_datasets) {
    path <- file.path(destination, datasets[i], var_name)
    dir.create(path, recursive = TRUE)
    startdate <- gsub("-", "", startdates)

    dim(startdate) <- c(sdate = length(startdate)) 
    Apply(list(data, startdate, Dates),
          target_dims = list(c('member', 'time', 'latitude', 'longitude'), 
                             NULL, 'time'),
          fun = .saveExp, var_name = var_name, units = units,
          dims_var = dims_var, cdo_grid_name = cdo_grid_name, projection = projection,
          destination = path, extra_string = extra_string)
  }
}

# data is an array with dimensions: member, time, lat, lon:
# Dates is a vector of the dates for the time dimension
# dims_var is a list with the ncdim_def of common variables in dataset: member, lat and lon:
# data <- 1:(3 * 4 * 5 * 6)
# dim(data) <- c(longitude = 3, latitude = 4, time = 5, member = 6)
# var_name <- 'tas'
# units <- 'K'
# lon <- 1:3
# lat <- 1:4
# sdate = '19001101'
# destination = '/esarchive/scratch/nperez/git/Flor/cstools/'
# dims_var = list(ncdim_def(name = 'lon', units = 'degrees',
#                          vals = as.vector(lon), longname = 'longitude'),
#                ncdim_def(name = 'lat', units = 'degrees_north',
#                         vals = as.vector(lat), longname = 'latitude'),
#                ncdim_def(name = 'ensemble', units = "adim",
#                          vals = 1 : 6,
#                          longname = 'ensemble', create_dimvar = TRUE))
#Dates <- as.Date(c("1900-11-01", "1900-12-01", "1901-01-01", "1901-02-01", "1901-03-01"))
#.saveExp(data, sdate, Dates, var_name, units, dims_var, cdo_grid_name = 'r360x181', projection = 'none', destination)
.saveExp <- function(data, sdate, Dates, var_name, units, dims_var, 
                     cdo_grid_name, projection, destination, extra_string) {
  dim_names <- names(dim(data))
  if (any(dim_names != c('longitude', 'latitude', 'member', 'time'))) {
    data <- Reorder(data, c('longitude', 'latitude', 'member', 'time'))
  }
  differ <- as.numeric((Dates - Dates[1])/3600)  
  dim_time <- ncdim_def(name = 'time', units = paste('hours since', Dates[1]),
                        vals = differ, calendar = 'proleptic_gregorian',
                        longname = 'time', unlim = TRUE)
  list_pos = length(dims_var) + 1
  dims_var[[list_pos]] <- dim_time
  datanc  <- ncvar_def(name = var_name,
                       units = units,
                       dim = dims_var, missval = -99999)
  if (is.null(extra_string)) {
    file_name <- paste0(var_name, "_", sdate, ".nc")
  } else {
    file_name <- paste0(var_name, "_", extra_string, "_", sdate, ".nc")
  }
  full_filename <- file.path(destination, file_name)
  file_nc <- nc_create(full_filename, datanc)
  ncvar_put(file_nc, datanc, data)
  ncatt_put(file_nc, datanc, 'coordinates', cdo_grid_name)
  ncatt_put(file_nc, datanc, 'projection', projection)
  nc_close(file_nc)
}
