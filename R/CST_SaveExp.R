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
#'
#'@seealso \code{\link{CST_Load}}, \code{\link{as.s2dv_cube}} and \code{\link{s2dv_cube}}
#'
#'@import s2dverification
#'@import ncdf4
#'@import abind
#'
#'@examples
#'\dontrun{
#'library(CSTools)
#'data <- lonlat_data$exp
#'destination <- "./path/"
#'CST_SaveExp(data = data, destination = destination)
#'}
#'
#'@export
CST_SaveExp <- function(data, destination = "./CST_Data") {
  if (!is.character(destination) & length(destination) > 1) {
    stop("Parameter 'destination' must be a character string of one element ",
         "indicating the name of the file (including the folder if needed) ",
         "where the data will be saved.")
  }
  if (!inherits(data, 's2dv_cube')) {
    stop("Parameter 'data' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load.")
  }
dimname <- names(dim(data$data))
  if (any(dimname == "time")) {
    dimname[which(dimname == "time")] <- "ftime"
    names(dim(data$data))[which(dimname == "time")] <- "ftime"
  }
  if (any(dimname == "memb")) {
    dimname[which(dimname == "memb")] <- "member"
    names(dim(data$data))[which(dimname == "memb")] <- "member"
  }
  if (any(dimname == "ensemble")) {
    dimname[which(dimname == "ensemble")] <- "member"
    names(dim(data$data))[which(dimname == "ensemble")] <- "member"
  }
  if (any(dimname == "longitude")) {
    dimname[which(dimname == "longitude")] <- "lon"
    names(dim(data$data))[which(dimname == "longitude")] <- "lon"
  }
  if (any(dimname == "latitude")) {
    dimname[which(dimname == "latitude")] <- "lat"
    names(dim(data$data))[which(dimname == "latitude")] <- "lat"
  }
names(dim(data$data)) <- dimname

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
n_sdates <- dim(data$data)[sdate_pos] # number of files to create

dataset_pos <- which(dimname == "dataset" | dimname == "dat")
dims <- dim(data$data)
  if (length(dataset_pos) == 0) {
    warning("Element 'data' in parameter 'data' hasn't 'dataset' dimension. ",
            "All data is stored in the same 'dataset' folder.")
    data$data <- InsertDim(var = data$data, posdim = 1, lendim = 1)
    names(dim(data$data))[1] <- "dataset"
    dimname <- c("dataset", dimname)
    dataset_pos = 1
  } else if (length(dataset_pos) > 1) {
    stop("Element 'data' in parameter 'data' has more than one 'dataset'",
         " dimension.")
  }
n_datasets <- dim(data$data)[dataset_pos] # number of folder by dataset
# dataset names:
datasets <- names(data$Datasets)
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
          data$data <- adrop(data$data, drop = var_pos)
          dimname <- names(dim(data$data))
      }
  }
var_name <- data$Variable$varName
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
                         vals = as.vector(data$lon), longname = 'longitude')
    dims_var[[list_pos]] <- dim_lon
    list_pos <- list_pos + 1
  }
    if (any(dimname == 'latitude') | any(dimname == 'lat')) {
    dim_lat <- ncdim_def(name = 'lat', units = 'degrees_north',
                         vals = as.vector(data$lat), longname = 'latitude')
    dims_var[[list_pos]] <- dim_lat
    list_pos <- list_pos + 1
  }
  if (any(!(dimname %in% known_dim_names))) {
    dims_member <- dimname[!(dimname %in% known_dim_names)]
    if (length(dims_member) > 1) {
      stop("Ask for saving realizations or further dimensions to the mantainer.")
    } else {
      dim_memb <- ncdim_def(name = 'ensemble', units = "adim",
                          vals = 1 : dim(data$data)[which(dimname == 'member')],
                          longname = 'ensemble', create_dimvar = TRUE)
    dims_var[[list_pos]] <- dim_memb
    list_pos <- list_pos + 1
    }
  }

  # Lead-time depends on the start date
  nlt <- length(data$Dates$start)/n_sdates

  if (any(dimname == 'level')) {
    stop("Ask for saving 3Dim fields to the mantainer.")
  }

  for (i in 1 : n_datasets) {
    path <- file.path(destination, datasets[i], var_name)
    dir.create(path, recursive = TRUE)
    startdate <- gsub("-", "", data$Datasets[[i]]$InitializationDates[[1]])
    file_name <- paste0(var_name, "_", startdate, ".nc")
    full_filename <- file.path(path, file_name)

    data_dataset <- Subset(data$data, along = which(dimname == 'dataset'), indices = i)
    standard_order <- c("lon", "lat", "member", "ftime")
    change_names <- c("lon", "lat", "ensemble", "ftime")
      for (j in 1 : n_sdates) {
        n_data <- s2dverification::Subset(data_dataset,
                                          along = which(dimname == 'sdate'),
                                          indices = j, drop = TRUE)
        pos_standard_order <- match( standard_order, names(dim(n_data)))
        n_data <- aperm(n_data, pos_standard_order)

        names(dim(n_data)) <- change_names

    # Lead-time depends on the start date
    # The correct time should be selected from $Dates$start
    time_values <- as.Date(substr(data$Dates$start[(j * nlt - nlt + 1):(j * nlt)],
                           1, 10))

    if (any(dimname == 'time') | any(dimname == 'ftime')) {
        dim_time <- ncdim_def(name = 'time', units = 'days since 1970-01-01',
                              vals = as.numeric(time_values),
                              longname = 'time', unlim = TRUE)
          if (i == 1 & j == 1) {
             dims_var[[list_pos]] <- dim_time
             list_pos <- list_pos + 1
          }
    }
    if (!is.character(attributes(data$Variable)$units)) {
        units = attributes(data$Variable)$variable$units
    } else {
        units = attributes(data$Variable)$units
    }
        datanc  <- ncvar_def(name = var_name,
                             units = units,
                             dim = dims_var, missval = -99999)
        file_nc <- nc_create(full_filename[j], datanc)
        ncvar_put(file_nc, datanc, n_data)
        ncatt_put(file_nc, datanc, 'coordinates', attr(data$lon, 'cdo_grid_name'))
        ncatt_put(file_nc, datanc, 'projection', attr(data$lon, 'projection'))
        nc_close(file_nc)
      }
  }
}
