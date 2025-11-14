#' Generate a Summary of the data and metadata in the s2dv_cube object
#'
#'This function prints the summary of the data and metadata of an object of 
#'class \code{s2dv_cube}.
#'   
#' @author Theertha Kariyathan, \email{theertha.kariyathan@bsc.es}
#'
#' @param data An \code{s2dv_cube} object containing:
#'   - \code{data}: N-dimensional array with named dimensions
#'   - \code{dim}: Dimensions, including \code{var} (variables).
#'   - \code{attrs}: Attributes such as \code{VarName} and \code{Metadata}.
#'   - \code{coords}: Named list with coordinates of dimensions.
#' @param show_NA A logical value indicating if details of NA values in the
#'   loaded object will be displayed in the output or not. Default = FALSE.
#' @param show_loaded_files A logical value indicating if the names of the 
#'   loaded files will be displayed in the output or not. Default = FALSE.
#' @param var_dim A character string indicating the name of the variable 
#'   dimension. Default = "var".
#' @return Does not return a value. This function prints a detailed summary of 
#' an \code{s2dv_cube} object to the console.
#' @details The function uses the data and metadata from the loaded 
#'   \code{s2dv_cube} object to generate a summary of the object.The summary    
#'   includes :
#'   - months: Months that have been loaded.
#'   - range: Range of the dates that have been loaded.
#'   - dimensions: Object dimensions.
#'   - Statistical summary of the data in data: Basic statistical 
#'     summary of the data.
#'   - Variable: Loaded Variables, along with their units (units:)
#'   - NA-Indices per Dimension: Index with NA values per dimension
#'   - Number of NAs in NA-Indices per Dimensions: Number of NAs, 
#'     in the Indices with NA values per dimension  
#'   - Loaded files: Successfully loaded Files
#'    
#' @examples
#' # Example 1:
#' CST_Summary(data = lonlat_temp_st$exp)
#' 
#' # Example 2:
#' \dontrun{
#' # s2dv cube paths
#' repos1 <- "/esarchive/exp/ecmwf/system4_m1/monthly_mean/$var$_f6h/$var$_$sdate$.nc"
#' repos2 <- "/esarchive/exp/ecmwf/system5_m1/monthly_mean/$var$_f6h/$var$_$sdate$.nc"
#'
#' # Create data cube
#' data <- CST_Start(dat = list(
#'                   list(name = 'system4_m1', path = repos1),
#'                   list(name = 'system5_m1', path = repos2)),
#'                   var = c('tas', 'sfcWind'),
#'                   sdate = '20170101',
#'                   ensemble = indices(1),
#'                   time = indices(1:3),
#'                   lat = indices(1:5),
#'                   lon = indices(1:5),
#'                   synonims = list(lat = c('lat', 'latitude'),
#'                                   lon = c('lon', 'longitude')),
#'                   return_vars = list(time = 'sdate',
#'                                      longitude = 'dat',
#'                                      latitude = 'dat'),
#'                   metadata_dims = c('dat', 'var'),
#'                   retrieve = TRUE)
#'
#' # Generate summary
#' CST_Summary(data)
#' }
#'
#' @seealso \link[CSTools]{CST_Start} or \link[CSTools]{s2dv_cube} for creating
#'   an s2dv cube object.
#' @importFrom utils capture.output
#' @importFrom dplyr %>%
#' @export

CST_Summary <- function(data, show_loaded_files = FALSE, show_NA = FALSE, 
                        var_dim = "var") {
  # Check 's2dv_cube'
  if (!inherits(data, "s2dv_cube")) {
    stop("Parameter 'data' must be of the class 's2dv_cube'.")
  }
  if (!is.logical(show_loaded_files)) {
    stop("Parameter 'show_loaded_files' must be logical.")
  }
  if (!is.logical(show_NA)) {
    stop("Parameter 'show_NA' must be logical.")
  }
  
  if (!is.character(var_dim) || length(var_dim) != 1) {
    stop("'var_dim' must be a single character string.")
  }
  
  if (!(var_dim %in% names(data$dims))) {
    warning(paste("Dimension", var_dim, "not found. Summary is not be split by variable."))
  }
  
  # Get name, leadtime months and date range
  object_name <- deparse(substitute(data)) 
  date_format <- "%b %d %Y"
  months <- unique(format(as.Date(data$attrs$Dates), format = "%B"))
  months <- paste(as.character(months), collapse = ", ")
  sdate_min <- format(min(as.Date(data$attrs$Dates), na.rm = TRUE),
                      format = date_format)
  sdate_max <- format(max(as.Date(data$attrs$Dates), na.rm = TRUE),
                      format = date_format)
  # Log the summary
  cat("DATA SUMMARY:\n")
  
  cat(paste(object_name, "months:", months), "\n")
  cat(paste(object_name, "range:", sdate_min, "to", sdate_max), "\n")
  cat(paste(object_name, "dimensions:"), "\n")
  
  # Use capture.output() and for loop to display results neatly
  output_string <- capture.output(dim(data$data))
  for (i in output_string) {
    cat(i, "\n")
  }
  
  # Print statistical summary of the data for every variable
  cat(paste("Statistical summary of the data in ", object_name, ":"), "\n")
  
  if (!(var_dim %in% names(data$dims))) {
    output_string <- capture.output(summary(data$data))
    for (i in output_string) {
      cat(i, "\n")
    }    
    
  } else {
    for (var_index in 1:data$dims[[var_dim]]) {
      variable_name <- data$attrs$Variable$varName[var_index]
      variable_units <- data$attrs$Variable$metadata[[variable_name]]$units
      cat(paste0("Variable: ", variable_name,
                 " (units: ", variable_units, ")"), "\n")
      output_string <- capture.output(summary(Subset(data$data,
                                                     along = var_dim,
                                                     indices = list(var_index))))
      for (i in output_string) {
        cat(i, "\n")
      }
    }
  }
  
  if (show_NA) {
    # Number of NAs per time dimension and latitude/longitude dimensions 
    list_na <- lapply(seq_along(dim(data$data)), function(dim) {
      apply(data$data, dim, function(x) sum(is.na(x)))
    })
    
    # Identify dimensions with NAs
    na_list <- sapply(list_na, function(x) which(x != 0) %>% paste(collapse = ",")) %>% unlist()
    names(na_list) <- names(dim(data$data))
    
    # Count the number of NAs per identified dimension
    num_nas <- sapply(list_na, function(x) x[which(x != 0)] %>% paste(collapse = ",")) %>% unlist()
    names(num_nas) <- names(dim(data$data))
    
    # Generate output strings
    output_na_list <- paste(names(na_list), na_list, sep = ": ", collapse = " ")
    output_num_nas <- paste(names(num_nas), num_nas, sep = ": ", collapse = " ")
    
    cat("NA-Indices per Dimension\n")
    cat(output_na_list, "\n")
    
    cat("Number of NAs in NA-Indices per Dimensions\n") 
    cat(output_num_nas, "\n") 
  }
  
  # Loaded files
  if (show_loaded_files) {
    all_files <- lapply(data$attrs$source_files, unlist) %>% unlist()
    loaded_files <- all_files[!is.na(all_files)]
    cat("Loaded files:\n")
    cat(paste(loaded_files, collapse = "\n"), "\n")
  }
  cat("---------------------------------------------", "\n")
}

