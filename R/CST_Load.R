#' CSTools Data Retreival Function
#'
#' This function aggregates, subsets and retrieves sub-seasonal, seasonal, decadal or climate projection data from NetCDF files in a local file system or on remote OPeNDAP servers, and arranges it for easy application of the CSTools functions.
#'
#' It receives any number of parameters (`...`) that are automatically forwarded to the `s2dverification::Load` function. See details in `?s2dverification::Load`.
#'
#' It is recommended to use this function in combination with the `zeallot::"%<-%"` operator, to directly assing the two returned 's2dv_cube's to two separate variables, which can then be sent independently to other functions in CSTools as needed. E.g.: `c(exp, obs) <- CST_Load(...)`.
#'
#' @param ... Parameters that are automatically forwarded to the `s2dverification::Load` function. See details in `?s2dverification::Load`.
#' @return A list with one or two S3 objects, named 'exp' and 'obs', of the class 's2dv_cube', containing experimental and date-corresponding observational data, respectively. These 's2dv_cube's can be ingested by other functions in CSTools. If the parameter `exp` in the call to `CST_Load` is set to `NULL`, then only the 'obs' component is returned, and viceversa.
#' @author Nicolau Manubens, \email{nicolau.manubens@bsc.es}
#' @import s2dverification
#' @importFrom utils glob2rx
#' @export
#' @examples
#' \dontrun{
#' library(zeallot)
#' startDates <- c('20001101', '20011101', '20021101', 
#'                 '20031101', '20041101', '20051101')
#' c(exp, obs) %<-% 
#'   CST_Load(
#'     var = 'tas', 
#'     exp = 'system5c3s', 
#'     obs = 'era5', 
#'     nmember = 15,
#'     sdates = startDates,
#'     leadtimemax = 3,
#'     latmin = 27, latmax = 48,
#'     lonmin = -12, lonmax = 40,
#'     output = 'lonlat',
#'     nprocs = 1
#'   )
#' }
#' \dontshow{
#' exp <- CSTools::lonlat_data$exp
#' obs <- CSTools::lonlat_data$obs
#' }
CST_Load <- function(...) {
  exp <- Load(...)

  if (is.null(exp) || (is.null(exp$mod) && is.null(exp$obs))) {
    stop("The s2dverification::Load call did not return any data.")
  }

  obs <- exp
  obs$mod <- NULL
  exp$obs <- NULL
  names(exp)[[1]] <- 'data'
  names(obs)[[1]] <- 'data'

  remove_matches <- function(v, patterns) {
    if (length(v) > 0) {
      matches <- c()
      for (pattern in patterns) {
        matches <- c(matches, which(grepl(pattern, v)))
      }
      if (length(matches) > 0) {
        v <- v[-matches]
      }
    }
    v
  }

  harmonize_patterns <- function(v) {
    matches <- grepl('.*\\.nc$', v)
    if (sum(!matches) > 0) {
      match_indices <- which(!matches)
      v[match_indices] <- sapply(v[match_indices], function(x) paste0(x, '*'))
    }
    v <- glob2rx(v)
    v <- gsub('\\$.*\\$', '*', v)
    v
  }

  if (!is.null(obs$data)) {
    obs$Datasets$exp <- NULL
    obs$Datasets <- obs$Datasets$obs
    obs_path_patterns <- sapply(obs$Datasets, function(x) attr(x, 'source'))
    obs_path_patterns <- harmonize_patterns(obs_path_patterns)
  }

  if (!is.null(exp$data)) {
    exp$Datasets$obs <- NULL
    exp$Datasets <- exp$Datasets$exp
    exp_path_patterns <- sapply(exp$Datasets, function(x) attr(x, 'source'))
    exp_path_patterns <- harmonize_patterns(exp_path_patterns)
  }

  if (!is.null(obs$data) && !is.null(exp$data)) {
    obs$source_files <- remove_matches(obs$source_files,
                                       exp_path_patterns)
    obs$not_found_files <- remove_matches(obs$not_found_files,
                                          exp_path_patterns)

    exp$source_files <- remove_matches(exp$source_files,
                                       obs_path_patterns)
    exp$not_found_files <- remove_matches(exp$not_found_files,
                                          obs_path_patterns)
  }

  result <- list()
  if (!is.null(exp$data)) {
    class(exp) <- 's2dv_cube'
    result$exp <- exp
  }
  if (!is.null(obs$data)) {
    class(obs) <- 's2dv_cube'
    result$obs <- obs
  }
  result
}
