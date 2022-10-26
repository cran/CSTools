#'Conversion of 'startR_array' or 'list' objects to 's2dv_cube'
#'
#'This function converts data loaded using startR package or s2dv Load function into a 's2dv_cube' object.
#'
#'@author Perez-Zanon Nuria, \email{nuria.perez@bsc.es}
#'@author Nicolau Manubens, \email{nicolau.manubens@bsc.es}
#'
#'@param object an object of class 'startR_array' generated from function \code{Start} from startR package (version 0.1.3 from earth.bsc.es/gitlab/es/startR) or a list output from function \code{Load} from s2dv package.
#'
#'@return The function returns a 's2dv_cube' object to be easily used with functions \code{CST} from CSTools package.
#'
#'@seealso \code{\link{s2dv_cube}}, \code{\link[s2dv]{Load}}, \code{\link[startR]{Start}} and \code{\link{CST_Load}}
#'@examples
#'\dontrun{
#'library(startR)
#'repos <- '/esarchive/exp/ecmwf/system5_m1/monthly_mean/$var$_f6h/$var$_$sdate$.nc'
#'data <- Start(dat = repos,
#'              var = 'tas',
#'              sdate = c('20170101', '20180101'),
#'              ensemble = indices(1:20),
#'              time = 'all',
#'              latitude = 'all',
#'              longitude = indices(1:40),
#'              return_vars = list(latitude = 'dat', longitude = 'dat', time = 'sdate'),
#'              retrieve = TRUE)
#'data <- as.s2dv_cube(data)
#'class(data)
#'startDates <- c('20001101', '20011101', '20021101',
#'                 '20031101', '20041101', '20051101')
#'data <- Load(var = 'tas', exp = 'system5c3s', 
#'             nmember = 15, sdates = startDates,
#'             leadtimemax = 3, latmin = 27, latmax = 48,
#'             lonmin = -12, lonmax = 40, output = 'lonlat')
#'data <- as.s2dv_cube(data)
#'class(data)
#'}
#'@export
as.s2dv_cube <- function(object) {
  if (is.list(object)) {
    if (is.null(object) || (is.null(object$mod) && is.null(object$obs))) {
      stop("The s2dv::Load call did not return any data.")
    }
    obs <- object
    obs$mod <- NULL
    object$obs <- NULL
    names(object)[[1]] <- 'data'
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

    if (!is.null(object$data)) {
      object$Datasets$obs <- NULL
      object$Datasets <- object$Datasets$exp
      exp_path_patterns <- sapply(object$Datasets, function(x) attr(x, 'source'))
      exp_path_patterns <- harmonize_patterns(exp_path_patterns)
    }

    if (!is.null(obs$data) && !is.null(object$data)) {
      obs$source_files <- remove_matches(obs$source_files,
                                         exp_path_patterns)
      obs$not_found_files <- remove_matches(obs$not_found_files,
                                            exp_path_patterns)
  
      object$source_files <- remove_matches(object$source_files,
                                            obs_path_patterns)
      object$not_found_files <- remove_matches(object$not_found_files,
                                               obs_path_patterns)
    }  

    result <- list()
    if (!is.null(object$data)) {
      class(object) <- 's2dv_cube'
      result$exp <- object
    }
    if (!is.null(obs$data)) {
      class(obs) <- 's2dv_cube'
      result$obs <- obs
    }
    if (is.list(result)) {
        if (is.null(result$exp)) {
            result <- result$obs
        } else if (is.null(result$obs)) {
            result <- result$exp
        } else {
            warning("The output is a list of two 's2dv_cube' objects",
                    " corresponding to 'exp' and 'obs'.")
        }
    }
 
  } else if (inherits(object, 'startR_array')) {
    result <- list()
    result$data <- as.vector(object)
    dim(result$data) <- dim(object)

    dat_attr_names <- names(attributes(object)$Variables$dat1)
    common_attr_names <- names(attributes(object)$Variables$common)
    # $lon
    known_lon_names <- utils::getFromNamespace(".KnownLonNames", "s2dv")()
    if (!is.null(dat_attr_names[which(dat_attr_names %in% known_lon_names)]) &
        !identical(dat_attr_names[which(dat_attr_names %in% known_lon_names)], character(0))) {
      result$lon <- attributes(object)$Variables$dat1[[dat_attr_names[which(dat_attr_names %in% known_lon_names)]]]
    } else if (!is.null(common_attr_names[which(common_attr_names %in% known_lon_names)]) &
               !identical(common_attr_names[which(common_attr_names %in% known_lon_names)], character(0))) {
      result$lon <- attributes(object)$Variables$common[[common_attr_names[which(common_attr_names %in% known_lon_names)]]]
    } else {
      warning("'lon' is not found in this object.")
      result$lon <- NULL
    }
    # $lat
    known_lat_names <- utils::getFromNamespace(".KnownLatNames", "s2dv")()
    if (!is.null(dat_attr_names[which(dat_attr_names %in% known_lat_names)]) &
        !identical(dat_attr_names[which(dat_attr_names %in% known_lat_names)], character(0))) {
      result$lat <- attributes(object)$Variables$dat1[[dat_attr_names[which(dat_attr_names %in% known_lat_names)]]]
    } else if (!is.null(common_attr_names[which(common_attr_names %in% known_lat_names)]) &
               !identical(common_attr_names[which(common_attr_names %in% known_lat_names)], character(0))) {
      result$lat <- attributes(object)$Variables$common[[common_attr_names[which(common_attr_names %in% known_lat_names)]]]
    } else {
      warning("'lat' is not found in this object.")
      result$lat <- NULL
    }

    vars <- which(!common_attr_names %in% c("time", known_lon_names, known_lat_names))

    if (length(vars) > 1) {
      warning("More than one variable has been provided and ",
              "only the first one '", common_attr_names[vars[1]],"' will be used.")
      vars <- vars[1]
    }

    Variable <- list()
    Variable$varName <- names(attributes(object)$Variables$common)[vars]
    attr(Variable, 'variable') <- attributes(object)$Variables$common[[vars]]
    result$Variable <- Variable
    dims <- dim(object)
    if (any(c('sdate', 'sdates') %in% names(dims))) {
        n_sdates <- dims[which(names(dims) == 'sdate' | names(dims) == 'sdates')]
        sdates <- attributes(object)$Variables$common$time[1 : n_sdates]
    } else {
        sdates <- attributes(object)$Variables$common$time[1]
    }
    Dataset <- list(list(InitializationDates = list(Member_1 = sdates)))
    names(Dataset) <- list(deparse(substitute(object)))
    result$Datasets <- Dataset
    result$Dates$start <- attributes(object)$Variables$common$time
    result$when <- Sys.time()
    result$source_files <- as.vector(attributes(object)$Files)
    result$load_parameters <- attributes(object)$FileSelectors
    class(result) <- 's2dv_cube'
  } else {
    stop("The class of parameter 'object' is not implemented",
         " to be converted into 's2dv_cube' class yet.")
  }
  result

}
