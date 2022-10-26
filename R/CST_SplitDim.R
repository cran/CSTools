#'Function to Split Dimension
#'
#'@author Nuria Perez-Zanon, \email{nuria.perez@bsc.es}
#'
#'@description This function split a dimension in two. The user can select the dimension to split and provide indices indicating how to split that dimension or dates and the frequency expected (monthly or by day, month and year). The user can also provide a numeric frequency indicating the length of each division.
#'
#'@param data a 's2dv_cube' object
#'@param split_dim a character string indicating the name of the dimension to split
#'@param indices a vector of numeric indices or dates. If left at NULL, the dates provided in the s2dv_cube object (element Dates) will be used.
#'@param freq a character string indicating the frequency: by 'day', 'month' and 'year' or 'monthly' (by default). 'month' identifies months between 1 and 12 independently of the year they belong to, while 'monthly' differenciates months from different years.
#'@param new_dim_name a character string indicating the name of the new dimension.
#'@param insert_ftime an integer indicating the number of time steps to add at the begining of the time series.
#'
#'@details Parameter 'insert_ftime' has been included for the case of using daily data, requiring split the temporal dimensions by months (or similar) and the first lead time doesn't correspondt to the 1st day of the month. In this case, the insert_ftime could be used, to get a final output correctly organized. E.g.: leadtime 1 is the 2nd of November and the input time series extend to the 31st of December. When requiring split by month with \code{inset_ftime = 1}, the 'monthly' dimension of length two will indicate the month (position 1 for November and position 2 for December), dimension 'time' will be length 31. For November, the position 1 and 31 will be NAs, while from positon 2 to 30 will be filled with the data provided. This allows to select correctly days trhough time dimension.
#'@import abind
#'@importFrom ClimProjDiags Subset
#'@examples
#'
#'data <- 1 : 20
#'dim(data) <- c(time = 10, lat = 2)
#'data <-list(data = data)
#'class(data) <- 's2dv_cube'
#'indices <- c(rep(1,5), rep(2,5))
#'new_data <- CST_SplitDim(data, indices = indices)
#'time <- c(seq(ISOdate(1903, 1, 1), ISOdate(1903, 1, 4), "days"),
#'          seq(ISOdate(1903, 2, 1), ISOdate(1903, 2, 4), "days"),
#'          seq(ISOdate(1904, 1, 1), ISOdate(1904, 1, 2), "days"))
#'data <- list(data = data$data, Dates = time)
#'class(data) <- 's2dv_cube'
#'new_data <- CST_SplitDim(data, indices = time)
#'dim(new_data$data)
#'new_data <- CST_SplitDim(data, indices = time, freq = 'day')
#'dim(new_data$data)
#'new_data <- CST_SplitDim(data, indices = time, freq = 'month')
#'dim(new_data$data)
#'new_data <- CST_SplitDim(data, indices = time, freq = 'year')
#'dim(new_data$data)
#'@export
CST_SplitDim <- function(data, split_dim = 'time', indices = NULL,
                         freq = 'monthly', new_dim_name = NULL, insert_ftime = NULL) {
    if (!inherits(data, 's2dv_cube')) {
      stop("Parameter 'data' must be of the class 's2dv_cube', ",
           "as output by CSTools::CST_Load.")
    }
    if (!is.null(insert_ftime)) {
      if (!is.numeric(insert_ftime)) {
         stop("Parameter 'insert_ftime' should be an integer.")
      } else {
        if (length(insert_ftime) > 1) {
          warning("Parameter 'insert_ftime' must be of length 1, and only the",
                  " first element will be used.")
          insert_ftime <- insert_ftime[1]
        }
        # adding NAs at the begining of the data in ftime dim
        ftimedim <- which(names(dim(data$data)) == 'ftime')
        dims <- dim(data$data)
        dims[ftimedim] <- insert_ftime
        empty_array <- array(NA, dims)
        data$data <- abind(empty_array, data$data, along = ftimedim)
        names(dim(data$data)) <- names(dims)
        # adding dates to Dates for the new NAs introduced
        if ((data$Dates[[1]][2] - data$Dates[[1]][1]) == 1) {
            timefreq <- 'days'
        } else {
            timefreq <- 'months'
            warning("Time frequency of forecast time is considered monthly.")
        }
        start <- data$Dates[[1]]
        dim(start) <- c(ftime = length(start)/dims['sdate'], sdate = dims['sdate'])
        #new <- array(NA, prod(dim(data$data)[c('ftime', 'sdate')]))
        # Pending fix transform to UTC when concatenaiting
        data$Dates$start <- do.call(c, lapply(1:dim(start)[2], function(x) {
                               seq(start[1,x] - as.difftime(insert_ftime, 
                               units = timefreq),
                               start[dim(start)[1],x], by = timefreq, tz = "UTC")}))
      }
    } 
    if (is.null(indices)) {
      if (any(split_dim %in% c('ftime', 'time', 'sdate'))) {
        if (is.list(data$Dates)) {
            indices <- data$Dates[[1]]
        } else {
            indices <- data$Dates
        }
        if (any(names(dim(data$data)) %in% 'sdate')) {
            if (!any(names(dim(data$data)) %in% split_dim)) {
                stop("Parameter 'split_dims' must be one of the dimension ",
                     "names in parameter 'data'.")
            }   
            indices <- indices[1 : dim(data$data)[which(names(dim(data$data)) ==
                                                        split_dim)]]     
        }
      }
    }
    data$data <- SplitDim(data$data, split_dim = split_dim, indices = indices,
                        freq = freq, new_dim_name = new_dim_name)
    return(data)
}
#'Function to Split Dimension
#'
#'@author Nuria Perez-Zanon, \email{nuria.perez@bsc.es}
#'
#'@description This function split a dimension in two. The user can select the dimension to split and provide indices indicating how to split that dimension or dates and the frequency expected (monthly or by day, month and year). The user can also provide a numeric frequency indicating the length of each division.
#'
#'@param data an n-dimensional array with named dimensions
#'@param split_dim a character string indicating the name of the dimension to split
#'@param indices a vector of numeric indices or dates
#'@param freq a character string indicating the frequency: by 'day', 'month' and 'year' or 'monthly' (by default). 'month' identifies months between 1 and 12 independetly of the year they belong to, while 'monthly' differenciates months from different years. Parameter 'freq' can also be numeric indicating the length in which to subset the dimension.
#'@param new_dim_name a character string indicating the name of the new dimension.
#'@import abind
#'@importFrom ClimProjDiags Subset
#'@examples
#'
#'data <- 1 : 20
#'dim(data) <- c(time = 10, lat = 2)
#'indices <- c(rep(1,5), rep(2,5))
#'new_data <- SplitDim(data, indices = indices)
#'time <- c(seq(ISOdate(1903, 1, 1), ISOdate(1903, 1, 4), "days"),
#'          seq(ISOdate(1903, 2, 1), ISOdate(1903, 2, 4), "days"),
#'          seq(ISOdate(1904, 1, 1), ISOdate(1904, 1, 2), "days"))
#'new_data <- SplitDim(data, indices = time)
#'new_data <- SplitDim(data, indices = time, freq = 'day')
#'new_data <- SplitDim(data, indices = time, freq = 'month')
#'new_data <- SplitDim(data, indices = time, freq = 'year')
#'@export
SplitDim <- function(data, split_dim = 'time', indices, freq = 'monthly',
                     new_dim_name = NULL) {
    # check data
    if (is.null(data)) {
        stop("Parameter 'data' cannot be NULL.")
    }
    if (is.null(dim(data))) {
        dim(data) = c(time = length(data))
    }
    if (is.null(names(dim(data)))) {
        stop("Parameter 'data' must have dimension names.")
    }
    dims <- dim(data)
    # check split_dim
    if (!is.character(split_dim)) {
        stop("Parameter 'split_dim' must be a character.")
    }
    if (length(split_dim) > 1) {
        split_dim <- split_dim[1]
        warning("Parameter 'split_dim' has length greater than ",
                "one and only the first element will be used.")
    }
    if (!any(names(dims) %in% split_dim)) {
        stop("Parameter 'split_dims' must be one of the dimension ",
             "names in parameter 'data'.")
    }
   pos_split <- which(names(dims) == split_dim)
    # check indices and freq
    if (is.null(indices)) {
        if (!is.numeric(freq)) {
            stop("Parameter 'freq' must be a integer number indicating ",
                 " the length of each chunk.")
        } else {
            if (!((dims[pos_split] / freq) %% 1 == 0)) {
                stop("Parameter 'freq' must be proportional to the ",
                     "length of the 'split_dim' in parameter 'data'.")
            }
            indices <- rep(1 : (dims[pos_split] / freq), freq)
            indices <- sort(indices)
            repited <- sort(unique(indices))
        }
    } else if (is.numeric(indices)) {
        if (!is.null(freq)) {
           if (freq != 'monthly') {
                warning("Parameter 'freq' is not being used since ",
                        "parameter 'indices' is numeric.")
           }
        }
        repited <- sort(unique(indices))
    } else {
        # Indices should be Dates and freq character
        if (!is.character(freq)) {
            stop("Parameter 'freq' must be a character indicating ",
                 "how to divide the dates provided in parameter 'indices'",
                 ", 'monthly', 'anually' or 'daily'.")
         }
         if (!(any(class(indices) %in% c('POSIXct')))) {
             indices <- try( {
                 if (is.character(indices)) {
                     as.POSIXct(indices)
                 } else {
                      as.POSIXct(indices)
                 }
             })
             if ('try-error' %in% class(indices) | 
                 sum(is.na(indices)) == length(indices)) {
                 stop("Dates provided in parameter 'indices' must be of class",
                      " 'POSIXct' or convertable to 'POSIXct'.")
             } 
         }
    }
    if (length(indices) != dims[pos_split]) {
        stop("Parameter 'indices' has different length of parameter ",
             "data in the dimension supplied in 'split_dim'.")
    }
    # check indices as dates:
    if (!is.numeric(indices)) {
        if (freq == 'day') {
            indices <- as.numeric(strftime(indices, format = "%d"))
            repited <- unique(indices)
        } else if (freq == 'month') {
            indices <- as.numeric(strftime(indices, format = "%m"))
            repited <- unique(indices)
        } else if (freq == 'year') {
            indices <- as.numeric(strftime(indices, format = "%Y"))
            repited <- unique(indices)
        } else if (freq == 'monthly' ) {
            indices <- as.numeric(strftime(indices, format = "%m%Y"))
            repited <- unique(indices)
        } else {
            stop("Parameter 'freq' must be numeric or a character: ",
             "by 'day', 'month', 'year' or 'monthly' (for ",
             "distinguishable month).")
       }
    }
    # check new_dim_name
    if (!is.null(new_dim_name)) {
        if (!is.character(new_dim_name)) {
            stop("Parameter 'new_dim_name' must be character string")
        }
        if (length(new_dim_name) > 1) {
            new_dim_name <- new_dim_name[1]
            warning("Parameter 'new_dim_name' has length greater than 1 ",
                    "and only the first elemenst is used.")
        }
    }
    max_times <- max(unlist(lapply(repited, 
                             function(x){sum(indices == x)})))
    data <- lapply(repited, function(x) {rebuild(x, data, along = split_dim,
                        indices = indices, max_times)})
    data <- abind(data, along = length(dims) + 1)
    if (is.character(freq)) {
        names(dim(data)) <- c(names(dims), freq)
    } else {
        names(dim(data)) <- c(names(dims), 'index')
    }
    if (!is.null(new_dim_name)) {
        names(dim(data)) <- c(names(dims), new_dim_name)
    }
return(data)
}

rebuild <- function(x, data, along, indices, max_times) {
    a <- Subset(data, along = along, indices = which(indices == x))
    pos_dim <- which(names(dim(a)) == along)
    if (dim(a)[pos_dim] != max_times) {
        adding <- max_times - dim(a)[pos_dim]
        new_dims <- dim(a)
        new_dims[pos_dim] <- adding
        extra <- array(NA, dim = new_dims)
        a <- abind(a, extra, along = pos_dim)
        names(dim(a)) <- names(dim(data))
    }
    return(a)
}

