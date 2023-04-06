#'Plots the observed weekly means and climatology of a timeseries data
#' 
#'@description This function plots the observed weekly means and climatology of 
#'a timeseries data using ggplot package. It compares the weekly climatology in 
#'a specified period (reference period) to the observed conditions during the 
#'target period analyzed in the case study (included in the reference period).
#' 
#'@param data A multidimensional array with named dimensions with at least sdate 
#'  and time dimensions containing observed daily data. It can also be a 
#'  dataframe with computed percentiles as input for ggplot. The target year 
#'  must be included in the input data.
#'@param first_date The first date of the target period of timeseries. It can be 
#'  of class 'Date', 'POSIXct' or a character string in the format 'yyyy-mm-dd'. 
#'  It must be a date included in the reference period.
#'@param ref_period_ini A numeric value indicating the first year of the 
#'  reference period.
#'@param ref_period_end A numeric value indicating the last year of the 
#'  reference period.
#'@param time_dim A character string indicating the daily time dimension name. 
#'  The default value is 'time'.
#'@param sdate_dim A character string indicating the start year dimension name. 
#'  The default value is 'sdate'.
#'@param title The text for the top title of the plot. 
#'@param palette A palette name from the R Color Brewer’s package. The default 
#'  value is 'Blues'.
#'@param fileout A character string indicating the file name where to save the 
#'  plot. If not specified (default) a graphics device will pop up.
#'@param device A character string indicating the device to use. Can either be 
#'  a device function (e.g. png), or one of "eps", "ps", "tex" (pictex), 
#'  "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only).
#'@param width A numeric value of the plot width in units ("in", "cm", "mm", or 
#'  "px"). It is set to 8 by default.
#'@param height A numeric value of the plot height in units ("in", "cm", "mm", 
#'  or "px"). It is set to 6 by default.
#'@param units Units of the size of the device (file or window) to plot in. 
#'  Inches (’in’) by default.
#'@param dpi A numeric value of the plot resolution. It is set to 300 by 
#'  default.
#' 
#'@return A ggplot object containing the plot.
#' 
#'@examples
#'data <- array(rnorm(49*20, 274, 7), dim = c(time = 49, sdate = 20))
#'PlotWeeklyClim(data = data, first_date = '2010-08-09', 
#'               ref_period_ini = 1998,
#'               ref_period_end = 2020)
#' 
#'@import multiApply
#'@import lubridate
#'@import ggplot2
#'@import RColorBrewer
#'@import scales
#'@importFrom ClimProjDiags Subset
#'@importFrom s2dv MeanDims
#'@export
PlotWeeklyClim <- function(data, first_date, ref_period_ini, ref_period_end, 
                           time_dim = 'time', sdate_dim = 'sdate',
                           title = "Observed weekly means and climatology",
                           palette = "Blues", fileout = NULL, 
                           device = NULL, width = 8, height = 6,
                           units = 'in', dpi = 300) {
  # Check input arguments
  # data
  if (is.array(data)) {
    if (is.null(names(dim(data)))) {
      stop("Parameter 'data' must have named dimensions.")
    }
    is_array <- TRUE
  } else if (is.data.frame(data)) {
    col_names <- c("week", "clim", "p10", "p90", "p33", "p66", 
                   "week_mean", "day", "data")
    if (!all(col_names %in% names(data))) {
      stop(paste0("If parameter 'data' is a data frame, it must contain the ",
                  "following column names: 'week', 'clim', 'p10', 'p90', 'p33', ", 
                  "'p66', 'week_mean', 'day' and 'data'."))
    }
    is_array <- FALSE
  } else {
    stop("Parameter 'data' must be an array or a data frame.")
  }
  if (is_array) {
    # time_dim
    if (!is.character(time_dim)) {
      stop("Parameter 'time_dim' must be a character string.")
    }
    if (!all(time_dim %in% names(dim(data)))) {
      stop("Parameter 'time_dim' is not found in 'data' dimension.")
    }
    if (dim(data)[time_dim] < 7) {
      stop(paste0("Parameter 'data' must have the dimension 'time_dim' of ",
                  "length equal or grater than 7 to compute the weekly means."))
    }
    # sdate_dim
    if (!is.character(sdate_dim)) {
      stop("Parameter 'sdate_dim' must be a character string.")
    }
    if (!sdate_dim %in% names(dim(data))) {
      warning(paste0("Parameter 'sdate_dim' is not found in 'data' dimension. ",
                     "A dimension of length 1 has been added."))
      data <- InsertDim(data, 1, lendim = 1, name = sdate_dim)
    }
    # ref_period_ini and ref_period_end
    if (!is.numeric(ref_period_ini) | !is.numeric(ref_period_end)) {
      stop("Parameters 'ref_period_ini' and 'ref_period_end' must be numeric.")
    }
    # first_date
    if ((!inherits(first_date, "POSIXct") & !inherits(first_date, "Date")) &&
        (!is.character(first_date) | nchar(first_date) != 10)) {
      stop(paste0("Parameter 'first_date' must be a character string ", 
                  "indicating the date in the format 'yyyy-mm-dd', 'POSIXct' ", 
                  "or 'Dates' class."))
    }
    first_date <- ymd(first_date)
    target_year <- year(first_date)
    if (target_year < ref_period_ini | target_year > ref_period_end) {
      stop("Parameter 'first_date' must be a date included in the reference period.")
    }

    # Dates creation
    dates <- seq(first_date, first_date + days(dim(data)[time_dim]-1), by = "1 day")
    index_first_date <- which(dates == first_date)
    index_last_date <- length(dates) - (length(dates) %% 7)
    last_date <- dates[index_last_date]

    # Weekly aggregations
    data_subset <- Subset(data, along = time_dim, 
                          indices = index_first_date:index_last_date)
    weekly_aggre <- SplitDim(data_subset, split_dim = time_dim, 
                             indices = sort(rep(1:(length(index_first_date:index_last_date)/7), 7)),
                             new_dim_name = 'week')
    weekly_means <- MeanDims(weekly_aggre, time_dim)
    weekly_clim <- MeanDims(weekly_means, sdate_dim)

    weekly_p10 <- Apply(weekly_means, target_dims = sdate_dim,
                        fun = function(x) {quantile(x, 0.10)})$output1
    weekly_p90 <- Apply(weekly_means, target_dims = sdate_dim,
                        fun = function(x) {quantile(x, 0.90)})$output1
    weekly_p33 <- Apply(weekly_means, target_dims = sdate_dim,
                        fun = function(x) {quantile(x, 0.33)})$output1
    weekly_p66 <- Apply(weekly_means, target_dims = sdate_dim,
                        fun = function(x) {quantile(x, 0.66)})$output1

    clim <- p10 <- p90 <- p33 <- p66 <- NULL
    weekly_data <- data.frame(clim = as.vector(weekly_clim),
                              p10 = as.vector(weekly_p10),
                              p90 = as.vector(weekly_p90), 
                              p33 = as.vector(weekly_p33),
                              p66 = as.vector(weekly_p66),
                              week = 1:(length(index_first_date:index_last_date)/7))

    daily <- Subset(data_subset, along = sdate_dim, 
                    indices = which(ref_period_ini:ref_period_end == target_year), 
                    drop = TRUE)

    dims_subset <- names(dim(daily))[which(!names(dim(daily)) %in% c(time_dim, sdate_dim))]

    if (!identical(dims_subset, character(0))) {
      daily <- Subset(daily, dims_subset, as.list(rep(1, length(dims_subset))), drop = TRUE)
    }
    
    daily_data <- data.frame(day = seq(first_date, last_date, by = "1 day"),
                             data = daily, 
                             week = sort(rep(1:(length(index_first_date:index_last_date)/7), 7)))
    week_mean <- aggregate(data ~ week, data = daily_data, mean)
    weekly_data <- cbind(weekly_data, week_mean$data)
    colnames(weekly_data)[7] <- 'week_mean'
    all <- merge(weekly_data, daily_data, by = 'week')
  } else {
    all <- data
  }

  # Create a ggplot object
  cols <- colorRampPalette(brewer.pal(9, palette))(6)

  p <- ggplot(all, aes(x = day)) +
       geom_ribbon(aes(ymin = p10, ymax = p90, group = week, fill = "p10-p90"), 
                   alpha = 0.7) + # extremes clim
       geom_ribbon(aes(ymin = p33, ymax = p66, group = week, fill = "p33-p66"),
                   alpha = 0.7) +  # terciles clim
       geom_line(aes(y = clim, group = week, color = "climatological mean",
                 linetype = "climatological mean"),
                 alpha = 1.0, size = 0.7) + # mean clim
       geom_line(aes(y = data, color = "observed daily mean",
                 linetype = "observed daily mean"),
                 alpha = 1, size = 0.2) + # daily evolution
       geom_line(aes(y = week_mean, group = week, color = "observed weekly mean",
                 linetype = "observed weekly mean"),
                 alpha = 1, size = 0.7) + 		 # weekly evolution
       theme_bw() + ylab(paste0('tas', " (", "deg.C", ")")) + xlab(NULL) + 
       ggtitle(title) +
       scale_x_date(breaks = seq(min(all$day), max(all$day), by = "7 days"),
                    minor_breaks = NULL, expand = c(0.03, 0.03),
                    labels = date_format("%d %b %Y")) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1),
             panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                             colour = "gray92"),
             panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                             colour = "gray92"),
             legend.spacing = unit(-0.2, "cm")) +
       scale_fill_manual(name = NULL,
                         values = c("p10-p90" = cols[3], "p33-p66" = cols[4])) +
       scale_color_manual(name = NULL, values = c("climatological mean" = cols[5],
                                                  "observed daily mean" = "grey20",
                                                  "observed weekly mean" = "black")) +
       scale_linetype_manual(name = NULL, values = c("climatological mean" = "solid",
                                                     "observed daily mean" = "dashed",
                                                     "observed weekly mean" = "solid"),
       guide = guide_legend(override.aes = list(lwd = c(0.7, 0.2, 0.7)))) +
       guides(fill = guide_legend(order = 1))

  # Return the ggplot object
  if (is.null(fileout)) {
    return(p)
  } else {
    ggsave(filename = fileout, plot = p, device = device, height = height, 
           width = width, units = units, dpi = dpi) 
  }
}












  