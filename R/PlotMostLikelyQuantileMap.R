#'Plot Maps of Most Likely Quantiles
#'
#'@author Veronica Torralba, \email{veronica.torralba@bsc.es}, Nicolau Manubens, \email{nicolau.manubens@bsc.es}
#'@description This function receives as main input (via the parameter \code{probs}) a collection of longitude-latitude maps, each containing the probabilities (from 0 to 1) of the different grid cells of belonging to a category. As many categories as maps provided as inputs are understood to exist. The maps of probabilities must be provided on a common rectangular regular grid, and a vector with the longitudes and a vector with the latitudes of the grid must be provided. The input maps can be provided in two forms, either as a list of multiple two-dimensional arrays (one for each category) or as a three-dimensional array, where one of the dimensions corresponds to the different categories.
#'
#'@param probs a list of bi-dimensional arrays with the named dimensions 'latitude' (or 'lat') and 'longitude' (or 'lon'), with equal size and in the same order, or a single tri-dimensional array with an additional dimension (e.g. 'bin') for the different categories. The arrays must contain probability values between 0 and 1, and the probabilities for all categories of a grid cell should not exceed 1 when added.
#'@param lon a numeric vector with the longitudes of the map grid, in the same order as the values along the corresponding dimension in \code{probs}.
#'@param lat a numeric vector with the latitudes of the map grid, in the same order as the values along the corresponding dimension in \code{probs}.
#'@param cat_dim the name of the dimension along which the different categories are stored in \code{probs}. This only applies if \code{probs} is provided in the form of 3-dimensional array. The default expected name is 'bin'.
#'@param bar_titles vector of character strings with the names to be drawn on top of the color bar for each of the categories. As many titles as categories provided in \code{probs} must be provided.
#'@param col_unknown_cat character string with a colour representation of the colour to be used to paint the cells for which no category can be clearly assigned. Takes the value 'white' by default.
#'@param ... additional parameters to be sent to \code{PlotCombinedMap} and \code{PlotEquiMap}.
#'@seealso \code{PlotCombinedMap} and \code{PlotEquiMap}
#' 
#'@import s2dverification
#'@importFrom maps map 
#'@importFrom graphics box image layout mtext par plot.new
#'@importFrom grDevices adjustcolor bmp colorRampPalette dev.cur dev.new dev.off hcl jpeg pdf png postscript svg tiff
#'@examples
#'# Simple example
#'x <- array(1:(20 * 10), dim = c(lat = 10, lon = 20)) / 200
#'a <- x * 0.6
#'b <- (1 - x) * 0.6
#'c <- 1 - (a + b)
#'lons <- seq(0, 359.5, length = 20)
#'lats <- seq(-89.5, 89.5, length = 10)
#'PlotMostLikelyQuantileMap(list(a, b, c), lons, lats, 
#'                          toptitle = 'Most likely tercile map',
#'                          bar_titles = paste('% of belonging to', c('a', 'b', 'c')), 
#'                          brks = 20, width = 10, height = 8)
#'
#'# More complex example
#'n_lons <- 40
#'n_lats <- 20
#'n_timesteps <- 100
#'n_bins <- 4
#'
#'# 1. Generation of sample data
#'lons <- seq(0, 359.5, length = n_lons)
#'lats <- seq(-89.5, 89.5, length = n_lats)
#'
#'# This function builds a 3-D gaussian at a specified point in the map.
#'make_gaussian <- function(lon, sd_lon, lat, sd_lat) {
#'  w <- outer(lons, lats, function(x, y) dnorm(x, lon, sd_lon) * dnorm(y, lat, sd_lat))
#'  min_w <- min(w)
#'  w <- w - min_w
#'  w <- w / max(w)
#'  w <- t(w)
#'  names(dim(w)) <- c('lat', 'lon')
#'  w
#'}
#'
#'# This function generates random time series (with values ranging 1 to 5) 
#'# according to 2 input weights.
#'gen_data <- function(w1, w2, n) {
#'  r <- sample(1:5, n, 
#'              prob = c(.05, .9 * w1, .05, .05, .9 * w2), 
#'              replace = TRUE)
#'  r <- r + runif(n, -0.5, 0.5)
#'  dim(r) <- c(time = n)
#'  r
#'}
#'
#'# We build two 3-D gaussians.
#'w1 <- make_gaussian(120, 80, 20, 30)
#'w2 <- make_gaussian(260, 60, -10, 40)
#'
#'# We generate sample data (with dimensions time, lat, lon) according 
#'# to the generated gaussians
#'sample_data <- multiApply::Apply(list(w1, w2), NULL, 
#'                                 gen_data, n = n_timesteps)$output1
#'
#'# 2. Binning sample data
#'prob_thresholds <- 1:n_bins / n_bins
#'prob_thresholds <- prob_thresholds[1:(n_bins - 1)]
#'thresholds <- quantile(sample_data, prob_thresholds)
#'
#'binning <- function(x, thresholds) {
#'  n_samples <- length(x)
#'  n_bins <- length(thresholds) + 1
#'
#'  thresholds <- c(thresholds, max(x))
#'  result <- 1:n_bins
#'  lower_threshold <- min(x) - 1
#'  for (i in 1:n_bins) {
#'    result[i] <- sum(x > lower_threshold & x <= thresholds[i]) / n_samples
#'    lower_threshold <- thresholds[i]
#'  }
#'
#'  dim(result) <- c(bin = n_bins)
#'  result
#'}
#'
#'bins <- multiApply::Apply(sample_data, 'time', binning, thresholds)$output1
#'
#'# 3. Plotting most likely quantile/bin
#'PlotMostLikelyQuantileMap(bins, lons, lats, 
#'                          toptitle = 'Most likely quantile map',
#'                          bar_titles = paste('% of belonging to', letters[1:n_bins]),
#'                          mask = 1 - (w1 + w2 / max(c(w1, w2))), 
#'                          brks = 20, width = 10, height = 8)
#'
#'@export
PlotMostLikelyQuantileMap <- function(probs, lon, lat, cat_dim = 'bin',
                                      bar_titles = NULL, 
                                      col_unknown_cat = 'white',
                                      ...) {
  # Check probs
  error <- FALSE
  if (is.list(probs)) {
    if (length(probs) < 1) {
      stop("Parameter 'probs' must be of length >= 1 if provided as a list.")
    }
    check_fun <- function(x) {
      is.numeric(x) && (length(dim(x)) == 2)
    }
    if (!all(sapply(probs, check_fun))) {
      error <- TRUE
    }
    ref_dims <- dim(probs[[1]])
    equal_dims <- all(sapply(probs, function(x) identical(dim(x), ref_dims)))
    if (!equal_dims) {
      stop("All arrays in parameter 'probs' must have the same dimension ",
           "sizes and names when 'probs' is provided as a list of arrays.")
    }
    num_probs <- length(probs)
    probs <- unlist(probs)
    dim(probs) <- c(ref_dims, map = num_probs)
    cat_dim <- 'map'
  }
  if (!is.numeric(probs)) {
    error <- TRUE
  }
  if (is.null(dim(probs))) {
    error <- TRUE
  }
  if (length(dim(probs)) != 3) {
    error <- TRUE
  }
  if (error) {
    stop("Parameter 'probs' must be either a numeric array with 3 dimensions ",
         " or a list of numeric arrays of the same size with the 'lon' and ",
         "'lat' dimensions.")
  }
  dimnames <- names(dim(probs))
  
  # Check cat_dim
  if (is.character(cat_dim)) {
    if (is.null(dimnames)) {
      stop("Specified a dimension name in 'cat_dim' but no dimension names provided ",
           "in 'probs'.")
    }
    cat_dim <- which(dimnames == cat_dim)
    if (length(cat_dim) < 1) {
      stop("Dimension 'cat_dim' not found in 'probs'.")
    }
    cat_dim <- cat_dim[1]
  } else if (!is.numeric(cat_dim)) {
    stop("Parameter 'cat_dim' must be either a numeric value or a ",
         "dimension name.")
  }
  if (length(cat_dim) != 1) {
    stop("Parameter 'cat_dim' must be of length 1.")
  }
  cat_dim <- round(cat_dim)
  nprobs <- dim(probs)[cat_dim]
  
  # Check bar_titles
  if (is.null(bar_titles)) {
    if (nprobs == 3) {
      bar_titles <- list("Below normal (%)", "Normal (%)", "Above normal (%)")
    } else if (nprobs == 5) {
      bar_titles <- list("Low (%)", "Below normal (%)", 
                         "Normal (%)", "Above normal (%)", "High (%)")
    } else {
      bar_titles <- paste0("Cat. ", 1:nprobs, " (%)")
    }
  }
  
  minimum_value <- ceiling(1 / nprobs * 10 * 1.1) * 10

  # By now, the PlotCombinedMap function is included below in this file. 
  # In the future, PlotCombinedMap will be part of s2dverification and will 
  # be properly imported.
  PlotCombinedMap(probs * 100, lon, lat, map_select_fun = max, 
                  display_range = c(minimum_value, 100),
                  map_dim = cat_dim, 
                  bar_titles = bar_titles, 
                  col_unknown_map = col_unknown_cat, ...)
}

PlotCombinedMap <- function(maps, lon, lat, 
                            map_select_fun, display_range, 
                            map_dim = 'map',
                            brks = NULL, cols = NULL,  
                            col_unknown_map = 'white',
                            mask = NULL, col_mask = 'grey',
                            bar_titles = NULL, legend_scale = 1,
                            fileout = NULL, width = 8, height = 5, 
                            size_units = 'in', res = 100, 
                            ...) {
  args <- list(...)

  # If there is any filenames to store the graphics, process them
  # to select the right device 
  if (!is.null(fileout)) {
    deviceInfo <- .SelectDevice(fileout = fileout, width = width, height = height, 
                                units = size_units, res = res)
    saveToFile <- deviceInfo$fun
    fileout <- deviceInfo$files
  }
  
  # Check probs
  error <- FALSE
  if (is.list(maps)) {
    if (length(maps) < 1) {
      stop("Parameter 'maps' must be of length >= 1 if provided as a list.")
    }
    check_fun <- function(x) {
      is.numeric(x) && (length(dim(x)) == 2)
    }
    if (!all(sapply(maps, check_fun))) {
      error <- TRUE
    }
    ref_dims <- dim(maps[[1]])
    equal_dims <- all(sapply(maps, function(x) identical(dim(x), ref_dims)))
    if (!equal_dims) {
      stop("All arrays in parameter 'maps' must have the same dimension ",
           "sizes and names when 'maps' is provided as a list of arrays.")
    }
    num_maps <- length(maps)
    maps <- unlist(maps)
    dim(maps) <- c(ref_dims, map = num_maps)
    map_dim <- 'map'
  }
  if (!is.numeric(maps)) {
    error <- TRUE
  }
  if (is.null(dim(maps))) {
    error <- TRUE
  }
  if (length(dim(maps)) != 3) {
    error <- TRUE
  }
  if (error) {
    stop("Parameter 'maps' must be either a numeric array with 3 dimensions ",
         " or a list of numeric arrays of the same size with the 'lon' and ",
         "'lat' dimensions.")
  }
  dimnames <- names(dim(maps))
  
  # Check map_dim
  if (is.character(map_dim)) {
    if (is.null(dimnames)) {
      stop("Specified a dimension name in 'map_dim' but no dimension names provided ",
           "in 'maps'.")
    }
    map_dim <- which(dimnames == map_dim)
    if (length(map_dim) < 1) {
      stop("Dimension 'map_dim' not found in 'maps'.")
    } else {
      map_dim <- map_dim[1]
    }
  } else if (!is.numeric(map_dim)) {
    stop("Parameter 'map_dim' must be either a numeric value or a ",
         "dimension name.")
  }
  if (length(map_dim) != 1) {
    stop("Parameter 'map_dim' must be of length 1.")
  }
  map_dim <- round(map_dim)
  
  # Work out lon_dim and lat_dim
  lon_dim <- NULL
  if (!is.null(dimnames)) {
    lon_dim <- which(dimnames %in% c('lon', 'longitude'))[1]
  }
  if (length(lon_dim) < 1) {
    lon_dim <- (1:3)[-map_dim][1]
  }
  lon_dim <- round(lon_dim)
  
  lat_dim <- NULL
  if (!is.null(dimnames)) {
    lat_dim <- which(dimnames %in% c('lat', 'latitude'))[1]
  }
  if (length(lat_dim) < 1) {
    lat_dim <- (1:3)[-map_dim][2]
  }
  lat_dim <- round(lat_dim)
  
  # Check lon
  if (!is.numeric(lon)) {
    stop("Parameter 'lon' must be a numeric vector.")
  }
  if (length(lon) != dim(maps)[lon_dim]) {
    stop("Parameter 'lon' does not match the longitude dimension in 'maps'.")
  }
  
  # Check lat
  if (!is.numeric(lat)) {
    stop("Parameter 'lat' must be a numeric vector.")
  }
  if (length(lat) != dim(maps)[lat_dim]) {
    stop("Parameter 'lat' does not match the longitude dimension in 'maps'.")
  }
  
  # Check map_select_fun
  if (is.numeric(map_select_fun)) {
    if (length(dim(map_select_fun)) != 2) {
      stop("Parameter 'map_select_fun' must be an array with dimensions ",
           "'lon' and 'lat' if provided as an array.")
    }
    if (!identical(dim(map_select_fun), dim(maps)[-map_dim])) {
      stop("The dimensions 'lon' and 'lat' in the 'map_select_fun' array must ",
           "have the same size, name and order as in the 'maps' parameter.")
    }
  }
  if (!is.function(map_select_fun)) {
    stop("The parameter 'map_select_fun' must be a function or a numeric array.")
  }
  
  # Check display_range
  if (!is.numeric(display_range) || length(display_range) != 2) {
    stop("Parameter 'display_range' must be a numeric vector of length 2.")
  }
  
  # Check brks
  if (is.null(brks) || (is.numeric(brks) && length(brks) == 1)) {
    num_brks <- 5
    if (is.numeric(brks)) {
      num_brks <- brks
    }
    brks <- seq(from = display_range[1], to = display_range[2], length.out = num_brks)
  }
  if (!is.numeric(brks)) {
    stop("Parameter 'brks' must be a numeric vector.")
  }
  
  # Check cols
  col_sets <- list(c("#A1D99B", "#74C476", "#41AB5D", "#238B45"),
                   c("#6BAED6FF", "#4292C6FF", "#2171B5FF", "#08519CFF"),
                   c("#FFEDA0FF", "#FED976FF", "#FEB24CFF", "#FD8D3CFF"),
                   c("#FC4E2AFF", "#E31A1CFF", "#BD0026FF", "#800026FF"),
                   c("#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497"))
  if (is.null(cols)) {
    if (length(col_sets) >= dim(maps)[map_dim]) {
      chosen_sets <- 1:(dim(maps)[map_dim])
      chosen_sets <- chosen_sets + floor((length(col_sets) - length(chosen_sets)) / 2)
    } else {
      chosen_sets <- array(1:length(col_sets), dim(maps)[map_dim])
    }
    cols <- col_sets[chosen_sets]
  } else {
    if (!is.list(cols)) {
      stop("Parameter 'cols' must be a list of character vectors.")
    }
    if (!all(sapply(cols, is.character))) {
      stop("Parameter 'cols' must be a list of character vectors.")
    }
    if (length(cols) != dim(maps)[map_dim]) {
      stop("Parameter 'cols' must be a list of the same length as the number of ",
           "maps in 'maps'.")
    }
  }
  for (i in 1:length(cols)) {
    if (length(cols[[i]]) != (length(brks) - 1)) {
      cols[[i]] <- colorRampPalette(cols[[i]])(length(brks) - 1)
    }
  }
  
  # Check bar_titles
  if (is.null(bar_titles)) {
    if (!is.null(names(cols))) {
      bar_titles <- names(cols)
    } else {
      bar_titles <- paste0("Map ", 1:length(cols))
    }
  } else {
    if (!is.character(bar_titles)) {
      stop("Parameter 'bar_titles' must be a character vector.")
    }
    if (length(bar_titles) != length(cols)) {
      stop("Parameter 'bar_titles' must be of the same length as the number of ",
           "maps in 'maps'.")
    }
  }
  
  # Check legend_scale
  if (!is.numeric(legend_scale)) {
    stop("Parameter 'legend_scale' must be numeric.")
  }
  
  # Check col_unknown_map
  if (!is.character(col_unknown_map)) {
    stop("Parameter 'col_unknown_map' must be a character string.")
  }
  
  # Check col_mask
  if (!is.character(col_mask)) {
    stop("Parameter 'col_mask' must be a character string.")
  }
  
  # Check mask
  if (!is.null(mask)) {
    if (!is.numeric(mask)) {
      stop("Parameter 'mask' must be numeric.")
    }
    if (length(dim(mask)) != 2) {
      stop("Parameter 'mask' must have two dimensions.")
    }
    if ((dim(mask)[1] != dim(maps)[lat_dim]) || 
        (dim(mask)[2] != dim(maps)[lon_dim])) {
      stop("Parameter 'mask' must have dimensions c(lat, lon).")
    }
  }
  
  #----------------------
  # Identify the most likely map
  #----------------------
  if (is.function(map_select_fun)) {
    range_width <- display_range[2] - display_range[1]
    ml_map <- apply(maps, c(lat_dim, lon_dim), function(x) {
      if (any(is.na(x))) {
        res <- NA
      } else {
        res <- which(x == map_select_fun(x))
        if (length(res) > 0) {
          res <- res[1]
          if (map_select_fun(x) < display_range[1] || map_select_fun(x) > display_range[2]) {
            res <- -0.5
          } else {
            res <- res + (map_select_fun(x) - display_range[1]) / range_width
          }
        } else {
          res <- -0.5
        }
      }
      res
    })
  } else {
    stop("Providing 'map_select_fun' as array not implemented yet.")
    ml_map <- map_select_fun
  }
  nmap <- dim(maps)[map_dim]
  nlat <- length(lat)
  nlon <- length(lon)
  
  #----------------------
  # Set latitudes from minimum to maximum
  #----------------------
  if (lat[1] > lat[nlat]){
    lat <- lat[nlat:1]
    indices <- list(nlat:1, TRUE)
    ml_map <- do.call("[", c(list(x = ml_map), indices))
    if (!is.null(mask)){
      mask <- mask[nlat:1, ]
    }
  }
  
  #----------------------
  # Set layout and parameters
  #----------------------
  # Open connection to graphical device
  if (!is.null(fileout)) {
    saveToFile(fileout)
  } else if (names(dev.cur()) == 'null device') {
    dev.new(units = size_units, res = res, width = width, height = height)
  }
  plot.new()
  par(font.main = 1)
  layout(matrix(c(rep(1, nmap),2:(nmap + 1)), 2, nmap, byrow = TRUE), heights = c(6, 1.5))
  
  #----------------------
  # Set colors and breaks and then PlotEquiMap
  #----------------------
  tcols <- c(col_unknown_map, cols[[1]])
  for (k in 2:nmap) {
    tcols <- append(tcols, c(col_unknown_map, cols[[k]]))
  }
  brks_norm <- seq(0, 1, length.out = length(brks))
  tbrks <- c(-1, brks_norm + rep(1:nmap, each = length(brks)))
  PlotEquiMap(var = ml_map, lon = lon, lat = lat, 
              brks = tbrks, cols = tcols, drawleg = FALSE, 
              filled.continents = FALSE, ...)
  
  #----------------------
  # Add overplot on top
  #----------------------
  if (!is.null(mask)) {
    cols_mask <- sapply(seq(from = 0, to = 1, length.out = 10), 
                        function(x) adjustcolor(col_mask, alpha.f = x))
    image(lon, lat, t(mask), axes = FALSE, col = cols_mask, 
          breaks = seq(from = 0, to = 1, by = 0.1), 
          xlab='', ylab='', add = TRUE, xpd = TRUE)
    if (!exists('coast_color')) {
      coast_color <- 'black'
    }
    if (min(lon) < 0) {
      map('world', interior = FALSE, add = TRUE, lwd = 1, col = coast_color) # Low resolution world map (lon -180 to 180).
    } else {
      map('world2', interior = FALSE, add = TRUE, lwd = 1, col = coast_color) # Low resolution world map (lon 0 to 360).
    }
    box()
  }
  
  #----------------------
  # Add colorbars 
  #----------------------
  if ('toptitle' %in% names(args)) {
    size_title <- 1
    if ('title_scale' %in% names(args)) {
      size_title <- args[['title_scale']]
    }
    old_mar <- par('mar')
    old_mar[3] <- old_mar[3] - (2 * size_title + 1)
    par(mar = old_mar)
  }
  for (k in 1:nmap){
    ColorBar(brks = brks, cols = cols[[k]], vertical = FALSE, 
             draw_separators = TRUE, extra_margin = c(2, 0, 2, 0), 
             label_scale = legend_scale * 1.5)
    if (!is.null(bar_titles)) {
      mtext(bar_titles[[k]], 3, line = -3, cex = 1.5)
    }
  }
  
  # If the graphic was saved to file, close the connection with the device
  if (!is.null(fileout)) dev.off()
}

# Once PlotCombined is included in s2dverification and removed from
# CSTools, this function will be removed from CSTools too.
.SelectDevice <- function(fileout, width, height, units, res) {
  # This function is used in the plot functions to check the extension of the 
  # files where the graphics will be stored and select the right R device to 
  # save them.
  # If the vector of filenames ('fileout') has files with different 
  # extensions, then it will only accept the first one, changing all the rest 
  # of the filenames to use that extension.

  # We extract the extension of the filenames: '.png', '.pdf', ...
  ext <- regmatches(fileout, regexpr("\\.[a-zA-Z0-9]*$", fileout))

  if (length(ext) != 0) {
    # If there is an extension specified, select the correct device
    ## units of width and height set to accept inches
    if (ext[1] == ".png") {
      saveToFile <- function(fileout) {
        png(filename = fileout, width = width, height = height, res = res, units = units)
      }
    } else if (ext[1] == ".jpeg") {
      saveToFile <- function(fileout) {
        jpeg(filename = fileout, width = width, height = height, res = res, units = units)
      }
    } else if (ext[1] %in% c(".eps", ".ps")) {
      saveToFile <- function(fileout) {
        postscript(file = fileout, width = width, height = height)
      }
    } else if (ext[1] == ".pdf") {
      saveToFile <- function(fileout) {
        pdf(file = fileout, width = width, height = height)
      }
    } else if (ext[1] == ".svg") {
      saveToFile <- function(fileout) {
        svg(filename = fileout, width = width, height = height)
      }
    } else if (ext[1] == ".bmp") {
      saveToFile <- function(fileout) {
        bmp(filename = fileout, width = width, height = height, res = res, units = units)
      }
    } else if (ext[1] == ".tiff") {
      saveToFile <- function(fileout) {
        tiff(filename = fileout, width = width, height = height, res = res, units = units)
      }
    } else {
      warning("file extension not supported, it will be used '.eps' by default.")
      ## In case there is only one filename
      fileout[1] <- sub("\\.[a-zA-Z0-9]*$", ".eps", fileout[1])
      ext[1] <- ".eps"
      saveToFile <- function(fileout) {
        postscript(file = fileout, width = width, height = height)
      }
    }
    # Change filenames when necessary
    if (any(ext != ext[1])) {
      warning(paste0("some extensions of the filenames provided in 'fileout' are not ", ext[1],". The extensions are being converted to ", ext[1], "."))
      fileout <- sub("\\.[a-zA-Z0-9]*$", ext[1], fileout)
    }
  } else {
    # Default filenames when there is no specification
    warning("there are no extensions specified in the filenames, default to '.eps'")
    fileout <- paste0(fileout, ".eps")
    saveToFile <- postscript
  }

  # return the correct function with the graphical device, and the correct 
  # filenames
  list(fun = saveToFile, files = fileout)
}
