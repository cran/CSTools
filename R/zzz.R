# Function to permute arrays of non-atomic elements (e.g. POSIXct)
.aperm2 <- function(x, new_order) {
  old_dims <- dim(x)
  attr_bk <- attributes(x)
  if ('dim' %in% names(attr_bk)) {
    attr_bk[['dim']] <- NULL
  }
  if (is.numeric(x)) {
    x <- aperm(x, new_order)
  } else {
    y <- array(1:length(x), dim = dim(x))
    y <- aperm(y, new_order)
    x <- x[as.vector(y)]
  }
  dim(x) <- old_dims[new_order]
  attributes(x) <- c(attributes(x), attr_bk)
  x
}

# verbose-only printing function
.printv <- function(value, verbosity = TRUE) {
  if (verbosity) {
    print(value)
  }
}

# normalize a time series
.standardize <- function(timeseries) {
  out <- (timeseries - mean(timeseries, na.rm = T)) / sd(timeseries, na.rm = T)
  return(out)
}

.selbox <- function(lon, lat, xlim = NULL, ylim = NULL) {
  if (!is.null(xlim)) {
    # This transforms c(-20, -10) to c(340, 350) but c(-20, 10) is unchanged
    # Bring them all to the same units in the 0:360 range
    xlim1 <- xlim[1] %% 360
    xlim2 <- xlim[2] %% 360
    lonm <- lon %% 360
    if (lonm[1] > tail(lonm, 1)) {
      lonm <- lon
    }
    if (xlim1 > xlim2) {
      # If box crosses 0
      ilonsel <- (lonm >= xlim1) | (lonm <= xlim2)
    } else {
      ilonsel <- (lonm >= xlim1) & (lonm <= xlim2)
    }
    if (!any(ilonsel)) {
      stop("No intersection between longitude bounds and data domain.")
    }
  } else {
    ilonsel <- 1:length(lon)
  }
  if (!is.null(ylim)) {
    ilatsel <- (lat >= ylim[1]) & (lat <= ylim[2])
  } else {
    ilatsel <- 1:length(lat)
  }
  return(list(ilon = ilonsel, ilat = ilatsel))
}

# produce a 2d matrix of area weights
.area.weight <- function(ics, ipsilon, root = T) {
  field <- array(NA, dim = c(length(ics), length(ipsilon)))
  if (root == T) {
    for (j in 1:length(ipsilon)) {
      field[, j] <- sqrt(cos(pi / 180 * ipsilon[j]))
    }
  }

  if (root == F) {
    for (j in 1:length(ipsilon)) {
      field[, j] <- cos(pi / 180 * ipsilon[j])
    }
  }

  return(field)
}

#Draws Color Bars for Categories
#A wrapper of s2dv::ColorBar to generate multiple color bars for different 
#categories, and each category has different color set.
GradientCatsColorBar <- function(nmap, brks = NULL, cols = NULL, vertical = TRUE, subsampleg = NULL,
                                 bar_limits, var_limits = NULL,
                                 triangle_ends = NULL, plot = TRUE,
                                 draw_separators = FALSE,
                                 bar_titles = NULL, title_scale = 1, label_scale = 1, extra_margin = rep(0, 4),
                                 ...) {
  # bar_limits
  if (!is.numeric(bar_limits) || length(bar_limits) != 2) {
    stop("Parameter 'bar_limits' must be a numeric vector of length 2.")
  }

  # Check brks
  if (is.null(brks) || (is.numeric(brks) && length(brks) == 1)) {
    num_brks <- 5
    if (is.numeric(brks)) {
      num_brks <- brks
    }
    brks <- seq(from = bar_limits[1], to = bar_limits[2], length.out = num_brks)
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
    if (length(col_sets) >= nmap) {
      chosen_sets <- 1:nmap
      chosen_sets <- chosen_sets + floor((length(col_sets) - length(chosen_sets)) / 2)
    } else {
      chosen_sets <- array(1:length(col_sets), nmap)
    }
    cols <- col_sets[chosen_sets]
  } else {
    if (!is.list(cols)) {
      stop("Parameter 'cols' must be a list of character vectors.")
    }
    if (!all(sapply(cols, is.character))) {
      stop("Parameter 'cols' must be a list of character vectors.")
    }
    if (length(cols) != nmap) {
      stop("Parameter 'cols' must be a list of the same length as the number of ",
           "maps in 'maps'.")
    }
  }
  for (i in 1:length(cols)) {
    if (length(cols[[i]]) != (length(brks) - 1)) {
      cols[[i]] <- grDevices::colorRampPalette(cols[[i]])(length(brks) - 1)
    }
  }

  # Check bar_titles
  if (is.null(bar_titles)) {
    if (nmap == 3) {
      bar_titles <- c("Below normal (%)", "Normal (%)", "Above normal (%)")
    } else if (nmap == 5) {
      bar_titles <- c("Low (%)", "Below normal (%)",
                         "Normal (%)", "Above normal (%)", "High (%)")
    } else {
      bar_titles <- paste0("Cat. ", 1:nmap, " (%)")
    }
  }

  if (plot) {
    for (k in 1:nmap) {
      s2dv::ColorBar(brks = brks, cols = cols[[k]], vertical = FALSE, subsampleg = subsampleg,
#                     bar_limits = bar_limits, var_limits = var_limits,
                     triangle_ends = triangle_ends, plot = TRUE,
                     draw_separators = draw_separators,
                     title = bar_titles[[k]], title_scale = title_scale,
                     label_scale = label_scale, extra_margin = extra_margin)
    }
  } else {
    #TODO: col_inf and col_sup
    return(list(brks = brks, cols = cols))
  }

}

