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
