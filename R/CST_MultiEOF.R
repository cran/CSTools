#' @rdname CST_MultiEOF
#' @title EOF analysis of multiple variables
#'
#' @author Jost von Hardenberg - ISAC-CNR, \email{j.vonhardenberg@isac.cnr.it}
#' @author Paolo Davini - ISAC-CNR, \email{p.davini@isac.cnr.it}
#'
#' @description This function performs EOF analysis over multiple variables,
#' accepting in input a list of CSTools objects. Based on Singular Value Decomposition. For each field the EOFs are computed and the corresponding PCs are standardized (unit variance, zero mean); the minimum number of principal components needed to reach the user-defined variance is retained. The function weights the input data for the latitude cosine square root.

#'
#' @param datalist A list of objects of the class 's2dv_cube', containing the variables to be analysed.
#' Each data object in the list is expected to have an element named \code{$data} with at least two
#' spatial dimensions named "lon" and "lat", a dimension "ftime" and a dimension "sdate".
#' @param neof_composed Number of composed eofs to return in output
#' @param minvar Minimum variance fraction to be explained in first decomposition
#' @param neof_max Maximum number of single eofs considered in the first decomposition
#' @param lon_lim Vector with longitudinal range limits for the EOF calculation for all input variables
#' @param lat_lim Vector with latitudinal range limits for the EOF calculation for all input variables
#' @return A list with elements \code{$coeff} (an array of time-varying principal component coefficients),
#'         \code{$variance} (a matrix of explained variances),
#'         \code{eof_pattern} (a matrix of EOF patterns obtained by regression for each variable).
#' @import abind
#' @examples
#' \donttest{
#' library(zeallot)
#' library(ClimProjDiags)
#' c(exp, obs) %<-% lonlat_data
#' # Create three datasets (from the members)
#' exp1 <- exp
#' exp2 <- exp
#' exp3 <- exp
#' exp1$data <- Subset(exp$data, along = 2, indices = 1 : 5)
#' exp2$data <- Subset(exp$data, along = 2, indices = 6 : 10)
#' exp3$data <- Subset(exp$data, along = 2, indices = 11 : 15)
#'
#' cal <- CST_MultiEOF(list(exp1, exp2, exp3), neof_max=5, neof_composed=2)
#' str(cal)
#' # List of 3
#' # $ coeff      : num [1:3, 1:6, 1:2, 1:5] -0.312 -0.588 0.724 1.202 1.181 ...
#' # $ variance   : num [1:2, 1:5] 0.413 0.239 0.352 0.27 0.389 ...
#' # $ eof_pattern: num [1:3, 1:53, 1:22, 1:2, 1:5] -1.47 -0.446 -0.656 -1.534 -0.464 ...
#' dim(cal$coeff)
#' #  ftime   sdate     eof  member
#' #      3       6       2       3
#'
#' cal <- CST_MultiEOF(list(exp1, exp2, exp3) , minvar=0.9)
#' str(cal)
#' # $ coeff      : num [1:3, 1:6, 1:5, 1:5] 0.338 0.603 -0.736 -1.191 -1.198 ...
#' # $ variance   : num [1:5, 1:5] 0.3903 0.2264 0.1861 0.1032 0.0379 ...
#' # $ eof_pattern: num [1:3, 1:53, 1:22, 1:5, 1:5] 1.477 0.454 0.651 1.541 0.47 ...
#'
#' cal <- CST_MultiEOF(list(exp1, exp2))
#' cal <- CST_MultiEOF(list(exp1, exp2, exp3), lon_lim=c(5, 30), lat_lim=c(35, 50), neof_composed=3)
#' }
#' @export

CST_MultiEOF <- function(datalist,
                         neof_max = 40, neof_composed = 5, minvar = 0.6,
                         lon_lim = NULL, lat_lim = NULL) {

  if (!(all(sapply(datalist, inherits, 's2dv_cube')))) {
     stop("Elements of the list in parameter 'datalist' must be of the class ",
          "'s2dv_cube', as output by CSTools::CST_Load.")
  }

  # Check if all dims equal
  adims=lapply(lapply(datalist, function(x) x$data), dim)
  if( !all(apply(apply(abind(adims, along = 0), 2, duplicated), 2, sum) ==
          (length(adims)-1))) {
    stop("Input data fields must all have the same dimensions.")
  }

  #print("Pasting data...")
  exp <- abind(lapply(datalist, '[[', 'data'), along = 0)
  dim(exp) <- c(var = length(datalist), dim(datalist[[1]]$data))
  #print("...done")

  if (any(is.na(exp)))  {
    stop("Input data contain NA values.")
  }

  result <- MultiEOF(exp, datalist[[1]]$lon, datalist[[1]]$lat,
                      datalist[[1]]$Dates$start, minvar = minvar,
                      neof_max = neof_max, neof_composed = neof_composed,
                      lon_lim = lon_lim, lat_lim = lat_lim)

  return(result)
}

#' @rdname MultiEOF
#' @title EOF analysis of multiple variables starting from an array (reduced version)
#'
#' @author Jost von Hardenberg - ISAC-CNR, \email{j.vonhardenberg@isac.cnr.it}
#' @author Paolo Davini - ISAC-CNR, \email{p.davini@isac.cnr.it}
#'
#' @description This function performs EOF analysis over multiple variables, accepting in input an array with a dimension \code{"var"} for each variable to analyse. Based on Singular Value Decomposition. For each field the EOFs are computed and the corresponding PCs are standardized (unit variance, zero mean); the minimum number of principal components needed to reach the user-defined variance is retained. The function weights the input data for the latitude cosine square root.
#'
#' @param data A multidimensional array with dimension \code{"var"},
#' containing the variables to be analysed. The other diemnsions follow the same structure as the
#' \code{"exp"} element of a 's2dv_cube' object.
#' @param lon Vector of longitudes.
#' @param lat Vector of latitudes.
#' @param time Vector or matrix of dates in POSIXct format.
#' @param lon_dim String with dimension name of longitudinal coordinate
#' @param lat_dim String with dimension name of latitudinal coordinate
#' @param neof_composed Number of composed eofs to return in output
#' @param minvar Minimum variance fraction to be explained in first decomposition
#' @param neof_max Maximum number of single eofs considered in the first decomposition
#' @param lon_lim Vector with longitudinal range limits for the calculation for all input variables
#' @param lat_lim Vector with latitudinal range limits for the calculation for all input variables
#' @return A list with elements \code{$coeff} (an array of time-varying principal component coefficients),
#'         \code{$variance} (a matrix of explained variances),
#'         \code{eof_pattern} (a matrix of EOF patterns obtained by regression for each variable).
#' @import multiApply
#' @export

MultiEOF <- function(data, lon, lat, time,
                       lon_dim = "lon", lat_dim = "lat",
                       neof_max = 40, neof_composed = 5, minvar = 0.6,
                       lon_lim = NULL, lat_lim = NULL) {
  # Check/detect time_dim

  # reorder and group ftime and sdate together at the end in that order
  cdim0 <- dim(data)
  imaskt <- names(cdim0) %in% "ftime"
  imasks <- names(cdim0) %in% "sdate"
  data <- .aperm2(data, c(which(!(imasks | imaskt)),
                          which(imaskt), which(imasks)))
  cdim <- dim(data)
  ind <- 1:length(which(!(imaskt | imasks)))
  # compact (multiply) time_dim dimensions
  dim(data) <- c(cdim[ind], samples = prod(cdim[-ind]))

  # Repeatedly apply .multi.eofs
  result <- Apply(data, c("var", lon_dim, lat_dim, "samples"),
                  .multi.eofs, lon, lat, time, neof_max = neof_max,
                  neof_composed = neof_composed, minvar = minvar,
                  xlim = lon_lim, ylim = lat_lim)

  # Expand back samples to compacted dims
  dim(result$coeff) <- c(cdim[-ind], dim(result$coeff)[-1])
  # Recover first lon and first lat list
  dd=dim(result$lon)[1]; m=matrix(1, nrow=dd, ncol=length(dim(result$lon)));
  m[1:dd]=1:dd; result$lon = result$lon[m]
  dd=dim(result$lat)[1]; m=matrix(1, nrow=dd, ncol=length(dim(result$lat)));
  m[1:dd]=1:dd; result$lat = result$lat[m]

  return(result)
}


#' Atomic MultiEOF
#' @param field_arr_raw an array of dimension: (n_field, lon, lat, time).
#' where n_field is the number of variables over which to calculate
#' the multi_eofs.
#' @param neof_composed Number of composed eofs to return in output
#' @param minvar Minimum variance fraction to be explained in first decomposition
#' @param neof_max Maximum number of single eofs considered in the first decomposition
#' @param xlim Vector with longitudinal range limits for the calculation
#' @param ylim Vector with latitudinal range limits for the calculation
#' @return A list with elements \code{$coeff} (an array of time-varying principal component coefficients),
#'         \code{$variance} (a matrix of explained variances),
#'         \code{eof_pattern} (a matrix of EOF patterns obtained by regression for each variable).
#' @noRd

.multi.eofs <- function(field_arr_raw, lon, lat, time,
                        neof_max = 40, neof_composed = 5, minvar = 0.6,
                        xlim = NULL, ylim = NULL) {

  if (exists(".lm.fit")) {
    lin.fit <- .lm.fit
  } else {
    lin.fit <- lm.fit
  }

  n_field <- dim(field_arr_raw)[1]

  etime <- .power.date(time)
  #print("Calculating anomalies...")
  field_arr <- array(dim = dim(field_arr_raw))
  for (k in seq(1, n_field, 1)) {
    field_arr[k, , , ] <- .daily.anom.mean(
                            lon, lat, field_arr_raw[k, , , ], etime
                          )
  }

  # area weighting, based on the root of cosine
  #print("Area Weighting...")
  ww <- .area.weight(lon, lat, root = T)

  for (k in seq(1, n_field, 1)) {
    field_orig <- field_arr[k, , , ]
    # calculate the area weight
    field <- sweep(field_orig, c(1, 2), ww, "*")

    idx <- .selbox(lon, lat, xlim, ylim)
    slon <- lon[idx$ilon]
    slat <- lat[idx$ilat]
    field <- field[idx$ilon, idx$ilat, ]

    # transform 3D field in a matrix
    field <- array(field, dim = c(dim(field)[1] * dim(field)[2], dim(field)[3]))

    # calling SVD
    SVD <- svd(field, nu = neof_max, nv = neof_max)

    # extracting EOFs (loading pattern), expansions coefficient
    # and variance explained
    pattern <- array(SVD$u, dim = c(dim(field)[1], dim(field)[2], neof_max))
    coefficient <- SVD$v
    variance <- (SVD$d[1:neof_max]) ^ 2 / sum((SVD$d) ^ 2)
    #print("Accumulated variance:")
    #print(cumsum(variance))
    reqPC <- which(cumsum(variance) > minvar)[1]
    #print("Number of EOFs needed for var:")
    variance <- variance[1:reqPC]
    coefficient <- coefficient[, 1:reqPC]
    if (reqPC == 1) {
      coefficient <- replicate(1, coefficient)
    }
    coefficient <- apply(coefficient, c(2), .standardize)
    regression <- array(NA, dim = c(length(lon), length(lat), neof_max))
    for (i in 1:reqPC) {
      regression[, , i] <- apply(field_orig, c(1, 2),
                                   function(x) lin.fit(as.matrix(coefficient[, i],
                                                       ncol = 1), x)$coefficients
                                )
    }

    assign(
      paste0("pc", k), list(coeff = coefficient, variance = variance,
                            wcoeff = sweep(coefficient, c(2), variance, "*"),
                            regression = regression)
    )
  }

  newpc <- NULL
  for (k in seq(1, n_field, 1)) {
    newpc <- cbind(newpc, get(paste0("pc", k))$wcoeff)
  }
  newpc <- t(newpc)

  #print("Calculating composed EOFs")
  SVD <- svd(newpc, nu = neof_composed, nv = neof_composed)
  # extracting EOFs, expansions coefficient and variance explained
  coefficient <- SVD$v
  variance <- (SVD$d[1:(neof_composed)]) ^ 2 / sum( (SVD$d) ^ 2)
  coefficient <- apply(coefficient, c(2), .standardize)

  # linear regressions on anomalies
  regression <- array(dim = c(n_field, length(lon),
                              length(lat), neof_composed))
  for (k in seq(1, n_field, 1)) {
    #print("Linear Regressions (it can take a while)... ")
    for (i in 1:neof_composed) {
      regression[k, , , i] <- apply(
                                field_arr[k, , , ], c(1, 2),
                                function(x) lin.fit(
                                              as.matrix(coefficient[, i],
                                                        ncol = 1),
                                            x)$coefficients
                              )
    }
  }

  #print("Finalize...")
  names(dim(coefficient)) <- c("time", "eof")
  variance <- array(variance)
  names(dim(variance)) <- "eof"
  names(dim(regression)) <- c("var", "lon", "lat", "eof")

  out <- list(coeff = coefficient, variance = variance, eof_pattern = regression, lon = slon, lat = slat)

  return(out)
}


# new function to create simple list with date values - Oct-18
# it needs a date or PCICt object, and returns also the season subdivision
.power.date <- function(datas, verbose = FALSE) {

  # create a "season" for continuous time, used by persistance tracking
  startpoints <- c(0, which(diff(datas) > 1))
  deltapoints <- diff(c(startpoints, length(datas)))
  seas <- inverse.rle(list(lengths = deltapoints,
                           values = seq(1, length(startpoints))))

  etime <- list(
    day = as.numeric(format(datas, "%d")),
    month = as.numeric(format(datas, "%m")),
    year = as.numeric(format(datas, "%Y")), data = datas, season = seas
  )

  .printv("Time Array Built", verbose)
  .printv(paste("Length:", length(seas)), verbose)
  .printv(paste("From", datas[1], "to", datas[length(seas)]), verbose)
  return(etime)
}


# function for daily anomalies, use array predeclaration and rowMeans (40 times faster!)
.daily.anom.mean <- function(ics, ipsilon, field, etime) {
  condition <- paste(etime$day, etime$month)
  daily <- array(NA, dim = c(length(ics), length(ipsilon),
                             length(unique(condition))))
  anom <- field * NA

  for (t in unique(condition)) {
    if (sum(t == condition) == 1) {
      print(paste0("Cannot compute a mean with only one value: ",
                   "using climatological mean"))
      anom[, , which(t == condition)] <- rowMeans(field, dims = 2)
    } else {
      daily[, , which(t == unique(condition))] <- rowMeans(field[, , t == condition], dims = 2)
      anom[, , which(t == condition)] <- sweep(field[, , which(t == condition)],
                                               c(1, 2),
                                               daily[, , which(t == unique(condition))], "-")
    }
  }
  return(anom)
}
