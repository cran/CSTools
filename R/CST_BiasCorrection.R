#' Bias Correction based on the mean and standard deviation adjustment
#'
#'@author Verónica Torralba, \email{veronica.torralba@bsc.es}
#'@description This function applies the simple bias adjustment technique described in Torralba et al. (2017). The adjusted forecasts have an equivalent standard deviation and mean to that of the reference dataset.
#'
#'@param exp an object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the seasonal forecast experiment data in the element named \code{$data}
#'@param obs an object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the observed data in the element named \code{$data}.
#'@param na.rm a logical value indicating whether missing values should be stripped before the computation proceeds, by default it is set to FALSE.
#'
#'@return an object of class \code{s2dv_cube} containing the bias corrected forecasts in the element called \code{$data} with the same dimensions of the experimental data.
#'
#'@references Torralba, V., F.J. Doblas-Reyes, D. MacLeod, I. Christel and M. Davis (2017). Seasonal climate prediction: a new source of information for the management of wind energy resources. Journal of Applied Meteorology and Climatology, 56, 1231-1247, doi:10.1175/JAMC-D-16-0204.1. (CLIM4ENERGY, EUPORIAS, NEWA, RESILIENCE, SPECS)
#'
#'@import s2dverification
#'@import multiApply
#'@examples
#'
#'# Example
#'# Creation of sample s2dverification objects. These are not complete
#'# s2dverification objects though. The Load function returns complete objects.
#'mod1 <- 1 : (1 * 3 * 4 * 5 * 6 * 7)
#'dim(mod1) <- c(dataset = 1, member = 3, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'obs1 <- 1 : (1 * 1 * 4 * 5 * 6 * 7)
#'dim(obs1) <- c(dataset = 1, member = 1, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'lon <- seq(0, 30, 5)
#'lat <- seq(0, 25, 5)
#'exp <- list(data = mod1, lat = lat, lon = lon)
#'obs <- list(data = obs1, lat = lat, lon = lon)
#'attr(exp, 'class') <- 's2dv_cube'
#'attr(obs, 'class') <- 's2dv_cube'
#'a <- CST_BiasCorrection(exp = exp, obs = obs)
#'str(a)
#'@export
CST_BiasCorrection <- function(exp, obs, na.rm = FALSE) {
  if (!inherits(exp, 's2dv_cube') || !inherits(obs, 's2dv_cube')) {
    stop("Parameter 'exp' and 'obs' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load.")
  }

  if (dim(obs$data)['member'] != 1) {
    stop("The length of the dimension 'member' in the component 'data' ",
         "of the parameter 'obs' must be equal to 1.")
  }
  dimnames <- names(dim(exp$data))
  BiasCorrected <- BiasCorrection(exp = exp$data, obs = obs$data, na.rm = na.rm)
  pos <- match(dimnames, names(dim(BiasCorrected)))
  BiasCorrected <- aperm(BiasCorrected, pos)
  names(dim(BiasCorrected)) <- dimnames
  exp$data <- BiasCorrected
  exp$Datasets <- c(exp$Datasets, obs$Datasets)
  exp$source_files <- c(exp$source_files, obs$source_files)
  return(exp)
}

BiasCorrection <- function (exp, obs , na.rm = FALSE) {
 
   if (!all(c('member', 'sdate') %in% names(dim(exp)))) {
    stop("Parameter 'exp' must have the dimensions 'member' and 'sdate'.")
  }

  if (!all(c('sdate') %in% names(dim(obs)))) {
    stop("Parameter 'obs' must have the dimension 'sdate'.")
  }

  if (any(is.na(exp)))  {
    warning("Parameter 'exp' contains NA values.")
  }
  
  if (any(is.na(obs)))  {
    warning("Parameter 'obs' contains NA values.")
  }

  if (!is.logical(na.rm))  {
    na.rm <- FALSE
    warning("Paramater 'na.rm' must be a logical, it has been set to FALSE.")
  }

  if (length(na.rm)>1) {
    na.rm <- na.rm[1]
    warning("Paramter 'na.rm' has length greater than 1, and only the fist element is used.")
  }

  target_dims_obs <- 'sdate'
  if ('member' %in% names(dim(obs))) {
    target_dims_obs <- c('member', target_dims_obs)
  }

  BiasCorrected <- Apply(data = list(var_obs = obs, var_exp = exp),
                      target_dims = list(target_dims_obs, c('member', 'sdate')),
                      fun = .sbc , na.rm = na.rm)$output1
  return(BiasCorrected)
}

.sbc <- function(var_obs, var_exp , na.rm = FALSE) {
  nmembers <- dim(var_exp)['member'][]
  ntime <- dim(var_exp)['sdate'][]
  if (all(names(dim(var_exp)) != c('member','sdate'))) {
    var_exp <- t(var_exp)
  }
  
  corrected <- NA * var_exp
  
  for (t in 1 : ntime) {
    # defining forecast,hindcast and observation in cross-validation
    fcst <- var_exp[ , t]
    hcst <- var_exp[ , -t]
    obs <- var_obs[-t]
    
    # parameters
    sd_obs <- sd(obs , na.rm = na.rm)
    sd_exp <- sd(hcst , na.rm = na.rm)
    clim_exp <- mean(hcst , na.rm = na.rm)
    clim_obs <- mean(obs , na.rm = na.rm)
    
    # bias corrected forecast
    corrected[ , t] <- ((fcst - clim_exp) * (sd_obs / sd_exp)) + clim_obs
  }
  names(dim(corrected)) <- c('member', 'sdate')
  return(corrected)
}
