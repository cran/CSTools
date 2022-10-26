#' Bias Correction based on the mean and standard deviation adjustment
#'
#'@author Verónica Torralba, \email{veronica.torralba@bsc.es}
#'@description This function applies the simple bias adjustment technique 
#'described in Torralba et al. (2017). The adjusted forecasts have an equivalent 
#'standard deviation and mean to that of the reference dataset.
#'
#'@param exp An object of class \code{s2dv_cube} as returned by \code{CST_Load} 
#'  function, containing the seasonal forecast experiment data in the element 
#'  named \code{$data}
#'@param obs An object of class \code{s2dv_cube} as returned by \code{CST_Load} 
#'  function, containing the observed data in the element named \code{$data}.
#'@param exp_cor An object of class \code{s2dv_cube} as returned by 
#'  \code{CST_Load} function, containing the seasonl forecast experiment to be 
#'  corrected. If it is NULL, the 'exp' forecast will be corrected.
#'@param na.rm A logical value indicating whether missing values should be 
#'  stripped before the computation proceeds, by default it is set to FALSE.
#'@param memb_dim A character string indicating the name of the member 
#'  dimension. By default, it is set to 'member'.
#'@param sdate_dim A character string indicating the name of the start date 
#'  dimension. By default, it is set to 'sdate'.
#'@param ncores An integer that indicates the number of cores for parallel 
#'  computations using multiApply function. The default value is NULL.
#'@return An object of class \code{s2dv_cube} containing the bias corrected 
#'forecasts with the same dimensions of the experimental data.
#'
#'@references Torralba, V., F.J. Doblas-Reyes, D. MacLeod, I. Christel and M. 
#'Davis (2017). Seasonal climate prediction: a new source of information for 
#'the management of wind energy resources. Journal of Applied Meteorology and 
#'Climatology, 56, 1231-1247, doi:10.1175/JAMC-D-16-0204.1. (CLIM4ENERGY, 
#'EUPORIAS, NEWA, RESILIENCE, SPECS)
#'
#'@examples
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
#'@import multiApply
#'@export
CST_BiasCorrection <- function(exp, obs, exp_cor = NULL, na.rm = FALSE, 
                               memb_dim = 'member', sdate_dim = 'sdate',
                               ncores = NULL) {
  if (!inherits(exp, 's2dv_cube') || !inherits(obs, 's2dv_cube')) {
    stop("Parameter 'exp' and 'obs' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load.")
  }
  if (!is.null(exp_cor)) {
    if (!inherits(exp_cor, 's2dv_cube')) {
      stop("Parameter 'exp_cor' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load.")
    }
    dimnames <- names(dim(exp_cor$data))
  } else {
    dimnames <- names(dim(exp$data))
  }

  BiasCorrected <- BiasCorrection(exp = exp$data, obs = obs$data, exp_cor = exp_cor$data, 
                                  memb_dim = memb_dim, sdate_dim = sdate_dim,
                                  na.rm = na.rm, ncores = ncores)
  
  pos <- match(dimnames, names(dim(BiasCorrected)))
  BiasCorrected <- aperm(BiasCorrected, pos)

  if (is.null(exp_cor)) {
    exp$data <- BiasCorrected
    exp$Datasets <- c(exp$Datasets, obs$Datasets)
    exp$source_files <- c(exp$source_files, obs$source_files)
    
    return(exp)

  } else {
    exp_cor$data <- BiasCorrected
    exp_cor$Datasets <- c(exp_cor$Datasets, exp$Datasets, obs$Datasets)
    exp_cor$source_files <- c(exp_cor$source_files, exp$source_files, obs$source_files)

    return(exp_cor)
  }
}

#' Bias Correction based on the mean and standard deviation adjustment
#'
#'@author Verónica Torralba, \email{veronica.torralba@bsc.es}
#'@description This function applies the simple bias adjustment technique 
#'described in Torralba et al. (2017). The adjusted forecasts have an equivalent 
#'standard deviation and mean to that of the reference dataset.
#'
#'@param exp A multidimensional array with named dimensions containing the 
#'  seasonal forecast experiment data with at least 'member' and 'sdate' 
#'  dimensions.
#'@param obs A multidimensional array with named dimensions containing the 
#'  observed data with at least 'sdate' dimension.
#'@param exp_cor A multidimensional array with named dimensions containing the 
#'  seasonl forecast experiment to be corrected. If it is NULL, the 'exp' 
#'  forecast will be corrected.
#'@param na.rm A logical value indicating whether missing values should be 
#'  stripped before the computation proceeds, by default it is set to FALSE.
#'@param memb_dim A character string indicating the name of the member 
#'  dimension. By default, it is set to 'member'.
#'@param sdate_dim A character string indicating the name of the start date 
#'  dimension. By default, it is set to 'sdate'.
#'@param ncores An integer that indicates the number of cores for parallel 
#'  computations using multiApply function. The default value is NULL.
#'
#'@return An array containing the bias corrected forecasts with the same 
#'dimensions of the experimental data.
#'
#'@references Torralba, V., F.J. Doblas-Reyes, D. MacLeod, I. Christel and M. 
#'Davis (2017). Seasonal climate prediction: a new source of information for the 
#'management of wind energy resources. Journal of Applied Meteorology and 
#'Climatology, 56, 1231-1247, doi:10.1175/JAMC-D-16-0204.1. (CLIM4ENERGY, 
#'EUPORIAS, NEWA, RESILIENCE, SPECS)
#'
#'@examples
#'mod1 <- 1 : (1 * 3 * 4 * 5 * 6 * 7)
#'dim(mod1) <- c(dataset = 1, member = 3, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'obs1 <- 1 : (1 * 1 * 4 * 5 * 6 * 7)
#'dim(obs1) <- c(dataset = 1, member = 1, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'a <- BiasCorrection(exp = mod1, obs = obs1)
#'@import multiApply
#'@export
BiasCorrection <- function(exp, obs, exp_cor = NULL, na.rm = FALSE,
                           memb_dim = 'member', sdate_dim = 'sdate',
                           ncores = NULL) {
  # Check inputs
  ## exp, obs
  if (!is.array(exp) || !is.numeric(exp)) {
    stop("Parameter 'exp' must be a numeric array.")
  }
  if (!is.array(obs) || !is.numeric(obs)) {
    stop("Parameter 'obs' must be a numeric array.")
  }
  obsdims <- names(dim(obs))
  expdims <- names(dim(exp))
  if (is.null(expdims)) {
    stop("Parameter 'exp' must have dimension names.")
  }
  if (is.null(obsdims)) {
    stop("Parameter 'obs' must have dimension names.")
  }
  if (any(is.na(exp)))  {
    warning("Parameter 'exp' contains NA values.")
  }
  if (any(is.na(obs)))  {
    warning("Parameter 'obs' contains NA values.")
  }
  ## exp_cor
  if (!is.null(exp_cor)) {
    exp_cordims <- names(dim(exp_cor))
    if (is.null(exp_cordims)) {
      stop("Parameter 'exp_cor' must have dimension names.")
    }
  }
  ## sdate_dim, memb_dim
  if (!is.character(sdate_dim) || length(sdate_dim) != 1) {
    stop("Parameter 'sdate_dim' must be a character string.")
  }
  if (!sdate_dim %in% expdims || !sdate_dim %in% obsdims) {
    stop("Parameter 'sdate_dim' is not found in 'exp' or 'obs' dimension.")
  }
  if (dim(exp)[sdate_dim] == 1) {
    stop("Parameter 'exp' must have dimension length of 'sdate_dim' bigger than 1.")
  }
  if (!all(c(memb_dim, sdate_dim) %in% expdims)) {
    stop("Parameter 'exp' requires 'sdate_dim' and 'memb_dim' dimensions.")  
  }
  if (memb_dim %in% obsdims) {
    if (dim(obs)[memb_dim] != 1) {
      stop("If parameter 'obs' has dimension 'memb_dim' its length must be equal to 1.")
    }
  }
  ## na.rm
  if (!is.logical(na.rm))  {
    na.rm <- FALSE
    warning("Paramater 'na.rm' must be a logical, it has been set to FALSE.")
  }
  if (length(na.rm) > 1) {
    na.rm <- na.rm[1]
    warning("Paramter 'na.rm' has length greater than 1, and only the fist element is used.")
  }
  ## ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) | ncores %% 1 != 0 | ncores <= 0 |
      length(ncores) > 1) {
      stop("Parameter 'ncores' must be either NULL or a positive integer.")
    }
  }

  target_dims_obs <- sdate_dim
  if (memb_dim %in% names(dim(obs))) {
    target_dims_obs <- c(memb_dim, target_dims_obs)
  }
  
  if (is.null(exp_cor)) {
    BiasCorrected <- Apply(data = list(var_obs = obs, var_exp = exp),
                           target_dims = list(target_dims_obs,
                                              c(memb_dim, sdate_dim)),
                           fun = .sbc, 
                           na.rm = na.rm, ncores = ncores)$output1
  } else {
    BiasCorrected <- Apply(data = list(var_obs = obs,
                                       var_exp = exp,
                                       var_cor = exp_cor),
                           target_dims = list(target_dims_obs,
                                              c(memb_dim, sdate_dim),
                                              c(memb_dim, sdate_dim)),
                           fun = .sbc, 
                           output_dims = c(memb_dim, sdate_dim),
                           na.rm = na.rm, ncores = ncores)$output1
  }
  return(BiasCorrected)
}

.sbc <- function(var_obs, var_exp , var_cor = NULL, na.rm = FALSE) {

  ntime <- dim(var_exp)[2]
  corrected <- NA * var_exp

  if (is.null(var_cor)) {
    for (t in 1:ntime) {
      # parameters
      sd_obs <- sd(var_obs[-t], na.rm = na.rm)
      sd_exp <- sd(var_exp[, -t], na.rm = na.rm)
      clim_exp <- mean(var_exp[, -t], na.rm = na.rm)
      clim_obs <- mean(var_obs[-t], na.rm = na.rm)
      
      # bias corrected forecast
      corrected[, t] <- ((var_exp[, t] - clim_exp) * (sd_obs / sd_exp)) + clim_obs
    }
  } else {
    # parameters
    sd_obs <- sd(var_obs, na.rm = na.rm)
    sd_exp <- sd(var_exp, na.rm = na.rm)
    clim_exp <- mean(var_exp, na.rm = na.rm)
    clim_obs <- mean(var_obs, na.rm = na.rm)
    
    # bias corrected forecast
    corrected <- ((var_cor - clim_exp) * (sd_obs / sd_exp)) + clim_obs
  }

  return(corrected)
}
