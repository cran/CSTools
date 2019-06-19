#'Forecast Calibration based on the ensemble inflation
#'
#'@author Verónica Torralba, \email{veronica.torralba@bsc.es}
#'@description This function applies a variance inflation technique described in Doblas-Reyes et al. (2005) in leave-one-out cross-validation. This bias adjustment method produces calibrated forecasts with equivalent mean and variance to that of the reference dataset, but at the same time preserve reliability. 
#'
#'@references Doblas-Reyes F.J, Hagedorn R, Palmer T.N. The rationale behind the success of multi-model ensembles in seasonal forecasting-II calibration and combination. Tellus A. 2005;57:234-252. doi:10.1111/j.1600-0870.2005.00104.x
#'
#'@param exp an object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the seasonal forecast experiment data in the element named \code{$data}.
#'@param obs an object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the observed data in the element named \code{$data}.
#'
#'@return an object of class \code{s2dv_cube} containing the calibrated forecasts in the element \code{$data} with the same dimensions of the experimental data.
#'
#'@import s2dverification
#'@import multiApply
#'
#'@seealso \code{\link{CST_Load}}
#'
#'@examples
#'# Example
#'# Load data using CST_Load or use the sample data provided:
#'library(zeallot)
#'c(exp, obs) %<-% areave_data
#'exp_calibrated <- CST_Calibration(exp = exp, obs = obs)
#'str(exp_calibrated)
#'@export
CST_Calibration <- function(exp, obs) {
  if (!inherits(exp, 's2dv_cube') || !inherits(exp, 's2dv_cube')) {
    stop("Parameter 'exp' and 'obs' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load.")
  }

  if (dim(obs$data)['member'] != 1) {
    stop("The length of the dimension 'member' in the component 'data' ",
         "of the parameter 'obs' must be equal to 1.")
  }
  dimnames <- names(dim(exp$data))
  Calibrated <- Calibration(exp = exp$data, obs = obs$data)
  pos <- match(dimnames, names(dim(Calibrated)))
  Calibrated <- aperm(Calibrated, pos)
  names(dim(Calibrated)) <- dimnames 
  exp$data <- Calibrated
  exp$Datasets <- c(exp$Datasets, obs$Datasets)
  exp$source_files <- c(exp$source_files, obs$source_files)
  return(exp)
}

Calibration <- function(exp, obs) {
  if (!all(c('member', 'sdate') %in% names(dim(exp)))) {
    stop("Parameter 'exp' must have the dimensions 'member' and 'sdate'.")
  }

  if (!all(c('sdate') %in% names(dim(obs)))) {
    stop("Parameter 'obs' must have the dimension 'sdate'.")
  }

  if (any(is.na(exp)))  {
    warning("Parameter 'exp' contains NA values.")
  }
  
  if (any(is.na(obs))) {
    warning("Parameter 'obs' contains NA values.")
  }

  target_dims_obs <- 'sdate'
  if ('member' %in% names(dim(obs))) {
    target_dims_obs <- c('member', target_dims_obs)
  }

  Calibrated <- Apply(data = list(var_obs = obs, var_exp = exp),
                      target_dims = list(target_dims_obs, c('member', 'sdate')),
                      fun = .cal)$output1

  return(Calibrated)
}

.cal <- function(var_obs, var_exp) {
  ntime <- dim(var_exp)[which(names(dim(var_exp)) == 'sdate')][]
  nmembers <- dim(var_exp)[which(names(dim(var_exp)) == 'member')][]

  if (all(names(dim(var_exp)) != c('member','sdate'))) {
    var_exp <- t(var_exp)
  }
  climObs <- NA * var_obs
  climPred <- NA * var_obs
  for (t in 1 : length(var_obs))
  {
    climObs[t] <- mean(var_obs[ , -t])
    climPred[t] <- mean(var_exp[ , -t])
  }
  var_obs <- var_obs - climObs
 
  calibrated <- NA * var_exp
  for (t in 1 : ntime) {
    # defining forecast,hindcast and observation in cross-validation
    fcst <- NA * var_exp[ , t]
    hcst <- NA * var_exp[ , -t]
    for (i in 1 : nmembers) {
      fcst[i] <- var_exp[i, t] - climPred[t]
      hcst[i, ] <- var_exp[i, -t]- climPred[t]
    }
    obs <- var_obs[-t]
    #coefficients
    em_fcst <- mean(fcst)
    em_hcst <- apply(hcst, c(2), mean)
    corr <- cor(em_hcst, obs)
    sd_obs <- sd(obs)
    sd_em_hcst <- sd(em_hcst)
    
    fcst_diff <- fcst - em_fcst
    hcst_diff <- NA * hcst
    for (n in 1 : nmembers) {
      hcst_diff[n,] <- hcst[n,] - em_hcst
    }
    sd_hcst_diff <- sd(hcst_diff)
    
    a <- corr * (sd_obs / sd_em_hcst)
    b <- (sd_obs / sd_hcst_diff) * sqrt(1 - (corr ^ 2))
    
    calibrated[, t] <- (a * em_fcst) + (b * fcst_diff) + climObs[t]
  }
  names(dim(calibrated)) <- c('member', 'sdate')
  return(calibrated)
}

