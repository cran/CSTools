#'Forecast Calibration 
#'
#'@author Verónica Torralba, \email{veronica.torralba@bsc.es} 
#'@author Bert Van Schaeybroeck, \email{bertvs@meteo.be}
#'@description Equivalent to function \code{Calibration} but for objects of class \code{s2dv_cube}.
#'
#'@param exp an object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the seasonal hindcast experiment data in the element named \code{$data}. The hindcast is used to calibrate the forecast in case the forecast is provided; if not, the same hindcast will be calibrated instead.
#'@param obs an object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the observed data in the element named \code{$data}.
#'@param exp_cor an optional object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the seasonal forecast experiment data in the element named \code{$data}. If the forecast is provided, it will be calibrated using the hindcast and observations; if not, the hindcast will be calibrated instead.
#'@param cal.method is the calibration method used, can be either \code{bias}, \code{evmos}, \code{mse_min}, \code{crps_min} or \code{rpc-based}. Default value is \code{mse_min}.
#'@param eval.method is the sampling method used, can be either \code{in-sample} or \code{leave-one-out}. Default value is the \code{leave-one-out} cross validation. In case the forecast is provided, any chosen eval.method is over-ruled and a third option is used.
#'@param multi.model is a boolean that is used only for the \code{mse_min} method. If multi-model ensembles or ensembles of different sizes are used, it must be set to \code{TRUE}. By default it is \code{FALSE}. Differences between the two approaches are generally small but may become large when using small ensemble sizes. Using multi.model when the calibration method is \code{bias}, \code{evmos} or \code{crps_min} will not affect the result.
#'@param na.fill is a boolean that indicates what happens in case calibration is not possible or will yield unreliable results. This happens when three or less forecasts-observation pairs are available to perform the training phase of the calibration. By default \code{na.fill} is set to true such that NA values will be returned. If \code{na.fill} is set to false, the uncorrected data will be returned. 
#'@param na.rm is a boolean that indicates whether to remove the NA values or not. The default value is \code{TRUE}. See Details section for further information about its use and compatibility with \code{na.fill}.
#'@param apply_to is a character string that indicates whether to apply the calibration to all the forecast (\code{"all"}) or only to those where the correlation between the ensemble mean and the observations is statistically significant (\code{"sign"}). Only useful if \code{cal.method == "rpc-based"}.
#'@param alpha is a numeric value indicating the significance level for the correlation test. Only useful if \code{cal.method == "rpc-based" & apply_to == "sign"}.
#'@param memb_dim is a character string indicating the name of the member dimension. By default, it is set to 'member'.
#'@param sdate_dim is a character string indicating the name of the start date dimension. By default, it is set to 'sdate'.
#'@param ncores is an integer that indicates the number of cores for parallel computations using multiApply function. The default value is one.
#'@return an object of class \code{s2dv_cube} containing the calibrated forecasts in the element \code{$data} with the same dimensions as the one in the exp object.
#'
#'@importFrom s2dv InsertDim
#'@import abind
#'
#'@seealso \code{\link{CST_Load}}
#'
#'@examples
#'# Example 1:
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
#'a <- CST_Calibration(exp = exp, obs = obs, cal.method = "mse_min", eval.method = "in-sample")
#'str(a)
#'
#'# Example 2:
#'mod1 <- 1 : (1 * 3 * 4 * 5 * 6 * 7)
#'mod2 <- 1 : (1 * 3 * 1 * 5 * 6 * 7)
#'dim(mod1) <- c(dataset = 1, member = 3, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'dim(mod2) <- c(dataset = 1, member = 3, sdate = 1, ftime = 5, lat = 6, lon = 7)
#'obs1 <- 1 : (1 * 1 * 4 * 5 * 6 * 7)
#'dim(obs1) <- c(dataset = 1, member = 1, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'lon <- seq(0, 30, 5)
#'lat <- seq(0, 25, 5)
#'exp <- list(data = mod1, lat = lat, lon = lon)
#'obs <- list(data = obs1, lat = lat, lon = lon)
#'exp_cor <- list(data = mod2, lat = lat, lon = lon)
#'attr(exp, 'class') <- 's2dv_cube'
#'attr(obs, 'class') <- 's2dv_cube'
#'attr(exp_cor, 'class') <- 's2dv_cube'
#'a <- CST_Calibration(exp = exp, obs = obs, exp_cor = exp_cor, cal.method = "evmos")
#'str(a)
#'@export

CST_Calibration <- function(exp, obs, exp_cor = NULL, cal.method = "mse_min", 
                            eval.method = "leave-one-out", multi.model = FALSE, 
                            na.fill = TRUE, na.rm = TRUE, apply_to = NULL, alpha = NULL,
                            memb_dim = 'member', sdate_dim = 'sdate', ncores = 1) {
  
  if(!missing(multi.model) & !(cal.method == "mse_min")){
	  warning(paste0("The multi.model parameter is ignored when using the calibration method ", cal.method))
  }
  
  if(is.null(exp_cor)){ #exp will be used to calibrate and will also be calibrated: "calibrate hindcast"
    if (!inherits(exp, "s2dv_cube") || !inherits(obs, "s2dv_cube")) {
        stop("Parameter 'exp' and 'obs' must be of the class 's2dv_cube', ",
            "as output by CSTools::CST_Load.")
    }
    exp$data <- Calibration(exp = exp$data, obs = obs$data, exp_cor = NULL,
      cal.method = cal.method, 
      eval.method = eval.method,  
      multi.model =  multi.model, 
      na.fill = na.fill, na.rm = na.rm, 
      apply_to = apply_to, alpha = alpha,
      memb_dim = memb_dim, sdate_dim = sdate_dim,
      ncores = ncores)
    exp$Datasets <- c(exp$Datasets, obs$Datasets)
    exp$source_files <- c(exp$source_files, obs$source_files)
   
    return(exp)  

  }else{ #if exp_cor is provided, it will be calibrated: "calibrate forecast instead of hindcast"
    eval.method = "hindcast-vs-forecast" #if exp_cor is provided, eval.method is overrruled (because if exp_cor is provided, the train data will be all data of "exp" and the evalutaion data will be all data of "exp_cor"; no need for "leave-one-out" or "in-sample")
    if (!inherits(exp, "s2dv_cube") || !inherits(obs, "s2dv_cube") || !inherits(exp_cor, "s2dv_cube")) {
        stop("Parameter 'exp', 'obs' and 'exp_cor' must be of the class 's2dv_cube', ",
            "as output by CSTools::CST_Load.")
    }
    exp_cor$data <- Calibration(exp = exp$data, obs = obs$data, exp_cor = exp_cor$data,
      cal.method = cal.method, 
      eval.method = eval.method,  
      multi.model =  multi.model, 
      na.fill = na.fill, na.rm = na.rm, 
      apply_to = apply_to, alpha = alpha,
      memb_dim = memb_dim, sdate_dim = sdate_dim,
      ncores = ncores)
    exp_cor$Datasets <- c(exp_cor$Datasets, obs$Datasets)
    exp_cor$source_files <- c(exp_cor$source_files, exp$source_files, obs$source_files)  
      
    return(exp_cor)

  }  
}



#'Forecast Calibration 
#'
#'@author Verónica Torralba, \email{veronica.torralba@bsc.es} 
#'@author Bert Van Schaeybroeck, \email{bertvs@meteo.be}
#'@description Five types of member-by-member bias correction can be performed. The \code{"bias"} method corrects the bias only, the \code{"evmos"} method applies a variance inflation technique to ensure the correction of the bias and the correspondence of variance between forecast and observation (Van Schaeybroeck and Vannitsem, 2011). The ensemble calibration methods \code{"mse_min"} and \code{"crps_min"} correct the bias, the overall forecast variance and the ensemble spread as described in Doblas-Reyes et al. (2005) and Van Schaeybroeck and Vannitsem (2015), respectively. While the \code{"mse_min"} method minimizes a constrained mean-squared error using three parameters, the \code{"crps_min"} method features four parameters and minimizes the Continuous Ranked Probability Score (CRPS). The \code{"rpc-based"} method adjusts the forecast variance ensuring that the ratio of predictable components (RPC) is equal to one, as in Eade et al. (2014).
#'@description Both in-sample or our out-of-sample (leave-one-out cross validation) calibration are possible.
#'@references Doblas-Reyes F.J, Hagedorn R, Palmer T.N. The rationale behind the success of multi-model ensembles in seasonal forecasting-II calibration and combination. Tellus A. 2005;57:234-252. doi:10.1111/j.1600-0870.2005.00104.x
#'@references Eade, R., Smith, D., Scaife, A., Wallace, E., Dunstone, N., Hermanson, L., & Robinson, N. (2014). Do seasonal-to-decadal climate predictions underestimate the predictability of the read world? Geophysical Research Letters, 41(15), 5620-5628. doi: 10.1002/2014GL061146
#'@references Van Schaeybroeck, B., & Vannitsem, S. (2011). Post-processing through linear regression. Nonlinear Processes in Geophysics, 18(2), 147. doi:10.5194/npg-18-147-2011
#'@references Van Schaeybroeck, B., & Vannitsem, S. (2015). Ensemble post-processing using member-by-member approaches: theoretical aspects. Quarterly Journal of the Royal Meteorological Society, 141(688), 807-818.  doi:10.1002/qj.2397
#'
#'@param exp a multidimensional array with named dimensions (at least 'sdate' and 'member') containing the seasonal hindcast experiment data. The hindcast is used to calibrate the forecast in case the forecast is provided; if not, the same hindcast will be calibrated instead.
#'@param obs a multidimensional array with named dimensions (at least 'sdate') containing the observed data.
#'@param exp_cor an optional multidimensional array with named dimensions (at least 'sdate' and 'member') containing the seasonal forecast experiment data. If the forecast is provided, it will be calibrated using the hindcast and observations; if not, the hindcast will be calibrated instead.
#'@param cal.method is the calibration method used, can be either \code{bias}, \code{evmos}, \code{mse_min}, \code{crps_min} or \code{rpc-based}. Default value is \code{mse_min}.
#'@param eval.method is the sampling method used, can be either \code{in-sample} or \code{leave-one-out}. Default value is the \code{leave-one-out} cross validation. In case the forecast is provided, any chosen eval.method is over-ruled and a third option is used.
#'@param multi.model is a boolean that is used only for the \code{mse_min} method. If multi-model ensembles or ensembles of different sizes are used, it must be set to \code{TRUE}. By default it is \code{FALSE}. Differences between the two approaches are generally small but may become large when using small ensemble sizes. Using multi.model when the calibration method is \code{bias}, \code{evmos} or \code{crps_min} will not affect the result.
#'@param na.fill is a boolean that indicates what happens in case calibration is not possible or will yield unreliable results. This happens when three or less forecasts-observation pairs are available to perform the training phase of the calibration. By default \code{na.fill} is set to true such that NA values will be returned. If \code{na.fill} is set to false, the uncorrected data will be returned. 
#'@param na.rm is a boolean that indicates whether to remove the NA values or not. The default value is \code{TRUE}.
#'@param apply_to is a character string that indicates whether to apply the calibration to all the forecast (\code{"all"}) or only to those where the correlation between the ensemble mean and the observations is statistically significant (\code{"sign"}). Only useful if \code{cal.method == "rpc-based"}.
#'@param alpha is a numeric value indicating the significance level for the correlation test. Only useful if \code{cal.method == "rpc-based" & apply_to == "sign"}.
#'@param memb_dim is a character string indicating the name of the member dimension. By default, it is set to 'member'.
#'@param sdate_dim is a character string indicating the name of the start date dimension. By default, it is set to 'sdate'.
#'@param ncores is an integer that indicates the number of cores for parallel computations using multiApply function. The default value is one.
#'@return an array containing the calibrated forecasts with the same dimensions as the \code{exp} array.
#'
#'@importFrom s2dv InsertDim MeanDims Reorder
#'@import abind
#'@import multiApply
#'@importFrom ClimProjDiags Subset
#'
#'@seealso \code{\link{CST_Load}}
#'
#'@details 
#'Both the \code{na.fill} and \code{na.rm} parameters can be used to indicate how the function has to handle the NA values. The \code{na.fill} parameter checks whether there are more than three forecast-observations pairs to perform the computation. In case there are three or less pairs, the computation is not carried out, and the value returned by the function depends on the value of this parameter (either NA if \code{na.fill == TRUE} or the uncorrected value if \code{na.fill == TRUE}). On the other hand, \code{na.rm} is used to indicate the function whether to remove the missing values during the computation of the parameters needed to perform the calibration.
#'
#'@examples
#'mod1 <- 1 : (1 * 3 * 4 * 5 * 6 * 7)
#'dim(mod1) <- c(dataset = 1, member = 3, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'obs1 <- 1 : (1 * 1 * 4 * 5 * 6 * 7)
#'dim(obs1) <- c(dataset = 1, member = 1, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'a <- Calibration(exp = mod1, obs = obs1)
#'str(a)
#'@export
Calibration <- function(exp, obs, exp_cor=NULL, cal.method = "mse_min",
                        eval.method = "leave-one-out",  
                        multi.model = FALSE, na.fill = TRUE, 
                        na.rm = TRUE, apply_to = NULL, alpha = NULL,
                        memb_dim = 'member', sdate_dim = 'sdate', ncores = 1) {
	
  dim.exp <- dim(exp)
  amt.dims.exp <- length(dim.exp)
  dim.obs <- dim(obs)
  amt.dims.obs <- length(dim.obs)
  dim.names.exp <- names(dim.exp)
  dim.names.obs <- names(dim.obs)
  if(!is.null(exp_cor)){
    dim.exp_cor <- dim(exp_cor)
    amt.dims.exp_cor <- length(dim.exp_cor)
    dim.names.exp_cor <- names(dim.exp_cor)
  }
  if (is.null(memb_dim) || !is.character(memb_dim)) {
    stop("Parameter 'memb_dim' should be a character string indicating the",
         "name of the dimension where members are stored in 'exp'.")
  }
  if (length(memb_dim) > 1) {
    memb_dim <- memb_dim[1]
    warning("Parameter 'memb_dim' has length greater than 1 and only",
            " the first element will be used.")
  } 
  
  if (is.null(sdate_dim) || !is.character(sdate_dim)) {
    stop("Parameter 'sdate_dim' should be a character string indicating the",
         "name of the dimension where start dates are stored in 'exp'.")    
  }
  if (length(sdate_dim) > 1) {
    sdate_dim <- sdate_dim[1]
    warning("Parameter 'sdate_dim' has length greater than 1 and only",
            " the first element will be used.")
  }
  target.dim.names.exp <- c(memb_dim, sdate_dim)
  target.dim.names.obs <- sdate_dim
  
  if (!all(target.dim.names.exp %in% dim.names.exp)) {
    stop("Parameter 'exp' must have the dimensions defined in memb_dim ",
         "and sdate_dim.")
  }

  if(!is.null(exp_cor)){
    if (!all(target.dim.names.exp %in% dim.names.exp_cor)) {
      stop("Parameter 'exp_cor' must have the dimensions defined in memb_dim ",
           "and sdate_dim.")
    }
  }

  if (!all(c(sdate_dim) %in% dim.names.obs)) {
    stop("Parameter 'obs' must have the dimension defined in sdate_dim ",
         "parameter.")
  }

  if (any(is.na(exp)))  {
    warning("Parameter 'exp' contains NA values.")
  }
  
  if(!is.null(exp_cor)){
    if (any(is.na(exp_cor)))  {
      warning("Parameter 'exp_cor' contains NA values.")
    }
  }

  if (any(is.na(obs))) {
    warning("Parameter 'obs' contains NA values.")
  }
  
  if (memb_dim %in% names(dim(obs))) {  
    obs <- Subset(obs, along = memb_dim, indices = 1, drop = "selected")
  }
  data.set.sufficiently.large.out <- 
    Apply(data = list(exp = exp, obs = obs),
      target_dims = list(exp = target.dim.names.exp, obs = target.dim.names.obs),
      ncores = ncores,
      fun = .data.set.sufficiently.large)$output1
  
  if(!all(data.set.sufficiently.large.out)){		
  	if(na.fill){
        warning("Some forecast data could not be corrected due to data lack",
                " and is replaced with NA values")
  	} else {
        warning("Some forecast data could not be corrected due to data lack",
                " and is replaced with uncorrected values")
  	 }
  }
  
  if (!na.rm %in% c(TRUE,FALSE)) {
    stop("Parameter 'na.rm' must be TRUE or FALSE.")
  }
  if (cal.method == 'rpc-based') {
    if (is.null(apply_to)) {
      apply_to <- 'sign'
      warning("'apply_to' cannot be NULL for 'rpc-based' method so it has been set to 'sign', as in Eade et al. (2014).")
    } else if (!apply_to %in% c('all','sign')) {
      stop("'apply_to' must be either 'all' or 'sign' when 'rpc-based' method is used.")
    }
    if (apply_to == 'sign') {
      if (is.null(alpha)) {
        alpha <- 0.1
        warning("'alpha' cannot be NULL for 'rpc-based' method so it has been set to 0.1, as in Eade et al. (2014).")
      } else if (!is.numeric(alpha) | alpha <= 0 | alpha >= 1) {
        stop("'alpha' must be a number between 0 and 1.")
      }
    }
  }
  
  if(is.null(exp_cor)){
    calibrated <- Apply(data = list(exp = exp, obs = obs),
        cal.method = cal.method,
        eval.method = eval.method,
        multi.model = multi.model,
        na.fill = na.fill, na.rm = na.rm, 
        apply_to = apply_to, alpha = alpha,
        target_dims = list(exp = target.dim.names.exp, obs = target.dim.names.obs),
        ncores = ncores, output_dims = target.dim.names.exp,
        fun = .cal)$output1
      dexes <- match(names(dim(exp)), names(dim(calibrated)))
      calibrated <- aperm(calibrated, dexes)
      dimnames(calibrated) <- dimnames(exp)[dexes]
  }else{
    calibrated <- Apply(data = list(exp = exp, obs = obs, exp_cor = exp_cor),
      cal.method = cal.method,
      eval.method = eval.method,
      multi.model = multi.model,
      na.fill = na.fill, na.rm = na.rm, 
      apply_to = apply_to, alpha = alpha,
      target_dims = list(exp = target.dim.names.exp, obs = target.dim.names.obs, exp_cor = target.dim.names.exp),
      ncores = ncores, output_dims = target.dim.names.exp,
      fun = .cal)$output1
    dexes <- match(names(dim(exp_cor)), names(dim(calibrated)))
    calibrated <- aperm(calibrated, dexes)
    dimnames(calibrated) <- dimnames(exp_cor)[dexes]
  }    

  return(calibrated)
}


.data.set.sufficiently.large <- function(exp, obs){
  amt.min.samples <- 3
  amt.good.pts <- sum(!is.na(obs) & !apply(exp, c(2), function(x) all(is.na(x))))
  return(amt.good.pts > amt.min.samples)
}

.make.eval.train.dexes <- function(eval.method, amt.points, amt.points_cor){ 
  if(eval.method == "leave-one-out"){
    dexes.lst <- lapply(seq(1, amt.points), function(x) return(list(eval.dexes = x,
                        train.dexes = seq(1, amt.points)[-x])))
  } else if (eval.method == "in-sample"){
    dexes.lst <- list(list(eval.dexes = seq(1, amt.points), 
                           train.dexes = seq(1, amt.points)))
  } else if (eval.method == "hindcast-vs-forecast"){
    dexes.lst <- list(list(eval.dexes = seq(1,amt.points_cor),
                           train.dexes = seq(1, amt.points)))
  } else {
    stop(paste0("unknown sampling method: ",eval.method))
  }
  return(dexes.lst)
}

.cal <- function(exp, obs, exp_cor = NULL, cal.method, eval.method, multi.model, na.fill, na.rm, apply_to, alpha) {
  if(is.null(exp_cor)){
    exp_cor <- exp ## generate a copy of exp so that the same function can run 
                   ## when exp_cor is provided and when it's not
  }
  obs <- as.vector(obs)
  dims.fc <- dim(exp)
  dims.fc_cor <- dim(exp_cor) ## new line
  amt.mbr <- dims.fc[1]
  amt.sdate <- dims.fc[2]
  amt.sdate_cor <- dims.fc_cor[2] ## new line
  var.cor.fc <- NA * exp_cor ## modified line (exp_cor instead of exp); 
                             ## in case of exp_cor not provided, at this point exp_cor 
                             ## is already the same as exp so no change
  names(dim(var.cor.fc)) <- dims.fc
  
  if(!.data.set.sufficiently.large(exp = exp, obs = obs)){
  	if(na.fill){
   	  return(var.cor.fc)
  	} else {
   	  var.cor.fc[] <- exp[]
  	  return(var.cor.fc)
  	}
  }
  eval.train.dexeses <- .make.eval.train.dexes(eval.method, amt.points = amt.sdate, 
                                               amt.points_cor = amt.sdate_cor)
  amt.resamples <- length(eval.train.dexeses)
  for (i.sample in seq(1, amt.resamples)) {
    # defining training (tr) and evaluation (ev) subsets 
    eval.dexes <- eval.train.dexeses[[i.sample]]$eval.dexes
    train.dexes <- eval.train.dexeses[[i.sample]]$train.dexes
    
    fc.ev <- exp_cor[ , eval.dexes, drop = FALSE] ## modified line (exp_cor instead of exp)
        ##  fc.ev is used to evaluate (not train; train should be done with exp (hindcast))
    fc.tr <- exp[ , train.dexes]
    obs.tr <- obs[train.dexes , drop = FALSE] 
    
    if(cal.method == "bias"){
	    var.cor.fc[ , eval.dexes] <- fc.ev + mean(obs.tr, na.rm = na.rm) - mean(fc.tr, na.rm = na.rm)
	  } else if(cal.method == "evmos"){ # forecast correction implemented
	    #calculate ensemble and observational characteristics
	    quant.obs.fc.tr <- .calc.obs.fc.quant(obs = obs.tr, fc = fc.tr, na.rm = na.rm)
      #calculate value for regression parameters
      init.par <- c(.calc.evmos.par(quant.obs.fc.tr, na.rm = na.rm))
	    #correct evaluation subset
      var.cor.fc[ , eval.dexes] <- .correct.evmos.fc(fc.ev , init.par, na.rm = na.rm)
	  } else if (cal.method == "mse_min"){
	    #calculate ensemble and observational characteristics
	    quant.obs.fc.tr <- .calc.obs.fc.quant(obs = obs.tr, fc = fc.tr, na.rm = na.rm)
      #calculate value for regression parameters
      init.par <- .calc.mse.min.par(quant.obs.fc.tr, multi.model, na.rm = na.rm)
	    #correct evaluation subset
      var.cor.fc[ , eval.dexes] <- .correct.mse.min.fc(fc.ev , init.par, na.rm = na.rm)      
    } else if (cal.method == "crps_min"){
	    #calculate ensemble and observational characteristics
	    quant.obs.fc.tr <- .calc.obs.fc.quant.ext(obs = obs.tr, fc = fc.tr, na.rm = na.rm)
      #calculate initial value for regression parameters
      init.par <- c(.calc.mse.min.par(quant.obs.fc.tr, na.rm = na.rm), 0.001)
      init.par[3] <- sqrt(init.par[3])
      #calculate regression parameters on training dataset
      optim.tmp <- optim(par = init.par, 
        fn = .calc.crps.opt, 
        gr = .calc.crps.grad.opt, 
        quant.obs.fc = quant.obs.fc.tr,
        na.rm = na.rm,
        method = "BFGS")
      
      mbm.par <- optim.tmp$par
	    #correct evaluation subset
      var.cor.fc[ , eval.dexes] <- .correct.crps.min.fc(fc.ev , mbm.par, na.rm = na.rm)
    } else if (cal.method == 'rpc-based') {
      ens_mean.ev <- multiApply::Apply(data = fc.ev, target_dims = names(amt.mbr), fun = mean, na.rm = na.rm)$output1 ## Ensemble mean
      ens_mean.tr <- multiApply::Apply(data = fc.tr, target_dims = names(amt.mbr), fun = mean, na.rm = na.rm)$output1 ## Ensemble mean
      ens_spread.tr <- multiApply::Apply(data = list(fc.tr, ens_mean.tr), target_dims = names(amt.sdate), fun = "-")$output1 ## Ensemble spread
      exp_mean.tr <- mean(fc.tr, na.rm = na.rm) ## Mean (climatology)
      var_signal.tr <- var(ens_mean.tr, na.rm = na.rm) ## Ensemble mean variance
      var_noise.tr <- var(as.vector(ens_spread.tr), na.rm = na.rm) ## Variance of ensemble members about ensemble mean (= spread)
      var_obs.tr <- var(obs.tr, na.rm = na.rm) ## Variance in the observations
      r.tr <- cor(x = ens_mean.tr, y = obs.tr, method = 'pearson', use = ifelse(test = isTRUE(na.rm), yes = "pairwise.complete.obs", no = "everything")) ## Correlation between observations and the ensemble mean
      if ((apply_to == 'all') || (apply_to == 'sign' && cor.test(ens_mean.tr, obs.tr, method = 'pearson', alternative = 'greater')$p.value < alpha)) {
        ens_mean_cal <- (ens_mean.ev - exp_mean.tr) * r.tr * sqrt(var_obs.tr) / sqrt(var_signal.tr) + exp_mean.tr
        var.cor.fc[ , eval.dexes] <- s2dv::Reorder(data = multiApply::Apply(data = list(exp = fc.ev, ens_mean = ens_mean.ev, ens_mean_cal = ens_mean_cal), target_dims = names(amt.sdate), fun = .CalibrationMembersRPC, var_obs = var_obs.tr, var_noise = var_noise.tr, r = r.tr)$output1, 
                                                   order = names(dims.fc))
        dim(var.cor.fc) <- dims.fc
      } else { ## no significant -> replacing with observed climatology
        var.cor.fc[ , eval.dexes] <- array(data = mean(obs.tr, na.rm = na.rm), dim = dim(fc.ev))
      }
    } else {
	    stop("unknown calibration method: ",cal.method)
    }
  }
  return(var.cor.fc) 
}

.calc.obs.fc.quant <- function(obs, fc, na.rm){ #function to calculate different quantities of a series of ensemble forecasts and corresponding observations
  amt.mbr <- dim(fc)[1]
  obs.per.ens <- InsertDim(obs, posdim = 1, lendim = amt.mbr)
  fc.ens.av <- apply(fc, c(2), mean, na.rm = na.rm)
  cor.obs.fc <- cor(fc.ens.av, obs, use = "complete.obs")
  obs.av <- mean(obs, na.rm = na.rm)
  obs.sd <- sd(obs, na.rm = na.rm)
  return(
    append(
      .calc.fc.quant(fc = fc, na.rm = na.rm),
      list(
        obs.per.ens = obs.per.ens,
        cor.obs.fc = cor.obs.fc,
        obs.av = obs.av,
        obs.sd = obs.sd
      )
    )
  )
}

.calc.obs.fc.quant.ext <- function(obs, fc, na.rm){ #extended function to calculate different quantities of a series of ensemble forecasts and corresponding observations
  amt.mbr <- dim(fc)[1]
  obs.per.ens <- InsertDim(obs, posdim = 1, lendim = amt.mbr)
  fc.ens.av <- apply(fc, c(2), mean, na.rm = na.rm)
  cor.obs.fc <- cor(fc.ens.av, obs, use = "complete.obs")
  obs.av <- mean(obs, na.rm = na.rm)
  obs.sd <- sd(obs, na.rm = na.rm)

  return(
    append(
      .calc.fc.quant.ext(fc = fc, na.rm = na.rm),
      list(
        obs.per.ens = obs.per.ens,
        cor.obs.fc = cor.obs.fc,
        obs.av = obs.av,
        obs.sd = obs.sd
      )
    )
  )
}


.calc.fc.quant <- function(fc, na.rm){ #function to calculate different quantities of a series of ensemble forecasts
  amt.mbr <- dim(fc)[1]
  fc.ens.av <- apply(fc, c(2), mean, na.rm = na.rm)
  fc.ens.av.av <- mean(fc.ens.av, na.rm = na.rm)
  fc.ens.av.sd <- sd(fc.ens.av, na.rm = na.rm)
  fc.ens.av.per.ens <- InsertDim(fc.ens.av, posdim = 1, lendim = amt.mbr)
  fc.ens.sd <- apply(fc, c(2), sd, na.rm = na.rm)
  fc.ens.var.av.sqrt <- sqrt(mean(fc.ens.sd^2, na.rm = na.rm))
  fc.dev <- fc - fc.ens.av.per.ens
  fc.dev.sd <- sd(fc.dev, na.rm = na.rm)
  fc.av <- mean(fc, na.rm = na.rm)
  fc.sd <- sd(fc, na.rm = na.rm)
  return(
    list(
      fc.ens.av = fc.ens.av,
      fc.ens.av.av = fc.ens.av.av,
      fc.ens.av.sd = fc.ens.av.sd,
      fc.ens.av.per.ens = fc.ens.av.per.ens,
      fc.ens.sd = fc.ens.sd,
      fc.ens.var.av.sqrt = fc.ens.var.av.sqrt,
      fc.dev = fc.dev,
      fc.dev.sd = fc.dev.sd,
      fc.av = fc.av,
      fc.sd = fc.sd
    )
  )
}

.calc.fc.quant.ext <- function(fc, na.rm){ #extended function to calculate different quantities of a series of ensemble forecasts

  amt.mbr <- dim(fc)[1]
  repmat1.tmp <- InsertDim(fc, posdim = 1, lendim = amt.mbr)
  repmat2.tmp <- aperm(repmat1.tmp, c(2, 1, 3))
  spr.abs <- apply(abs(repmat1.tmp - repmat2.tmp), c(3), mean, na.rm = na.rm)
  spr.abs.per.ens <- InsertDim(spr.abs, posdim = 1, lendim = amt.mbr)

  return(
    append(.calc.fc.quant(fc, na.rm = na.rm),
	  list(spr.abs = spr.abs, spr.abs.per.ens = spr.abs.per.ens))
  )
}

#Below are the core or elementary functions to calculate the regression parameters for the different methods
.calc.mse.min.par <- function(quant.obs.fc, multi.model = F, na.rm){
  par.out <- rep(NA, 3)
  
  if(multi.model){
    par.out[3] <- with(quant.obs.fc, obs.sd * sqrt(1. - cor.obs.fc^2) / fc.ens.var.av.sqrt)
  } else {
    par.out[3] <- with(quant.obs.fc, obs.sd * sqrt(1. - cor.obs.fc^2) / fc.dev.sd)
  }
  par.out[2] <- with(quant.obs.fc, abs(cor.obs.fc) * obs.sd / fc.ens.av.sd)
  par.out[1] <- with(quant.obs.fc, obs.av - par.out[2] * fc.ens.av.av, na.rm = na.rm)
  
  return(par.out)
}
.calc.evmos.par <- function(quant.obs.fc, na.rm){
  par.out <- rep(NA, 2)
  par.out[2] <- with(quant.obs.fc, obs.sd / fc.sd)
  par.out[1] <- with(quant.obs.fc, obs.av - par.out[2] * fc.ens.av.av, na.rm = na.rm)
  return(par.out)
}
#Below are the core or elementary functions to calculate the functions necessary for the minimization of crps
.calc.crps.opt <- function(par, quant.obs.fc, na.rm){
  return( 
    with(quant.obs.fc, 
      mean(abs(obs.per.ens - (par[1] + par[2] * fc.ens.av.per.ens +
	    ((par[3])^2 + par[4] / spr.abs.per.ens) * fc.dev)), na.rm = na.rm) -
        mean(abs((par[3])^2 * spr.abs + par[4]) / 2., na.rm = na.rm)
    )
  )
}

.calc.crps.grad.opt <- function(par, quant.obs.fc, na.rm){
  sgn1 <- with(quant.obs.fc,sign(obs.per.ens - 
    (par[1] + par[2] * fc.ens.av.per.ens +
    ((par[3])^2 + par[4] / spr.abs.per.ens) * fc.dev)))
  sgn2 <- with(quant.obs.fc, sign((par[3])^2 + par[4] / spr.abs.per.ens))
  sgn3 <- with(quant.obs.fc,sign((par[3])^2 * spr.abs + par[4]))
  deriv.par1 <- mean(sgn1, na.rm = na.rm)
  deriv.par2 <- with(quant.obs.fc, mean(sgn1 * fc.dev, na.rm = na.rm))
  deriv.par3 <- with(quant.obs.fc, 
    mean(2* par[3] * sgn1 * sgn2 * fc.ens.av.per.ens, na.rm = na.rm) -
    mean(spr.abs * sgn3, na.rm = na.rm) / 2.)
  deriv.par4 <- with(quant.obs.fc,
    mean(sgn1 * sgn2 * fc.ens.av.per.ens / spr.abs.per.ens, na.rm = na.rm) -
    mean(sgn3, na.rm = na.rm) / 2.)
  return(c(deriv.par1, deriv.par2, deriv.par3, deriv.par4))
}

#Below are the core or elementary functions to correct the evaluation set based on the regression parameters
.correct.evmos.fc <- function(fc, par, na.rm){
  quant.fc.mp <- .calc.fc.quant(fc = fc, na.rm = na.rm)
  return(with(quant.fc.mp, par[1] + par[2] * fc))
}
.correct.mse.min.fc <- function(fc, par, na.rm){
  quant.fc.mp <- .calc.fc.quant(fc = fc, na.rm = na.rm)
  return(with(quant.fc.mp, par[1] + par[2] * fc.ens.av.per.ens + fc.dev * par[3]))
}
.correct.crps.min.fc <- function(fc, par, na.rm){
  quant.fc.mp <- .calc.fc.quant.ext(fc = fc, na.rm = na.rm)
  return(with(quant.fc.mp, par[1] + par[2] * fc.ens.av.per.ens + fc.dev * abs((par[3])^2 + par[4] / spr.abs)))
}

# Function to calibrate the individual members with the RPC-based method
.CalibrationMembersRPC <- function(exp, ens_mean, ens_mean_cal, var_obs, var_noise, r){
  member_cal <- (exp - ens_mean) * sqrt(var_obs) * sqrt(1 - r^2) / sqrt(var_noise) + ens_mean_cal
  return(member_cal)
}
