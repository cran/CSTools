#'Forecast Calibration 
#'
#'@author Verónica Torralba, \email{veronica.torralba@bsc.es} 
#'@author Bert Van Schaeybroeck, \email{bertvs@meteo.be}
#'@description Four types of member-by-member bias correction can be performed. The \code{bias} method corrects the bias only, the \code{evmos} method applies a variance inflation technique to ensure the correction of the bias and the correspondence of variance between forecast and observation (Van Schaeybroeck and Vannitsem, 2011). The ensemble calibration methods \code{"mse_min"} and \code{"crps_min"} correct the bias, the overall forecast variance and the ensemble spread as described in Doblas-Reyes et al. (2005) and Van Schaeybroeck and Vannitsem (2015), respectively. While the \code{"mse_min"} method minimizes a constrained mean-squared error using three parameters, the \code{"crps_min"} method features four parameters and minimizes the Continuous Ranked Probability Score (CRPS). 
#'@description Both in-sample or our out-of-sample (leave-one-out cross validation) calibration are possible.
#'@references Doblas-Reyes F.J, Hagedorn R, Palmer T.N. The rationale behind the success of multi-model ensembles in seasonal forecasting-II calibration and combination. Tellus A. 2005;57:234-252. doi:10.1111/j.1600-0870.2005.00104.x
#'@references Van Schaeybroeck, B., & Vannitsem, S. (2011). Post-processing through linear regression. Nonlinear Processes in Geophysics, 18(2), 147. doi:10.5194/npg-18-147-2011
#'@references Van Schaeybroeck, B., & Vannitsem, S. (2015). Ensemble post-processing using member-by-member approaches: theoretical aspects. Quarterly Journal of the Royal Meteorological Society, 141(688), 807-818.  doi:10.1002/qj.2397
#'
#'@param exp an object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the seasonal forecast experiment data in the element named \code{$data}.
#'@param obs an object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the observed data in the element named \code{$data}.
#'@param cal.method is the calibration method used, can be either \code{bias}, \code{evmos}, \code{mse_min} or \code{crps_min}. Default value is \code{mse_min}.
#'@param eval.method is the sampling method used, can be either \code{in-sample} or \code{leave-one-out}. Default value is the \code{leave-one-out} cross validation.
#'@param multi.model is a boolean that is used only for the \code{mse_min} method. If multi-model ensembles or ensembles of different sizes are used, it must be set to \code{TRUE}. By default it is \code{FALSE}. Differences between the two approaches are generally small but may become large when using small ensemble sizes. Using multi.model when the calibration method is \code{bias}, \code{evmos} or \code{crps_min} will not affect the result.
#'@return an object of class \code{s2dv_cube} containing the calibrated forecasts in the element \code{$data} with the same dimensions of the experimental data.
#'
#'@import s2dverification
#'@import abind
#'
#'@seealso \code{\link{CST_Load}}
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
#'a <- CST_Calibration(exp = exp, obs = obs, cal.method = "mse_min", eval.method = "in-sample")
#'str(a)
#'@export


CST_Calibration <- function(exp, obs, cal.method = "mse_min", eval.method = "leave-one-out",  multi.model = F) {
  if (!inherits(exp, "s2dv_cube") || !inherits(obs, "s2dv_cube")) {
    stop("Parameter 'exp' and 'obs' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load.")
  }

  if (dim(obs$data)["member"] != 1) {
    stop("The length of the dimension 'member' in the component 'data' ",
         "of the parameter 'obs' must be equal to 1.")
  }
  
  if(!missing(multi.model) & !(cal.method == "mse_min")){
	  warning(paste0("The multi.model parameter is ignored when using the calibration method ", cal.method))
  }
  exp$data <- .calibration.wrap(exp = exp$data, obs = obs$data,
    cal.method = cal.method, eval.method = eval.method,  
    multi.model =  multi.model)
  
  exp$Datasets <- c(exp$Datasets, obs$Datasets)
  exp$source_files <- c(exp$source_files, obs$source_files)
  return(exp)
  
}

.calibration.wrap <- function(exp, obs, cal.method, eval.method, multi.model) {
	
  dim.exp <- dim(exp)
  amt.dims.exp <- length(dim.exp)
  dim.obs <- dim(obs)
  amt.dims.obs <- length(dim.obs)
  dim.names.exp <- names(dim.exp)
  dim.names.obs <- names(dim.obs)
  
  target.dim.names.exp <- c("member", "sdate")
  target.dim.names.obs <- c("sdate")
  
  if (!all(target.dim.names.exp %in% dim.names.exp)) {
    stop("Parameter 'exp' must have the dimensions 'member' and 'sdate'.")
  }

  if (!all(c("sdate") %in% dim.names.obs)) {
    stop("Parameter 'obs' must have the dimension 'sdate'.")
  }

  if (any(is.na(exp)))  {
    warning("Parameter 'exp' contains NA values.")
  }
  
  if (any(is.na(obs))) {
    warning("Parameter 'obs' contains NA values.")
  }
  
  if (dim(obs)['member']!=1){
    stop("Parameter 'obs' must have a member dimension with length=1")
  }
  
  obs <- Subset(obs, along = "member", indices = 1, drop = "selected")
    
  calibrated <- Apply(data = list(mod = exp, obs = obs),
    cal.method = cal.method,
    eval.method = eval.method,
    multi.model = multi.model,
    target_dims = list(mod = target.dim.names.exp, obs = target.dim.names.obs),
    fun = .cal)$output1

  dexes <- match(names(dim(exp)), names(dim(calibrated)))
  calibrated <- aperm(calibrated, dexes)
  dimnames(calibrated) <- dimnames(exp)[dexes]
  

  return(calibrated)
}

.make.eval.train.dexes <- function(eval.method, amt.points){ 
  if(eval.method == "leave-one-out"){
    dexes.lst <- lapply(seq(1, amt.points), function(x) return(list(eval.dexes = x, train.dexes = seq(1, amt.points)[-x])))
  } else if (eval.method == "in-sample"){
    dexes.lst <- list(list(eval.dexes = seq(1, amt.points), train.dexes = seq(1, amt.points)))
  } else {
    stop(paste0("unknown sampling method: ",eval.method))
  }
  return(dexes.lst)
}

.cal <- function(mod, obs, cal.method, eval.method, multi.model) {

  var.obs <- as.vector(obs)
  var.fc <- mod
  dims.fc <- dim(var.fc)
  amt.mbr <- dims.fc[1]
  amt.sdate <- dims.fc[2]
  var.cor.fc <- NA * var.fc
  
  eval.train.dexeses <- .make.eval.train.dexes(eval.method, amt.points = amt.sdate)
  amt.resamples <- length(eval.train.dexeses)
  for (i.sample in seq(1, amt.resamples)) {
    # defining training (tr) and evaluation (ev) subsets 
    eval.dexes <- eval.train.dexeses[[i.sample]]$eval.dexes
    train.dexes <- eval.train.dexeses[[i.sample]]$train.dexes
    
    fc.ev <- var.fc[ , eval.dexes, drop = FALSE]
    fc.tr <- var.fc[ , train.dexes]
    obs.tr <- var.obs[train.dexes , drop = FALSE] 
    
    if(cal.method == "bias"){
	  var.cor.fc[ , eval.dexes] <- fc.ev + mean(obs.tr, na.rm = TRUE) - mean(fc.tr, na.rm = TRUE)
	} else if(cal.method == "evmos"){
	  #calculate ensemble and observational characteristics
	  quant.obs.fc.tr <- .calc.obs.fc.quant(obs = obs.tr, fc = fc.tr)
      #calculate value for regression parameters
      init.par <- c(.calc.evmos.par(quant.obs.fc.tr))
	  #correct evaluation subset
      var.cor.fc[ , eval.dexes] <- .correct.evmos.fc(fc.ev , init.par)
	} else if (cal.method == "mse_min"){
	  #calculate ensemble and observational characteristics
	  quant.obs.fc.tr <- .calc.obs.fc.quant(obs = obs.tr, fc = fc.tr)
      #calculate value for regression parameters
      init.par <- .calc.mse.min.par(quant.obs.fc.tr, multi.model)
	  #correct evaluation subset
      var.cor.fc[ , eval.dexes] <- .correct.mse.min.fc(fc.ev , init.par)      
    } else if (cal.method == "crps_min"){
	  #calculate ensemble and observational characteristics
	  quant.obs.fc.tr <- .calc.obs.fc.quant.ext(obs = obs.tr, fc = fc.tr)
      #calculate initial value for regression parameters
      init.par <- c(.calc.mse.min.par(quant.obs.fc.tr), 0.001)
      init.par[3] <- sqrt(init.par[3])
      #calculate regression parameters on training dataset
      optim.tmp <- optim(par = init.par, fn = .calc.crps.opt, 
        quant.obs.fc = quant.obs.fc.tr, method = "Nelder-Mead")
      
      mbm.par <- optim.tmp$par
	  #correct evaluation subset
      var.cor.fc[ , eval.dexes] <- .correct.crps.min.fc(fc.ev , mbm.par)
    } else {
	  stop("unknown calibration method: ",cal.method)
    }
  }
  names(dim(var.cor.fc)) <- c("member", "sdate")
  return(var.cor.fc) 
}

.calc.obs.fc.quant <- function(obs, fc){ #function to calculate different quantities of a series of ensemble forecasts and corresponding observations
  amt.mbr <- dim(fc)[1]
  obs.per.ens <- InsertDim(var = obs, posdim = 1, lendim = amt.mbr)
  fc.ens.av <- apply(fc, c(2), mean, na.rm = TRUE)
  cor.obs.fc <- cor(fc.ens.av, obs, use = "complete.obs")
  obs.av <- mean(obs, na.rm = TRUE)
  obs.sd <- sd(obs, na.rm = TRUE)
  return(
    append(
      .calc.fc.quant(fc = fc),
      list(
        obs.per.ens = obs.per.ens,
        cor.obs.fc = cor.obs.fc,
        obs.av = obs.av,
        obs.sd = obs.sd
      )
    )
  )
}

.calc.obs.fc.quant.ext <- function(obs, fc){ #extended function to calculate different quantities of a series of ensemble forecasts and corresponding observations
  amt.mbr <- dim(fc)[1]
  obs.per.ens <- InsertDim(var = obs, posdim = 1, lendim = amt.mbr)
  fc.ens.av <- apply(fc, c(2), mean, na.rm = TRUE)
  cor.obs.fc <- cor(fc.ens.av, obs, use = "complete.obs")
  obs.av <- mean(obs, na.rm = TRUE)
  obs.sd <- sd(obs, na.rm = TRUE)

  return(
    append(
      .calc.fc.quant.ext(fc = fc),
      list(
        obs.per.ens = obs.per.ens,
        cor.obs.fc = cor.obs.fc,
        obs.av = obs.av,
        obs.sd = obs.sd
      )
    )
  )
}


.calc.fc.quant <- function(fc){ #function to calculate different quantities of a series of ensemble forecasts
  amt.mbr <- dim(fc)[1]
  fc.ens.av <- apply(fc, c(2), mean, na.rm = TRUE)
  fc.ens.av.av <- mean(fc.ens.av, na.rm = TRUE)
  fc.ens.av.sd <- sd(fc.ens.av, na.rm = TRUE)
  fc.ens.av.per.ens <- InsertDim(var = fc.ens.av, posdim = 1, lendim = amt.mbr)
  fc.ens.sd <- apply(fc, c(2), sd, na.rm = TRUE)
  fc.ens.var.av.sqrt <- sqrt(mean(fc.ens.sd^2, na.rm = TRUE))
  fc.dev <- fc - fc.ens.av.per.ens
  fc.dev.sd <- sd(fc.dev, na.rm = TRUE)
  fc.av <- mean(fc, na.rm = TRUE)
  fc.sd <- sd(fc, na.rm = TRUE)
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

.calc.fc.quant.ext <- function(fc){ #extended function to calculate different quantities of a series of ensemble forecasts

  amt.mbr <- dim(fc)[1]
  repmat1.tmp <- InsertDim(var = fc, posdim = 1, lendim = amt.mbr)
  repmat2.tmp <- aperm(repmat1.tmp, c(2, 1, 3))
  spr.abs <- apply(abs(repmat1.tmp - repmat2.tmp), c(3), mean, na.rm = TRUE)
  spr.abs.per.ens <- InsertDim(var = spr.abs, posdim = 1, lendim = amt.mbr)

  return(
    append(.calc.fc.quant(fc),
	  list(spr.abs = spr.abs, spr.abs.per.ens = spr.abs.per.ens))
  )
}

#Below are the core or elementary functions to calculate the regression parameters for the different methods
.calc.mse.min.par <- function(quant.obs.fc, multi.model = F){
  par.out <- rep(NA, 3)
  
  if(multi.model){
    par.out[3] <- with(quant.obs.fc, obs.sd * sqrt(1. - cor.obs.fc^2) / fc.ens.var.av.sqrt)
  } else {
    par.out[3] <- with(quant.obs.fc, obs.sd * sqrt(1. - cor.obs.fc^2) / fc.dev.sd)
  }
  par.out[2] <- with(quant.obs.fc, cor.obs.fc * obs.sd / fc.ens.av.sd)
  par.out[1] <- with(quant.obs.fc, obs.av - par.out[2] * fc.ens.av.av, na.rm = TRUE)
  
  return(par.out)
}
.calc.evmos.par <- function(quant.obs.fc){
  par.out <- rep(NA, 2)
  par.out[2] <- with(quant.obs.fc, obs.sd / fc.sd)
  par.out[1] <- with(quant.obs.fc, obs.av - par.out[2] * fc.ens.av.av, na.rm = TRUE)
  return(par.out)
}
#Below are the core or elementary functions to calculate the functions necessary for the minimization of crps
.calc.crps.opt <- function(par, quant.obs.fc){
  return( 
    with(quant.obs.fc, 
      mean(abs(obs.per.ens - (par[1] + par[2] * fc.ens.av.per.ens +
	    ((par[3])^2 + par[4] / spr.abs.per.ens) * fc.dev)), na.rm = TRUE) -
        mean(abs((par[3])^2 * spr.abs + par[4]) / 2., na.rm = TRUE)
    )
  )
}

#Below are the core or elementary functions to correct the evaluation set based on the regression parameters
.correct.evmos.fc <- function(fc, par){
  quant.fc.mp <- .calc.fc.quant(fc = fc)
  return(with(quant.fc.mp, par[1] + par[2] * fc))
}
.correct.mse.min.fc <- function(fc, par){
  quant.fc.mp <- .calc.fc.quant(fc = fc)
  return(with(quant.fc.mp, par[1] + par[2] * fc.ens.av.per.ens + fc.dev * par[3]))
}
.correct.crps.min.fc <- function(fc, par){
  quant.fc.mp <- .calc.fc.quant.ext(fc = fc)
  return(with(quant.fc.mp, par[1] + par[2] * fc.ens.av.per.ens + fc.dev * abs((par[3])^2 + par[4] / spr.abs)))
}
