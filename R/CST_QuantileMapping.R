#'Quantiles Mapping for seasonal or decadal forecast data
#'
#'@description This function is a wrapper from fitQmap and doQmap from package 'qmap'to be applied in CSTools objects of class 's2dv_cube'. The quantile mapping adjustment between an experiment, tipically a hindcast, and observations is applied to the experiment itself or to a provided forecast.
#'
#'@author Nuria Perez-Zanon, \email{nuria.perez@bsc.es}
#'@param exp an object of class \code{s2dv_cube}
#'@param obs an object of class \code{s2dv_cube}
#'@param exp_cor an object of class \code{s2dv_cube} in which the quantile mapping correction will be applied. If it is not specified, the correction is applied in object \code{exp}.
#'@param sample_dims a character vector indicating the dimensions that can be used as sample for the same distribution
#'@param sample_length a numeric value indicating the length of the timeseries window to be used as sample for the sample distribution and correction. By default, NULL, the total length of the timeseries will be used. 
#'@param method a character string indicating the method to be used: 'PTF','DIST','RQUANT','QUANT','SSPLIN'. By default, the empirical quantile mapping 'QUANT' is used.
#'@param ncores an integer indicating the number of parallel processes to spawn for the use for parallel computation in multiple cores.
#'@param ... additional arguments passed to the method specified by \code{method}.
#'
#'@details The different methods are:
#'\itemize{ 
#'\item{'PTF'} {fits a parametric transformations to the quantile-quantile relation of observed and modelled values. See \code{?qmap::fitQmapPTF}.}
#' \item{'DIST'} {fits a theoretical distribution to observed and to modelled time series. See \code{?qmap::fitQmapDIST}.}
#'\item{'RQUANT'} {estimates the values of the quantile-quantile relation of observed and modelled time series for regularly spaced quantiles using local linear least square regression. See \code{?qmap::fitQmapRQUANT}.}
#'\item{'QUANT'} {estimates values of the empirical cumulative distribution function of observed and modelled time series for regularly spaced quantiles. See \code{?qmap::fitQmapQUANT}.}
#'\item{'SSPLIN'} {fits a smoothing spline to the quantile-quantile plot of observed and modelled time series. See \code{?qmap::fitQmapSSPLIN}.}}
#'All methods accepts some common arguments:
#'\itemize{
#'\item{wet.day} {logical indicating whether to perform wet day correction or not.(Not available in 'DIS' method)}
#'\item{qstep} {NULL or a numeric value between 0 and 1.}}
#' When providing a forecast to be corrected through the pararmeter \code{exp_cor}, some inputs might need to be modified. The quantile correction is compute by comparing objects passed through 'exp' and 'obs' parameters, this correction will be later applied to the forecast provided in 'exp_cor'. Imaging the case of 'exp' and 'obs' having several start dates, stored using a dimension e.g. 'sdate', 'sample_dims' include this dimension 'sdate' and 'exp_cor' has forecasts for several sdates but different from the ones in 'exp'. In this case, the correction computed with 'exp' and 'obs' would be applied for each 'sdate' of 'exp_cor' separately. This example corresponds to a case of split a dataset in training set and validation set.
#' 
#'@return an oject of class \code{s2dv_cube} containing the experimental data after applyingthe quantile mapping correction. 
#') <- c(dataset = 1, member = 10, sdate = 20, ftime = 60 , 
#'@import qmap
#'@import multiApply
#'@import abind
#'
#'@seealso \code{qmap::fitQmap} and \code{qmap::doQmap}
#'@examples
#'library(qmap)
#'exp <- 1 : (1 * 5 * 10 * 6 * 2 * 3)
#'dim(exp) <- c(dataset = 1, member = 10, sdate = 5, ftime = 6 , 
#'              lat = 2, lon = 3)
#'exp <- list(data = exp)
#'class(exp) <- 's2dv_cube'
#'obs <- 101 : (100 + 1 * 1 * 5 * 6 * 2 * 3)
#'dim(obs) <- c(dataset = 1, member = 1, sdate = 5, ftime = 6 ,
#'              lat = 2, lon = 3)
#'obs <- list(data = obs)
#'class(obs) <- 's2dv_cube'
#'res <- CST_QuantileMapping(exp, obs, method = 'RQUANT')
#'\donttest{
#'exp <- lonlat_data$exp
#'obs <- lonlat_data$obs
#'res <- CST_QuantileMapping(exp, obs)
#'
#'exp_cor <- exp
#'exp_cor$data <- exp_cor$data[,,1,,,]
#'dim(exp_cor$data) <- c(dataset = 1, member = 15, sdate = 1, ftime = 3, 
#'                       lat = 22, lon = 53)
#'res <- CST_QuantileMapping(exp, obs, exp_cor,
#'                           sample_dims = c('sdate', 'ftime', 'member'))
#'res <- CST_QuantileMapping(exp, obs, exp_cor,
#'                           sample_dims = c('ftime', 'member'))
#'data(obsprecip)
#'data(modprecip)
#'exp <- modprecip$MOSS[1:10000]
#'dim(exp) <- c(time = length(exp))
#'exp <- list(data = exp)
#'class(exp) <- 's2dv_cube'
#'obs <- obsprecip$MOSS[1:10000]
#'dim(obs) <- c(time = length(obs))
#'obs <- list(data = obs)
#'class(obs) <- 's2dv_cube'
#'res <- CST_QuantileMapping(exp = exp, obs = obs, sample_dims = 'time',
#'                           method = 'DIST')
#'# Example using different lenght of members and sdates:
#'exp <- lonlat_data$exp
#'exp$data <- exp$data[,,1:4,,,]
#'dim(exp$data) <- c(dataset = 1, member = 15, sdate = 4, ftime = 3, 
#'                       lat = 22, lon = 53)
#'obs <- lonlat_data$obs
#'obs$data <- obs$data[,,1:4, ,,]
#'dim(obs$data) <- c(dataset = 1, member = 1, sdate = 4, ftime = 3, 
#'                       lat = 22, lon = 53)
#'exp_cor <- lonlat_data$exp
#'exp_cor$data <- exp_cor$data[,1:5,5:6,,,]
#'dim(exp_cor$data) <- c(dataset = 1, member = 5, sdate = 2, ftime = 3, 
#'                      lat = 22, lon = 53)
#'res <- CST_QuantileMapping(exp, obs, exp_cor,
#'                           sample_dims = c('sdate', 'ftime', 'member'))
#'exp_cor <- lonlat_data$exp
#'exp_cor$data <- exp_cor$data[,,5:6,,,]
#'dim(exp_cor$data) <- c(dataset = 1, member = 15, sdate = 2, ftime = 3, 
#'                       lat = 22, lon = 53)
#'res <- CST_QuantileMapping(exp, obs, exp_cor,
#'                           sample_dims = c('sdate', 'ftime', 'member'))
#'}
#'@export
CST_QuantileMapping <- function(exp, obs, exp_cor = NULL,
                                sample_dims = c('sdate', 'ftime', 'member'), 
                                sample_length = NULL, method = 'QUANT', ncores = NULL, ...) {
    if (!inherits(exp, 's2dv_cube') || !inherits(obs, 's2dv_cube')) {
        stop("Parameter 'exp' and 'obs' must be of the class 's2dv_cube', ",
             "as output by CSTools::CST_Load.")
    }
    if (!is.null(exp_cor)) {
        if (!inherits(exp_cor, 's2dv_cube')) {
            stop("Parameter 'exp_cor' must be of the class 's2dv_cube', ",
                 "as output by CSTools::CST_Load.")
        }
    }
    if (!(method %in% c('PTF','DIST','RQUANT','QUANT','SSPLIN'))) {
        stop("Parameter 'method' must be one of the following methods: ",
             "'PTF','DIST','RQUANT','QUANT','SSPLIN'.")
    }
    QMapped <- QuantileMapping(exp = exp$data, obs = obs$data, exp_cor = exp_cor$data,
                               sample_dims = sample_dims, sample_length = sample_length,
                               method = method, ncores = ncores, ...)
    if (is.null(exp_cor)) {
        exp$data <- QMapped
        exp$source_files <- c(exp$source_files, obs$source_files)
    } else {
        exp_cor$data <- QMapped
        exp_cor$source_files <- c(exp$source_files, obs$source_files, exp_cor$source_files)
        exp <- exp_cor
    }
    return(exp)
}
#'Quantiles Mapping for seasonal or decadal forecast data
#'
#'@description This function is a wrapper from fitQmap and doQmap from package 'qmap'to be applied in CSTools objects of class 's2dv_cube'. The quantile mapping adjustment between an experiment, tipically a hindcast, and observations is applied to the experiment itself or to a provided forecast.
#'
#'@author Nuria Perez-Zanon, \email{nuria.perez@bsc.es}
#'@param exp a multi-dimensional array with named dimensions containing the hindcast.
#'@param obs a multi-dimensional array with named dimensions (the same as the provided in 'exp') containing the reference dataset.
#'@param exp_cor a multi-dimensional array with named dimensions in which the quantile mapping correction will be applied. If it is not specified, the correction is applied in object \code{exp}.
#'@param sample_dims a character vector indicating the dimensions that can be used as sample for the same distribution
#'@param sample_length a numeric value indicating the length of the timeseries window to be used as sample for the sample distribution and correction. By default, NULL, the total length of the timeseries will be used. 
#'@param method a character string indicating the method to be used: 'PTF','DIST','RQUANT','QUANT','SSPLIN'. By default, the empirical quantile mapping 'QUANT' is used.
#'@param ncores an integer indicating the number of parallel processes to spawn for the use for parallel computation in multiple cores.
#'@param ... additional arguments passed to the method specified by \code{method}.
#'
#'@details The different methods are:
#'\itemize{ 
#'\item{'PTF'} {fits a parametric transformations to the quantile-quantile relation of observed and modelled values. See \code{?qmap::fitQmapPTF}.}
#' \item{'DIST'} {fits a theoretical distribution to observed and to modelled time series. See \code{?qmap::fitQmapDIST}.}
#'\item{'RQUANT'} {estimates the values of the quantile-quantile relation of observed and modelled time series for regularly spaced quantiles using local linear least square regression. See \code{?qmap::fitQmapRQUANT}.}
#'\item{'QUANT'} {estimates values of the empirical cumulative distribution function of observed and modelled time series for regularly spaced quantiles. See \code{?qmap::fitQmapQUANT}.}
#'\item{'SSPLIN'} {fits a smoothing spline to the quantile-quantile plot of observed and modelled time series. See \code{?qmap::fitQmapSSPLIN}.}}
#'All methods accepts some common arguments:
#'\itemize{
#'\item{wet.day} {logical indicating whether to perform wet day correction or not.(Not available in 'DIS' method)}
#'\item{qstep} {NULL or a numeric value between 0 and 1.}}
#'@return an oject of class \code{s2dv_cube} containing the experimental data after applyingthe quantile mapping correction. 
#') <- c(dataset = 1, member = 10, sdate = 20, ftime = 60 , 
#'@import qmap
#'@import multiApply
#'@import abind
#'
#'@seealso \code{qmap::fitQmap} and \code{qmap::doQmap}
#'@export
QuantileMapping <- function(exp, obs, exp_cor = NULL, sample_dims = 'ftime', 
                            sample_length = NULL, method = 'QUANT', ncores = NULL, ...) {
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
  if (any(is.na(obs))) {
    warning("Parameter 'obs' contains NA values.")
  }
  if (!is.null(exp_cor)) {
    exp_cordims <- names(dim(exp_cor))
    if (is.null(exp_cordims)) {
      stop("Parameter 'exp_cor' must have dimension names.")
    }
  }
  if (!all(sample_dims %in% expdims)) {
    stop("Parameter 'sample_dims' must be a vector of string character ",
         "indicating names of exiting dimension in parameter 'exp'.")
  }
  ## The sample_dims could be of any length (even different between exp and obs)
  ## but the repeated dims that aren't in sample_dims should be drop:
  commondims <- obsdims[obsdims %in% expdims]
  commondims <- names(which(unlist(lapply(commondims, function(x) {
                             dim(obs)[obsdims == x] != dim(exp)[expdims == x]}))))
  if (any(!(commondims %in% sample_dims))) {
      todrop <- commondims[!(commondims %in% sample_dims)]
      todrop <- match(todrop, obsdims)
        if (all(dim(obs)[todrop] != 1)) {
            stop("Review parameter 'sample_dims' or the data dimensions",
                 "since multiple dimensions with different length have ",
                 "being found in the data inputs that don't match with ",
                 "'sample_dims' parameter.")
         } else {
             obs <- adrop(obs, drop = todrop)
         }
  }
  if (!all(sample_dims %in% obsdims)) {
    newobsdims <- sample_dims[!sample_dims %in% obsdims]
    dim(obs) <- c(dim(obs), 1 : length(newobsdims))
    names(dim(obs))[-c(1:length(obsdims))] <- newobsdims
  }

  if (!is.null(exp_cor)) {
    commondims <- exp_cordims[exp_cordims %in% expdims]
    commondims <- names(which(unlist(lapply(commondims, function(x) {
                             dim(exp_cor)[exp_cordims == x] != dim(exp)[expdims == x]}))))
    if (any(commondims %in% sample_dims)) {
      todrop <- commondims[(commondims %in% sample_dims)]
      todroppos <- match(todrop, sample_dims)
      if (all(dim(exp_cor)[todrop] != 1)) {
        warning(paste("The sample_dims", paste(todrop, collapse = " "), 
                      "are not used when applying the",
                      "correction to 'exp_cor'"))
        sample_dims <- list(sample_dims, sample_dims, sample_dims[-todroppos]) 
      } else {
        exp_cor <- adrop(exp_cor, drop = todroppos)
      }
    } else {
      todrop <- commondims[!(commondims %in% sample_dims)]
      todrop <- match(todrop, obsdims)
      if (all(dim(exp_cor)[todrop] != 1)) {
        stop("Review parameter 'sample_dims' or the data dimensions ",
             "since multiple dimensions with different length have ",
             "being found in the data inputs that don't match with ",
             "'sample_dims' parameter.")
      } else {
        exp_cor <- adrop(exp_cor, drop = todrop)
      }
    }
  } 

  if (!is.null(sample_length) & !is.numeric(sample_length)) {
    warning("Parameter 'sample_length' has not been correctly defined and ",
            "the whole length of the timeseries will be used.")
    sample_length <- NULL
  }
  if (length(sample_length) > 1) {
    warning("Parameter 'sample_length' has length > 1 and only the first ",
            "element will be used.")
    sample_length <- sample_length[1]
  }
  if (!is.character(method)) {
    warning("Parameter 'method' must be a character string indicating ",
         "one of the following methods: 'PTF', 'DIST', 'RQUANT', 
         'QUANT', 'SSPLIN'. Method 'QUANT' is being used.")
    method = 'QUANT'
  }
  if (length(method) > 1) {
    warning("Parameter 'method' has length > 1 and only the first element",
            " is used.")
    method <- method[1]
  }

#  qmaped <- Apply(list(exp, obs), target_dims = sample_dims, fun = qmapcor, ...,  
#                  exp_cor = exp_cor, sample_length = sample_length,
  if (is.null(exp_cor)) {
    qmaped <- Apply(list(exp, obs), target_dims = sample_dims,
                  fun = qmapcor, ..., sample_length = sample_length,
                  method = method, ncores = ncores)$output1
  } else {
    qmaped <- Apply(list(exp, obs, exp_cor), target_dims = sample_dims, 
                    fun = qmapcor, ..., sample_length = sample_length,
                    method = method, ncores = ncores)$output1
  }
  pos <- match(expdims, names(dim(qmaped)))
  out_names <- names(dim(exp))
    if (length(pos) < length(dim(qmaped))) {
      toadd <- length(dim(qmaped)) - length(pos)
      toadd <- seq(max(pos) + 1, max(pos) + toadd, 1)
      pos <- c(pos, toadd)
      new <- names(dim(qmaped))[names(dim(qmaped)) %in% out_names == FALSE]
      out_names <- c(out_names, new)
    } 
  qmaped <- aperm(qmaped, pos)
  names(dim(qmaped)) <- out_names
  return(qmaped)
}
qmapcor <- function(exp, obs, exp_cor = NULL, sample_length = NULL, method = 'QUANT',
                    ...) {
    dimnames_exp <- names(dim(exp))
    dimnames_obs <- names(dim(obs))
    if (length(dimnames_exp) != length(dimnames_obs)) {
        stop("Parameters 'exp' and 'obs' must have the same number of dimensions.")
    }
    if (!all(dimnames_exp %in% dimnames_obs)) {
        stop("Parameters 'exp' and 'obs' must have the same dimension names.")
    } 
    if (!(any(names(dim(exp)) %in% 'ftime') | any(names(dim(exp)) %in% 'time') | 
          any(names(dim(exp)) %in% 'sdate'))) {
        stop("Parameters 'exp' and 'obs' must have a temporal dimension named ",
             "'time', 'ftime' or 'sdate'.")
    }
        
    dimensions <- dim(exp)
    if (!is.null(exp_cor)) {
        dimensions <- dim(exp_cor)
    }        
    if (any(names(dim(exp)) %in% 'ftime') | any(names(dim(exp)) %in% 'time')) {
        time_dim <- which(names(dim(exp)) == 'ftime' | names(dim(exp)) %in% 'time')
    } else {
        time_dim <- which(names(dim(exp)) == 'sdate')
    }
    if (any(names(dim(obs)) %in% 'ftime') | any(names(dim(obs)) %in% 'time')) {
        time_dim_obs <- which(names(dim(obs)) == 'ftime' | names(dim(obs)) %in% 'time')
    } else {
        time_dim_obs <- which(names(dim(obs)) == 'sdate')
    }
    if (is.null(sample_length)) {
        sample_length <- dim(exp)[time_dim]
    }
    nsamples <- dim(exp)[time_dim]/sample_length
    if (nsamples %% 1 != 0) { 
        # add NA to complete the last sample
        nsamples <- ceiling(nsamples)
        fillsample1D <- rep(NA, nsamples * sample_length -  dim(exp)[time_dim])
        if (length(dim(exp)) > 1) {
            fillsample <- rep(fillsample1D, prod(dim(exp)[-time_dim]))
            dims <- dim(exp)
            exp <- c(exp, fillsample)
            dim(exp) <- c(dims[-time_dim], sample = sample_length, 
                          ceiling(dims[time_dim]/sample_length))
            fillsample <- rep(fillsample1D, prod(dim(obs)[-time_dim_obs]))
            dims <- dim(obs)
            obs <- c(obs, fillsample)
            dim(obs) <- c(dims[-time_dim_obs], sample = sample_length,
                          ceiling(dims[time_dim_obs]/sample_length))
         } else {
            exp <- abind(exp, fillsample1D, along = time_dim)
            names(dim(exp)) <- dimnames_exp
            obs <- abind(obs, fillsample1D, along = time_dim_obs)
            names(dim(obs)) <- dimnames_obs
            dim(exp) <- c(dim(exp)[-time_dim], sample = sample_length,
                          dim(exp)[time_dim]/sample_length)
            dim(obs) <- c(dim(obs)[-time_dim_obs], sample = sample_length,
                          dim(obs)[time_dim_obs]/sample_length)
         }
    } else {
            dim(exp) <- c(dim(exp)[-time_dim], sample = sample_length,
                          dim(exp)[time_dim]/sample_length)
            dim(obs) <- c(dim(obs)[-time_dim_obs], sample = sample_length,
                          dim(obs)[time_dim_obs]/sample_length)          
    }
    if (any(names(dim(exp)) %in% 'ftime') | any(names(dim(exp)) %in% 'time')) {
        new_time_dim_exp <- which(names(dim(exp)) == 'ftime' | names(dim(exp)) %in% 'time')
    } else {
        new_time_dim_exp <- which(names(dim(exp)) == 'sdate')
    }
    if (any(names(dim(obs)) %in% 'ftime') | any(names(dim(obs)) %in% 'time')) {
        new_time_dim_obs <- which(names(dim(obs)) == 'ftime' | names(dim(obs)) %in% 'time')
    } else {
        new_time_dim_obs <- which(names(dim(obs)) == 'sdate')
    }

    if (!is.null(exp_cor)) {
        if (any(names(dim(exp_cor)) %in% 'ftime') | any(names(dim(exp_cor)) %in% 'time')) {
            time_dim_cor <- which(names(dim(exp_cor)) == 'ftime' | names(dim(exp_cor)) %in% 'time')
        } else {
            time_dim_cor <- which(names(dim(exp_cor)) == 'sdate')
        }

        nsamples <- dimensions[time_dim_cor]/sample_length
        if (nsamples %% 1 != 0) {
        nsamples <- ceiling(nsamples)
        fillsample1D <- rep(NA, nsamples * sample_length -  dimensions[time_dim_cor])
            if (length(dimensions) > 1) {
                fillsample <- rep(fillsample1D, prod(dimensions[-time_dim_cor]))
                 exp_cor <- c(exp_cor, fillsample)
                 dim(exp_cor) <- c(dim(exp_cor)[-time_dim_cor], sample = sample_length,
                               ceiling(dim(exp_cor)[time_dim_cor]/sample_length))
             } else {
                exp_cor <- abind(exp_cor, fillsample1D, along = time_dim_cor)
                names(dim(exp_cor)) <- names(dimensions)
             }
        }  
        dim(exp_cor) <- c(dim(exp_cor)[-time_dim_cor], sample = sample_length,
                          dim(exp_cor)[time_dim_cor]/sample_length)
        if (any(names(dim(exp_cor)) %in% 'ftime') | any(names(dim(exp_cor)) %in% 'time')) {
            new_time_dim_cor <- which(names(dim(exp_cor)) == 'ftime' |
                                      names(dim(exp_cor)) %in% 'time')
        } else {
            new_time_dim_cor <- which(names(dim(exp_cor)) == 'sdate')
        }

    } else {
        time_dim_cor <- time_dim
        exp_cor <- exp
        new_time_dim_cor  <- new_time_dim_exp
    }
    applied <- NULL
    for (i in 1 : nsamples) {
        if (i <= dim(obs)[new_time_dim_obs]) {
            sample_obs <- as.vector(asub(obs, idx = i, dims = new_time_dim_obs))
            sample_exp <- as.vector(asub(exp, idx = i, dims = new_time_dim_exp))
        } else {
            sample_obs <- as.vector(asub(obs, idx = dim(obs)[new_time_dim_obs],
                                    dims = new_time_dim_obs))
            sample_exp <- as.vector(asub(exp, idx = dim(exp)[new_time_dim_exp],
                                    dims = new_time_dim_exp))
        }
        if (i >= dim(obs)[new_time_dim_obs]) {
            sample_obs <- sample_obs[!is.na(sample_obs)]
            sample_exp <- sample_exp[!is.na(sample_exp)]
        }
        if (sum(sample_obs) == 0) {
            warning("The total sum of observed data sample in the sample number ",
                    i, ", is zero and the function may crash.")
        }
        if (sum(sample_exp) == 0) {
            warning("The total sum of experimental data sample in the sample number ",
                    i, ", is zero and the function may crash.")
        }
        if (length(sample_exp) < sample_length) {
            warning("The length of the sample used, ", length(sample_exp),
                    ", in the sample number ", i,
                    ", is smaller than the defined in parameter 'sample_length'.")
        }
        adjust <- fitQmap(sample_obs, sample_exp, method = method,
                          ...)
        sample_cor <- as.vector(asub(exp_cor, idx = i, dims = new_time_dim_cor))
        if (i == nsamples) {
            sample_cor <- sample_cor[!is.na(sample_cor)]
        }
        applied <- c(applied, doQmap(x = sample_cor, fobj = adjust, ...))
    }
    if (any(is.na(exp_cor))) {
        pos <- which(!is.na(exp_cor))
        exp_cor[pos] <- applied
     applied <- exp_cor
    }
    dim(applied) <- dimensions
    return(applied) 
}
