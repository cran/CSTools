#'Quantile Mapping for seasonal or decadal forecast data
#'
#'@description This function is a wrapper of fitQmap and doQmap from package
#''qmap' to be applied on the object of class 's2dv_cube'. The quantile mapping
#'adjustment between an experiment, typically a hindcast, and observation is 
#'applied to the experiment itself or to a provided forecast.

#'@author Nuria Perez-Zanon, \email{nuria.perez@bsc.es}
#'@param exp An object of class \code{s2dv_cube}.
#'@param obs An object of class \code{s2dv_cube}.
#'@param exp_cor A multidimensional array with named dimensions in which the
#'  quantile mapping correction should be applied. If it is not specified, the
#'  correction is applied to object 'exp' using leave-one-out cross-validation.
#'  This is useful to correct a forecast when the hindcast is provided in parameter 'exp'.
#'@param sdate_dim A character string indicating the dimension name in which 
#'  cross-validation would be applied when exp_cor is not provided. 'sdate' by 
#'  default.
#'@param memb_dim A character string indicating the dimension name where
#'  ensemble members are stored in the experimental arrays. It can be NULL if 
#'  there is no ensemble member dimension. It is set as 'member' by default.
#'@param window_dim  A character string indicating the dimension name in which extra
#'  samples are stored. This dimension is joined to the 'member' dimension.
#'  This is useful to correct daily data, for which robust statistics can be obtained
#'  by creating a window of dates around the target date.
#'@param method A character string indicating the method to be used:'PTF', 
#'  'DIST', 'RQUANT', 'QUANT', 'SSPLIN'. By default, the empirical quantile 
#'  mapping 'QUANT' is used.
#'@param na.rm A logical value indicating if missing values should be removed   
#'  (FALSE by default).
#'@param ncores An integer indicating the number of cores for parallel 
#'  computation using multiApply function. The default value is NULL (1). 
#'@param ... Additional parameters to be used by the method choosen. See qmap 
#'  package for details.
#'@param eval.method A character string indicating the evaluation method for cross-validaton.
#'the default method is 'leave-k-out', other available methods are 
#''retrospective', 'in-sample', 'hindcast-vs-forecast'.
#'@param k Positive integer. Default = 1.
#' In method \code{leave-k-out}, \code{k} is expected to be odd integer, 
#' indicating the number of points to leave out.
#' In method \code{retrospective}, \code{k} can be any positive integer greater than 1, 
#' indicating when to start.
#'
#'@return An object of class \code{s2dv_cube} containing the experimental data
#'after applying the quantile mapping correction.
#'
#'@seealso \code{\link[qmap]{fitQmap}} and \code{\link[qmap]{doQmap}} 
#'@examples
#'# Use synthetic data
#'exp <- NULL
#'exp$data <- 1 : c(1 * 3 * 5 * 4 * 3 * 2)
#'dim(exp$data) <- c(dataset = 1, member = 3, sdate = 5, ftime = 4,
#'                   lat = 3, lon = 2)
#'class(exp) <- 's2dv_cube'
#'obs <- NULL
#'obs$data <- 101 : c(100 + 1 * 1 * 5 * 4 * 3 * 2)
#'dim(obs$data) <- c(dataset = 1, member = 1, sdate = 5, ftime = 4,
#'                   lat = 3, lon = 2)
#'class(obs) <- 's2dv_cube'
#'res <- CST_QuantileMapping(exp, obs)
#'
#'@import qmap 
#'@import multiApply 
#'@import s2dv
#'@export
CST_QuantileMapping <- function(exp, obs, exp_cor = NULL, sdate_dim = 'sdate',
                                memb_dim = 'member', window_dim = NULL, 
                                method = 'QUANT', na.rm = FALSE, 
                                eval.method = "leave-k-out", k = 1, 
                                ncores = NULL, ...) {
  # Check 's2dv_cube'
  if (!inherits(exp, 's2dv_cube') || !inherits(obs, 's2dv_cube')) {
    stop("Parameter 'exp' and 'obs' must be of the class 's2dv_cube'.")
  }
  if (!is.null(exp_cor)) {
    if (!inherits(exp_cor, 's2dv_cube')) {
      stop("Parameter 'exp_cor' must be of the class 's2dv_cube'.")
    }
  }

  QMapped <- QuantileMapping(exp = exp$data, obs = obs$data, 
                             exp_cor = exp_cor$data,
                             sdate_dim = sdate_dim, memb_dim = memb_dim,
                             window_dim = window_dim, method = method,
                             eval.method = eval.method, k = k,
                             na.rm = na.rm, ncores = ncores, ...)
  if (is.null(exp_cor)) {
    exp$data <- QMapped
    exp$attrs$Datasets <- c(exp$attrs$Datasets, obs$attrs$Datasets)
    exp$attrs$source_files <- c(exp$attrs$source_files, obs$attrs$source_files)
    return(exp)
  } else {
    exp_cor$data <- QMapped
    exp_cor$attrs$Datasets <- c(exp_cor$attrs$Datasets, exp$attrs$Datasets, 
                                obs$attrs$Datasets)
    exp_cor$attrs$source_files <- c(exp_cor$attrs$source_files, exp$attrs$source_files, 
                                    obs$attrs$source_files)
    return(exp_cor)
  }
}

#'Quantile Mapping for seasonal or decadal forecast data
#'
#'@description This function is a wrapper of fitQmap and doQmap from package
#''qmap' to be applied on multi-dimensional arrays. The quantile mapping 
#'adjustment between an experiment, typically a hindcast, and observation is 
#'applied to the experiment itself or to a provided forecast.
#'
#'@author Nuria Perez-Zanon, \email{nuria.perez@bsc.es}
#'@param exp A multidimensional array with named dimensions containing the 
#'  hindcast. 
#'@param obs A multidimensional array with named dimensions containing the 
#'  reference dataset.
#'@param exp_cor A multidimensional array with named dimensions in which the
#'  quantile mapping correction should be applied. If it is not specified, the
#'  correction is applied to object 'exp' using leave-one-out cross-validation.
#'  This is useful to correct a forecast when the hindcast is provided in parameter 'exp'.
#'@param sdate_dim A character string indicating the dimension name in which
#'  cross-validation would be applied when exp_cor is not provided. 'sdate' by
#'  default.
#'@param memb_dim A character string indicating the dimension name where
#'  ensemble members are stored in the experimental arrays. It can be NULL if 
#'  there is no ensemble member dimension. It is set as 'member' by default.
#'@param window_dim A character string indicating the dimension name in which extra
#'  samples are stored. This dimension is joined to the 'member' dimension.
#'  This is useful to correct daily data, for which robust statistics can be obtained
#'  by creating a window of dates around the target date.
#'@param method A character string indicating the method to be used: 'PTF',
#'  'DIST', 'RQUANT', 'QUANT', 'SSPLIN'. By default, the empirical quantile 
#'  mapping 'QUANT' is used. 
#'@param na.rm A logical value indicating if missing values should be removed  
#'  (FALSE by default). 
#'@param ncores An integer indicating the number of cores for parallel 
#'  computation using multiApply function. The default value is NULL (1). 
#'@param ... Additional parameters to be used by the method choosen. See qmap 
#'  package for details.
#'@param eval.method A character string indicating the evaluation method for cross-validaton.
#'the default method is 'leave-k-out', other available methods are 
#''retrospective', 'in-sample', 'hindcast-vs-forecast'.
#'@param k Positive integer. Default = 1.
#' In method 'leave-k-out', 'k' is expected to be odd integer, 
#' indicating the number of points to leave out.
#' In method 'retrospective', 'k' can be any positive integer greater than 1, 
#' indicating when to start.
#'
#'@return An array containing the experimental data after applying the quantile
#'mapping correction.
#' 
#'@seealso \code{\link[qmap]{fitQmap}} and \code{\link[qmap]{doQmap}} 
#'@examples
#'# Use synthetic data
#'set.seed(123)
#'exp <- as.numeric(1:prod(6,10,15))
#'dim(exp) <- c(member = 6, syear = 10, window = 15)
#'obs <- as.numeric(rnorm(prod(1,10,15), 50))
#'dim(obs) <- c(member = 1, syear = 10, window = 15)
#'fcst <- 100*(1:prod(8,1,1))
#'dim(fcst) <- c(member = 8, syear = 1, swindow = 1)
#'res <- QuantileMapping(exp = exp, obs = obs, exp_cor = fcst, 
#'                       memb_dim = 'member', sdate_dim = 'syear', window_dim = 'window')
#'@import qmap 
#'@import multiApply 
#'@import s2dv
#'@export
QuantileMapping <- function(exp, obs, exp_cor = NULL, sdate_dim = 'sdate',
                            memb_dim = 'member', window_dim = NULL, 
                            method = 'QUANT', na.rm = FALSE, 
                            eval.method = "leave-k-out", k = 1, 
                            ncores = NULL, ...) {
  # exp and obs
  obsdims <- names(dim(obs))
  expdims <- names(dim(exp))
  if (!is.array(exp) || !is.numeric(exp)) {
    stop("Parameter 'exp' must be a numeric array.")
  }
  if (!is.array(obs) || !is.numeric(obs)) {
    stop("Parameter 'obs' must be a numeric array.")
  }
  if (is.null(expdims)) {
    stop("Parameter 'exp' must have dimension names.")
  }
  if (is.null(obsdims)) {
    stop("Parameter 'obs' must have dimension names.")
  }
  # sdate_dim
  if (!is.character(sdate_dim) | length(sdate_dim) != 1) {
    stop("Parameter 'sdate_dim' must be a character string.")
  }
  if (!sdate_dim %in% expdims | !sdate_dim %in% obsdims) {
    stop("Parameter 'sdate_dim' is not found in 'exp' or 'obs' dimension.")
  }
  if (dim(exp)[sdate_dim] == 1 || dim(obs)[sdate_dim] == 1) {
    stop("Parameter 'exp' and 'obs' must have dimension length of 'sdate_dim' bigger than 1.")
  }
  # exp_cor
  if (!is.null(exp_cor)) {
    if (is.null(names(dim(exp_cor)))) {
      stop("Parameter 'exp_cor' must have dimension names.")
    }
    if (!sdate_dim %in% names(dim(exp_cor))) {
      stop("Parameter 'sdate_dim' is not found in 'exp_cor' dimension.")
    }
  }
  # method
  if (!(method %in% c('PTF', 'DIST', 'RQUANT', 'QUANT', 'SSPLIN')) | length(method) != 1) {
    stop("Parameter 'method' must be one of the following methods: ",
         "'PTF', 'DIST', 'RQUANT', 'QUANT', 'SSPLIN'.")
  }
  # memb_dim
  if (is.null(memb_dim)) {
    remove_member <- TRUE
    memb_dim <- "temp_memb_dim"
    exp <- InsertDim(exp, posdim = 1, lendim = 1, name = "temp_memb_dim")
    obs <- InsertDim(obs, posdim = 1, lendim = 1, name = "temp_memb_dim")
    obsdims <- names(dim(obs))
    expdims <- names(dim(exp))
    if (!is.null(exp_cor)) {
      exp_cor <- InsertDim(exp_cor, posdim = 1, lendim = 1, name = "temp_memb_dim")
    }
  } else {
    remove_member <- FALSE
    if (!all(memb_dim %in% obsdims)) {
      obs <- InsertDim(obs, posdim = 1, lendim = 1,
                       name = memb_dim[!(memb_dim %in% obsdims)])
      obsdims <- names(dim(obs))
    }
    if (any(!memb_dim %in% expdims)) {
      stop(paste0("Parameter 'memb_dim' is not found in 'exp' dimensions. ", 
                  "Set it as NULL if there is no member dimension."))
    }
  }
  sample_dims <- c(memb_dim, sdate_dim)
  # window_dim
  if (!is.null(window_dim)) {
    if (!(window_dim %in% obsdims)) {
      stop("Parameter 'window_dim' is not found in 'obs'.")
    }
    obs <- CSTools::MergeDims(obs, c(memb_dim, window_dim))
    if (window_dim %in% expdims) {
      exp <- CSTools::MergeDims(exp, c(memb_dim, window_dim))
      warning("Parameter 'window_dim' is found in exp and is merged to 'memb_dim'.")
    }
  }
  # na.rm
  if (!is.logical(na.rm) | length(na.rm) > 1) {
    stop("Parameter 'na.rm' must be one logical value.")
  }
  # ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) | ncores %% 1 != 0 | ncores <= 0 |
      length(ncores) > 1) {
      stop("Parameter 'ncores' must be either NULL or a positive integer.")
    }
  }

  ###############################
  if (!is.null(exp_cor)) {
    qmaped <- Apply(list(exp, obs, exp_cor), target_dims = sample_dims, 
                    fun = .qmapcor, method = method, sdate_dim = sdate_dim,
                    na.rm = na.rm, ..., 
                    ncores = ncores)$output1
  } else {
    qmaped <- Apply(list(exp, obs), target_dims = sample_dims,
                    fun = .qmapcor, exp_cor = NULL, method = method,
                    sdate_dim = sdate_dim, na.rm = na.rm, ...,                
                    ncores = ncores)$output1
  }
  # remove added 'temp_memb_dim'
  if (remove_member) {
    dim(qmaped) <- dim(qmaped)[-which(names(dim(qmaped)) == "temp_memb_dim")]
  }

  return(qmaped)
}

# Example
# exp <- as.numeric(1:prod(6,10,15))
# dim(exp) <- c(member = 6, syear = 10, window = 15)
# obs <- as.numeric(rnorm(prod(1,10,15), 50))
# dim(obs) <- c(member = 1, syear = 10, window = 15)
# fcst <- 100*(1:prod(8,1,1))
# dim(fcst) <- c(member = 8, syear = 1, swindow = 1)
# 
# 
# res_sample <- QuantileMapping(exp = exp, obs = obs, exp_cor = NULL, 
#                               memb_dim = 'member', sdate_dim = 'syear', 
#                               window_dim = 'window', eval.method = "in-sample")  
# 
# res_leave <- QuantileMapping(exp = exp, obs = obs, exp_cor = NULL, 
#                              memb_dim = 'member', sdate_dim = 'syear', 
#                              window_dim = 'window', eval.method = "leave-k-out", k = 1)  
# 
# res_retro <- QuantileMapping(exp = exp, obs = obs, exp_cor = NULL, 
#                              memb_dim = 'member', sdate_dim = 'syear', 
#                              window_dim = 'window', eval.method = "retrospective", 
#                              k = 3) 
# 
# res_hind <- QuantileMapping(exp = exp, obs = obs, exp_cor = fcst, 
#                             memb_dim = 'member', sdate_dim = 'syear', 
#                             window_dim = 'window', eval.method = "hindcast-vs-forecast") 
#' 

.qmapcor <- function(exp, obs, exp_cor = NULL, sdate_dim = 'sdate', 
                     eval.method = 'leave-k-out', k = 1,
                     method = 'QUANT', na.rm = FALSE, ...) {
  
  
  if(!is.null(exp_cor) & eval.method != "hindcast-vs-forcast"){
    eval.method = 'hindcast-vs-forecast'
   # warning("Defaulting to 'hindcast-vs-forcast' as exp_cor is not 'NULL'")
  }
  
  if(eval.method == "retrospective" & k == 1){
    stop("k = 1 not expected, trainindex need at least two non-NA values to interpolate ")
  } #can replace with tryCatch for smooth fail
  
  
  index <- EvalTrainIndices(eval.method = eval.method, sample.length = dim(exp)[sdate_dim], k = k,
                            sample.length_cor = dim(exp_cor)[sdate_dim])
  applied <- exp * NA
  
  for(x in index){
  
    
    nas_pos <- which(!is.na(exp[, x$eval.dexes]))
    obs2 <- as.vector(obs[, x$train.dexes]) 
    exp2 <- as.vector(exp[, x$train.dexes])
    exp_cor2 <- as.vector(exp[, x$eval.dexes])
    # remove NAs
    obs2 <- obs2[!is.na(obs2)]
    exp2 <- exp2[!is.na(exp2)]   
    exp_cor2 <- exp_cor2[!is.na(exp_cor2)]
    
    adjust <- fitQmap(obs2, exp2, method = 'QUANT', ...)
    if(is.null(exp_cor)){
      if(na.rm){
        if(eval.method == "in-sample"){
          applied[nas_pos] <- doQmap(exp_cor2, adjust, ...) 
        }else{
          applied[nas_pos, x$eval.dexes] <- doQmap(exp_cor2, adjust, ...)
        }
      }else{
        
        if (anyNA(obs[, x$train.dexes]) | anyNA(exp[, x$train.dexes])) { 
          applied[, x$eval.dexes] <- NA
        }else{
          if(eval.method == "in-sample"){
            applied[nas_pos] <- doQmap(exp_cor2, adjust, ...)
          }else{
            applied[nas_pos, x$eval.dexes] <- doQmap(exp_cor2, adjust, ...)
          }
        }
      }
    }else{
      applied <- exp_cor * NA
      # add check if method is hindcast else give warning that since exp_cor is not null hindcast is selected
      #hindcast vs forecast
      if (na.rm) {
        tryCatch({
          
          adjust <- fitQmap(obs2, exp2,
                            method = method, ...)
          
          applied[!is.na(exp_cor)] <- doQmap(exp_cor[!is.na(exp_cor)],
                                             adjust, ...)
         
        },
        error = function(error_message) {
         return(applied)
        })
      } else { 
        adjust <- fitQmap(as.vector(obs), as.vector(exp), method = method, ...)
        applied <- doQmap(as.vector(exp_cor), adjust, ...)
      }
      dim(applied) <- dim(exp_cor)
    }
  }
  applied
}




