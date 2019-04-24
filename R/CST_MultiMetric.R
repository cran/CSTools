#'Multiple Metrics applied in Multiple Model Anomalies
#'
#'@author Mishra Niti, \email{niti.mishra@bsc.es}
#'@author Perez-Zanon Nuria, \email{nuria.perez@bsc.es}
#'@description This function calculates correlation (Anomaly Correlation Coefficient; ACC), root mean square error (RMS) and the root mean square error skill score (RMSSS) of individual anomaly models and multi-models mean (if desired) with the observations.
#'
#'@param exp an object of class \code{s2dv_cube} as returned by \code{CST_Anomaly} function, containing the anomaly of the seasonal forecast experiment data in the element named \code{$data}.
#'@param obs an object of class \code{s2dv_cube} as returned by \code{CST_Anomaly} function, containing the anomaly of observed data in the element named \code{$data}.
#'@param metric a character string giving the metric for computing the maximum skill. This must be one of the strings 'correlation', 'rms' or 'rmsss.
#'@param multimodel a logical value indicating whether a Multi-Model Mean should be computed.
#'
#'@return an object of class \code{s2dv_cube} containing the statistics of the selected metric in the element \code{$data} which is an array with two datset dimensions equal to the 'dataset' dimension in the \code{exp$data} and \code{obs$data} inputs. If \code{multimodel} is TRUE, the greatest first dimension correspons to the Multi-Model Mean. The third dimension contains the statistics selected. For metric \code{correlation}, the third dimension is of length four and they corresponds to the lower limit of the 95\% confidence interval, the statistics itselfs, the upper limit of the 95\% confidence interval and the 95\% significance level. For metric \code{rms}, the third dimension is length three and they corresponds to the lower limit of the 95\% confidence interval, the RMSE and the upper limit of the 95\% confidence interval. For metric \code{rmsss}, the third dimension is length two and they corresponds to the statistics itselfs and the p-value of the one-sided Fisher test with Ho: RMSSS = 0.
#'@seealso \code{\link[s2dverification]{Corr}}, \code{\link[s2dverification]{RMS}}, \code{\link[s2dverification]{RMSSS}} and \code{\link{CST_Load}}
#'@references 
#'Mishra, N., Prodhomme, C., & Guemas, V. (n.d.). Multi-Model Skill Assessment of Seasonal Temperature and Precipitation Forecasts over Europe, 29–31.\url{http://link.springer.com/10.1007/s00382-018-4404-z}
#' 
#'@import s2dverification
#'@import stats
#'@examples
#'library(zeallot)
#'mod <- 1 : (2 * 3 * 4 * 5 * 6 * 7)
#'dim(mod) <- c(dataset = 2, member = 3, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'obs <- 1 : (1 * 1 * 4 * 5 * 6 * 7)
#'dim(obs) <- c(dataset = 1, member = 1, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'lon <- seq(0, 30, 5)
#'lat <- seq(0, 25, 5)
#'exp <- list(data = mod, lat = lat, lon = lon)
#'obs <- list(data = obs, lat = lat, lon = lon)
#'attr(exp, 'class') <- 's2dv_cube'
#'attr(obs, 'class') <- 's2dv_cube'
#'c(ano_exp, ano_obs) %<-% CST_Anomaly(exp = exp, obs = obs, cross = TRUE, memb = TRUE)
#'a <- CST_MultiMetric(exp = ano_exp, obs = ano_obs)
#'str(a)
#'@export
CST_MultiMetric <- function(exp, obs, metric = "correlation", multimodel = TRUE) {
  if (!inherits(exp, 's2dv_cube') || !inherits(obs, 's2dv_cube')) {
    stop("Parameter 'exp' and 'obs' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load.")
  }

  if (dim(obs$data)['member'] != 1) {
    stop("The length of the dimension 'member' in the component 'data' ",
         "of the parameter 'obs' must be equal to 1.")
  }

  if (!is.null(names(dim(exp$data))) & !is.null(names(dim(obs$data)))) {
    if (all(names(dim(exp$data)) %in% names(dim(obs$data)))) {
      dimnames <- names(dim(exp$data))
    } else {
      stop("Dimension names of element 'data' from parameters 'exp'",
           " and 'obs' should have the same name dimmension.")
    }
  } else {
    stop("Element 'data' from parameters 'exp' and 'obs'",
         " should have dimmension names.")
  }
  
  if (!is.logical(multimodel)) {
    stop("Parameter 'multimodel' must be a logical value.")
  }
  if (length(multimodel) > 1) {
    multimodel <- multimodel[1]
    warning("Parameter 'multimodel' has length > 1 and only the first ",
            "element will be used.")
  }
  if (length(metric) > 1) {
    metric <- metric[1]
    warning("Parameter 'multimodel' has length > 1 and only the first ",
            "element will be used.")
  }
  
  # seasonal average of anomalies per model
  AvgExp <- MeanListDim(exp$data, narm = T, c(2, 4))
  AvgObs <- MeanListDim(obs$data, narm = T, c(2, 4))
  
  # indv model correlation
  if (metric == 'correlation') {
    corr <- Corr(AvgExp, AvgObs, posloop = 1, poscor = 2)
  } else if (metric == 'rms') {
    corr <- RMS(AvgExp, AvgObs, posloop = 1, posRMS = 2)
  } else if (metric == 'rmsss') {
    corr <- RMSSS(AvgExp, AvgObs, posloop = 1, posRMS = 2)
  } else {
    stop("Parameter 'metric' must be a character string indicating ",
         "one of the options: 'correlation', 'rms' or 'rmse'.")
  }
  if (multimodel == TRUE) {
    # seasonal avg of anomalies for multi-model
    AvgExp_MMM <- MeanListDim(AvgExp, narm = TRUE, 1)
    AvgObs_MMM <- MeanListDim(AvgObs, narm = TRUE, 1)
    # multi model correlation
    if (metric == 'correlation') {
      corr_MMM <- Corr(var_exp = InsertDim(AvgExp_MMM, 1, 1),
                       var_obs = InsertDim(AvgObs_MMM, 1, 1),
                       posloop = 1, poscor = 2)
    } else if (metric == 'rms') {
      corr_MMM <- RMS(var_exp = InsertDim(AvgExp_MMM, 1, 1),
                      var_obs = InsertDim(AvgObs_MMM, 1, 1),
                      posloop = 1, posRMS = 2)
    } else if (metric == 'rmsss') {
      corr_MMM <- RMSSS(var_exp = InsertDim(AvgExp_MMM, 1, 1),
                        var_obs = InsertDim(AvgObs_MMM, 1, 1), 
                        posloop = 1, posRMS = 2)
    } 
    corr <- abind::abind(corr, corr_MMM, along = 1)
  }
  
  names(dim(corr)) <- c(dimnames[1], dimnames[1], 'statistics', dimnames[5 : 6])
  
  #exp$data <- ano$ano_exp
  #obs$data <- ano$ano_obs
  exp$data <- corr
  exp$Datasets <- c(exp$Datasets, obs$Datasets)
  exp$source_files <- c(exp$source_files, obs$source_files)
  return(exp)
}
