#'Multivariate Root Mean Square Error (RMSE) 
#'
#'@author Deborah Verfaillie, \email{deborah.verfaillie@bsc.es}
#'@description This function calculates the RMSE from multiple variables, as the mean of each variable's RMSE scaled by its observed standard deviation. Variables can be weighted based on their relative importance (defined by the user).
#'
#'@param exp a list of objects, one for each variable, of class \code{s2dv_cube} as returned by \code{CST_Anomaly} function, containing the anomaly of the seasonal forecast experiment data in the element named \code{$data}.
#'@param obs a list of objects, one for each variable (in the same order than the input in 'exp') of class \code{s2dv_cube} as returned by \code{CST_Anomaly} function, containing the observed anomaly data in the element named \code{$data}.
#'
#'@param weight (optional) a vector of weight values to assign to each variable. If no weights are defined, a value of 1 is assigned to every variable.
#'
#'@return an object of class \code{s2dv_cube} containing the RMSE in the element \code{$data} which is an array with two datset dimensions equal to the 'dataset' dimension in the \code{exp$data} and \code{obs$data} inputs.  An array with dimensions: c(number of exp, number of obs, 1 (the multivariate RMSE value), number of lat, number of lon)
#'
#'@seealso \code{\link[s2dverification]{RMS}} and \code{\link{CST_Load}}
#'@import s2dverification
#'@examples
#'# Creation of sample s2dverification objects. These are not complete
#'# s2dverification objects though. The Load function returns complete objects.
#'# using package zeallot is optional:
#' library(zeallot)
#'# Example with 2 variables
#'mod1 <- 1 : (1 * 3 * 4 * 5 * 6 * 7)
#'mod2 <- 1 : (1 * 3 * 4 * 5 * 6 * 7)
#'dim(mod1) <- c(dataset = 1, member = 3, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'dim(mod2) <- c(dataset = 1, member = 3, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'obs1 <- 1 : (1 * 1 * 4 * 5 * 6 * 7)
#'obs2 <- 1 : (1 * 1 * 4 * 5 * 6 * 7)
#'dim(obs1) <- c(dataset = 1, member = 1, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'dim(obs2) <- c(dataset = 1, member = 1, sdate = 4, ftime = 5, lat = 6, lon = 7)
#'lon <- seq(0, 30, 5)
#'lat <- seq(0, 25, 5)
#'exp1 <- list(data = mod1, lat = lat, lon = lon, Datasets = "EXP1", 
#'             source_files = "file1", Variable = list('pre'))
#'attr(exp1, 'class') <- 's2dv_cube'
#'exp2 <- list(data = mod2, lat = lat, lon = lon, Datasets = "EXP2", 
#'             source_files = "file2", Variable = list('tas'))
#'attr(exp2, 'class') <- 's2dv_cube'
#'obs1 <- list(data = obs1, lat = lat, lon = lon, Datasets = "OBS1", 
#'             source_files = "file1", Variable = list('pre'))
#'attr(obs1, 'class') <- 's2dv_cube'
#'obs2 <- list(data = obs2, lat = lat, lon = lon, Datasets = "OBS2", 
#'             source_files = "file2", Variable = list('tas'))
#'attr(obs2, 'class') <- 's2dv_cube'
#'
#'c(ano_exp1, ano_obs1) %<-% CST_Anomaly(exp1, obs1, cross = TRUE, memb = TRUE)
#'c(ano_exp2, ano_obs2) %<-% CST_Anomaly(exp2, obs2, cross = TRUE, memb = TRUE)
#'ano_exp <- list(exp1, exp2)
#'ano_obs <- list(ano_obs1, ano_obs2)
#'weight <- c(1, 2)
#'a <- CST_MultivarRMSE(exp = ano_exp, obs = ano_obs, weight = weight)
#'str(a)
#'@export
CST_MultivarRMSE <- function(exp, obs, weight = NULL) {
  if (!is.list(exp) | !is.list(obs)) {
    stop("Parameters 'exp' and 'obs' must be lists of 's2dv_cube' objects")
  }

  if (!(all(sapply(exp, inherits, 's2dv_cube')))) {
     stop("Elements of the list in parameter 'exp' must be of the class ",
          "'s2dv_cube', as output by CSTools::CST_Load.")
  }

  if (!(all(sapply(obs, inherits, 's2dv_cube')))) {
     stop("Elements of the list in parameter 'obs' must be of the class ",
          "'s2dv_cube', as output by CSTools::CST_Load.")
  }

  if (length(exp) != length(obs)) {
    stop("Parameters 'exp' and 'obs' must be of the same length.")
  }
  
  nvar <- length(exp)
  
  if (nvar < 2) {
    stop("Parameters 'exp' and 'obs'  must contain at least two", 
         " s2dverification objects for two different variables.")
  }
  
  for (j in 1 : nvar) {
    if (!is.null(names(dim(exp[[j]]$data))) & !is.null(names(dim(obs[[j]]$data)))) {
      if (all(names(dim(exp[[j]]$data)) %in% names(dim(obs[[j]]$data)))) {
        dimnames <- names(dim(exp[[j]]$data))
      } else {
        stop("Dimension names of element 'data' from parameters 'exp'",
             " and 'obs' should have the same name dimmension.")
      }
    } else {
       stop("Element 'data' from parameters 'exp' and 'obs'",
           " should have dimmension names.")
    }
  }
  
  if (is.null(weight)) {
    weight <- c(rep(1, nvar))
  } else if (length(weight) != nvar) {
    stop("Parameter 'weight' must have a length equal to the number ",
         "of variables.")
  }
  obs_var <- unlist(lapply(obs, function(x) {
                    x[[which(names(x) == 'Variable')]]}))

  exp_var <- unlist(lapply(exp, function(x) {
                    x[[which(names(x) == 'Variable')]]}))
  
  if (all(exp_var != obs_var)) {
    stop("Variables in parameters 'exp' and 'obs' must be in the same order.")
  }
  mvrmse <- 0
  sumweights <- 0
  for (j in 1 : nvar) {
    # seasonal average of anomalies
    AvgExp <- MeanListDim(exp[[j]]$data, narm = TRUE, c(2, 4))
    AvgObs <- MeanListDim(obs[[j]]$data, narm = TRUE, c(2, 4))
    # multivariate RMSE (weighted) 
    rmse <- RMS(AvgExp, AvgObs, posloop = 1, posRMS = 2, conf = FALSE)
    stdev <- sd(AvgObs)
    mvrmse <- mvrmse + (rmse / stdev * as.numeric(weight[j]))
    sumweights <- sumweights + as.numeric(weight[j])
  }
  mvrmse <- mvrmse / sumweights 
  
  names(dim(mvrmse)) <- c(dimnames[1], dimnames[1], 'statistics', dimnames[5 : 6])
  exp_Datasets <- unlist(lapply(exp, function(x) {
                         x[[which(names(x) == 'Datasets')]]}))
  exp_source_files <- unlist(lapply(exp, function(x) {
                         x[[which(names(x) == 'source_files')]]}))
  obs_Datasets <- unlist(lapply(obs, function(x) {
                         x[[which(names(x) == 'Datasets')]]}))
  obs_source_files <- unlist(lapply(obs, function(x) {
                         x[[which(names(x) == 'source_files')]]}))

  exp <- exp[[1]]
  exp$data <- mvrmse
  exp$Datasets <- c(exp_Datasets, obs_Datasets)
  exp$source_files <- c(exp_source_files, obs_source_files)
  exp$Variable <- c(exp_var)
  return(exp)
}
