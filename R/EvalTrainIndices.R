#' Generate Training and Evaluation Indices for Cross-Validation
#'
#' This function generates training and evaluation indices based on different
#' cross-validation methods.
#' @author Theertha Kariyathan, \email{theertha.kariyathan@bsc.es}
#' 
#' @param eval.method Character. The cross-validation method. Options include:
#'   -\code{leave-k-out}: Leaves out \code{k} points at a time for evaluation.
#'   -\code{retrospective}: Uses past data for training and a future point for evaluation.
#'   -\code{in-sample}: Uses the entire dataset for both training and evaluation.
#'   -\code{hindcast-vs-forecast}: Uses all years from the hindcast sdate dimension
#'    as training and all years from the forecast sdate dimension will be corrected.
#'   
#' @param sample.length Integer. Length of the sample (in years).
#' @param sample.length_cor Integer. Length of forecast sample (in years) in
#'  \code{hindcast-vs-forecast} method.
#' @param k Positive integer. Default = 1.
#' In method \code{leave-k-out}, \code{k} is expected to be odd integer, 
#' indicating the number of points to leave out.
#' In method \code{retrospective}, \code{k} can be any positive integer, 
#' indicating when to start. 
#' @param tail.out Logical for method \code{leave-k-out}. Default = TRUE.  
#'  TRUE to remove both extremes keeping the same sample size for all k-folds 
#'  (e.g. sample.length=50, k=3, eval.dexes=1, train.dexes=(3,49)). 
#'  FALSE to remove only the corresponding tail 
#'  (e.g. sample.length=50, k=3, eval.dexes=1, train.dexes=(3,50))
#' @return A list of lists, where each element contains:
#'   \code{eval.dexes}: Indices of evaluation points.
#'   \code{train.dexes}: Indices of training points.
#'
#' @examples
#' # Leave-k-out cross-validation
#' EvalTrainIndices("leave-k-out", sample.length = 10, sample.length_cor = 5, k = 3)
#' # Retrospective cross-validation
#' EvalTrainIndices("retrospective", sample.length = 10, sample.length_cor = 5, k = 3)
#' # In-sample validation
#' EvalTrainIndices("in-sample", sample.length = 10, sample.length_cor = 5)
#' # Hindcast vs. Forecast validation
#' EvalTrainIndices("hindcast-vs-forecast", sample.length = 10, sample.length_cor = 5)
#'
#' @export
EvalTrainIndices <- function(eval.method, sample.length, sample.length_cor, 
                           k = 1, tail.out = TRUE) {
  
  if (eval.method == "leave-k-out") {
    
    if (is.null(k) || k >= sample.length || k < 1 || k %% 2 == 0 ) { 
      warning("k is expected to be a positive odd integer less than the sample.length")
    }
    
    if (k == 1) {
      tail.out == TRUE
    }
    amt.seq <- 1:sample.length
    strike = (k - 1) / 2
    alt_amt.seq <- c(tail(amt.seq, strike), amt.seq, head(amt.seq, strike))
    
    dexes.lst <- lapply(seq(1, sample.length), function(x, kfold = k) {
      if (tail.out) {
        ind = alt_amt.seq[x:(x + (kfold - 1))]
      } else{
        ind = ((x - strike):(x + strike))
        ind = ind[ind > 0 & ind <= sample.length]
      }
      return(list(eval.dexes = x,  train.dexes = amt.seq[!(amt.seq %in% ind)]))
    })
  } else if (eval.method == "retrospective") {
    # k can be any integer indicating the when to start
    if(k >= sample.length){
      warning("k is expected to be less than the sample.length")
    }
    dexes.lst <- base::Filter(length, lapply(seq(1, sample.length),
                                             function(x, mindata = k) {
                                               if (x > k) { 
                                                 eval.dexes <- x
                                                 train.dexes <- 1:(x-1)
                                                 return(list(eval.dexes = x,
                                                             train.dexes = 1:(x-1)))
                                               }})) 

  } else if (eval.method == "in-sample") {
    dexes.lst <- list(list(eval.dexes = seq(1, sample.length), 
                           train.dexes = seq(1, sample.length)))
  } else if (eval.method == "hindcast-vs-forecast") {
    dexes.lst <- list(list(eval.dexes = seq(1, sample.length_cor),
                           train.dexes = seq(1, sample.length)))
  } else {
    stop(paste0("unknown sampling method: ", eval.method))
  }
  
  return(dexes.lst)
} 


.make.eval.train.dexes <- function(eval.method, amt.points, amt.points_cor) {
  warning("This function is deprecated. Please use 'EvalTrainIndices()' instead.")
  if (eval.method == "leave-one-out") {
    dexes.lst <- lapply(seq(1, amt.points), function(x) return(list(eval.dexes = x,
                        train.dexes = seq(1, amt.points)[-x])))
  } else if (eval.method == "in-sample") {
    dexes.lst <- list(list(eval.dexes = seq(1, amt.points),
                           train.dexes = seq(1, amt.points)))
  } else if (eval.method == "hindcast-vs-forecast") {
    dexes.lst <- list(list(eval.dexes = seq(1,amt.points_cor),
                           train.dexes = seq(1, amt.points)))
  } else {
    stop(paste0("unknown sampling method: ", eval.method))
  }
  return(dexes.lst)
}
