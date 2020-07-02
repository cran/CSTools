#'Anomalies relative to a climatology along selected dimension with or without cross-validation
#'
#'@author Perez-Zanon Nuria, \email{nuria.perez@bsc.es}
#'@author Pena Jesus, \email{jesus.pena@bsc.es}
#'@description This function computes the anomalies relative to a climatology computed along the 
#'selected dimension (usually starting dates or forecast time) allowing the application or not of
#'crossvalidated climatologies. The computation is carried out independently for experimental and 
#'observational data products.
#'
#'@param exp an object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the seasonal forecast experiment data in the element named \code{$data}.
#'@param obs an object of class \code{s2dv_cube} as returned by \code{CST_Load} function, containing the observed data in the element named \code{$data}.'
#'@param cross A logical value indicating whether cross-validation should be applied or not. Default = FALSE.
#'@param memb A logical value indicating whether Clim() computes one climatology for each experimental data 
#'product member(TRUE) or it computes one sole climatology for all members (FALSE). Default = TRUE.
#'@param filter_span a numeric value indicating the degree of smoothing. This option is only available if parameter \code{cross} is set to FALSE.  
#'@param dim_anom An integer indicating the dimension along which the climatology will be computed. It 
#'usually corresponds to 3 (sdates) or 4 (ftime). Default = 3.
#'
#' @return A list with two S3 objects, 'exp' and 'obs', of the class 's2dv_cube', containing experimental and date-corresponding observational anomalies, respectively. These 's2dv_cube's can be ingested by other functions in CSTools.
#'
#'@import s2dverification
#'
#'@examples
#'# Example 1:
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
#'
#'anom1 <- CST_Anomaly(exp = exp, obs = obs, cross = FALSE, memb = TRUE)
#'str(anom1)
#'anom2 <- CST_Anomaly(exp = exp, obs = obs, cross = TRUE, memb = TRUE)
#'str(anom2)
#'
#'anom3 <- CST_Anomaly(exp = exp, obs = obs, cross = TRUE, memb = FALSE)
#'str(anom3)
#'
#'anom4 <- CST_Anomaly(exp = exp, obs = obs, cross = FALSE, memb = FALSE)
#'str(anom4)
#'
#'anom5 <- CST_Anomaly(lonlat_data$exp)
#'
#'anom6 <- CST_Anomaly(obs = lonlat_data$obs)
#'
#'@seealso \code{\link[s2dverification]{Ano_CrossValid}}, \code{\link[s2dverification]{Clim}} and \code{\link{CST_Load}}
#'
#'
#'@export
CST_Anomaly <- function(exp = NULL, obs = NULL, cross = FALSE, memb = TRUE,
                        filter_span = NULL, dim_anom = 3) {
 
  if (!inherits(exp, 's2dv_cube') & !is.null(exp) || 
      !inherits(obs, 's2dv_cube') & !is.null(obs)) {
    stop("Parameter 'exp' and 'obs' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load.")
  }

 if (!is.null(obs)) {
    if (dim(obs$data)['member'] != 1) {
    stop("The length of the dimension 'member' in the component 'data' ",
         "of the parameter 'obs' must be equal to 1.")
    }
  }
  case_exp = case_obs = 0
  if (is.null(exp) & is.null(obs)) {
    stop("One of the parameter 'exp' or 'obs' cannot be NULL.")
  }
  if (is.null(exp)) {
    exp <- obs 
    case_obs = 1
    warning("Parameter 'exp' is not provided and will be recycled.")
  }  
  if (is.null(obs)) {
    obs <- exp 
    case_exp = 1
    warning("Parameter 'obs' is not provided and will be recycled.")
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
  dim_exp <- dim(exp$data)
  dim_obs <- dim(obs$data)

  dimnames_data <- names(dim_exp)
  if (dim_exp[dim_anom] == 1 | dim_obs[dim_anom] == 1) {
    stop("The length of dimension 'dim_anom' in label 'data' of the parameter ",
         "'exp' and 'obs' must be greater than 1.")
  }
  if (!any(names(dim_exp)[dim_anom] == c('sdate', 'time', 'ftime'))) {
    warning("Parameter 'dim_anom' correspond to a position name different ",
            "than 'sdate', 'time' or 'ftime'.")
  } 
  if (!is.logical(cross) | !is.logical(memb) ) {
    stop("Parameters 'cross' and 'memb' must be logical.")
  }
  if (length(cross) > 1 | length(memb) > 1 ) {
    cross <- cross[1]
    warning("Parameter 'cross' has length greater than 1 and only the first element",
            "will be used.")
  }
  if (length(memb) > 1) {
     memb <- memb[1]
     warning("Parameter 'memb' has length greater than 1 and only the first element",
             "will be used.")
  }
  
  # Selecting time dimension through dimensions permutation
  if (dim_anom != 3) {
    dimperm <- 1 : length(dim_exp)
    dimperm[3] <- dim_anom
    dimperm[dim_anom] <- 3
    
    var_exp <- aperm(exp$data, perm = dimperm)
    var_obs <- aperm(obs$data, perm = dimperm)
    
    #Updating permuted dimensions
    dim_exp <- dim(exp$data)
    dim_obs <- dim(obs$data)
  }
  
  
  # Computating anomalies
  #----------------------
  
  # With cross-validation
  if (cross) {
    ano <- Ano_CrossValid(var_exp = exp$data, var_obs = obs$data, memb = memb)
   
    #  Without cross-validation 
  } else {
    tmp <- Clim(var_exp = exp$data, var_obs = obs$data, memb = memb)
    if (!is.null(filter_span)) {
        if (is.numeric(filter_span)) {
            pos_dims <- names(dim(tmp$clim_exp))
            reorder <- match(pos_dims, c('ftime',
                             pos_dims[-which(pos_dims == 'ftime')]))
            tmp$clim_obs <- aperm(apply(tmp$clim_obs, c(1 : 
               length(dim(tmp$clim_obs)))[-which(names(dim(tmp$clim_obs)) == 'ftime')],
               .Loess, loess_span = filter_span), reorder)
            tmp$clim_exp <- aperm(apply(tmp$clim_exp, c(1 :
               length(dim(tmp$clim_exp)))[-which(names(dim(tmp$clim_exp)) == 'ftime')],
               .Loess, loess_span = filter_span), reorder)
        } else {
            warning("Paramater 'filter_span' is not numeric and any filter",
                    " is being applied.")
        }
    }
    if (memb) { 
      clim_exp <- tmp$clim_exp
      clim_obs <- tmp$clim_obs
    } else {
      clim_exp <- InsertDim(tmp$clim_exp, 2, dim_exp[2]) 
      clim_obs <- InsertDim(tmp$clim_obs, 2, dim_obs[2]) 
    }
   
    clim_exp <- InsertDim(clim_exp, 3, dim_exp[3]) 
    clim_obs <- InsertDim(clim_obs, 3, dim_obs[3]) 
    ano <- NULL    
    ano$ano_exp <- exp$data - clim_exp
    ano$ano_obs <- obs$data - clim_obs 
  }
  
  # Permuting back dimensions to original order
  if  (dim_anom != 3) {
    
    if (case_obs == 0) {
      ano$ano_exp <- aperm(ano$ano_exp, perm = dimperm)
    } 
    if (case_exp == 0) {
      ano$ano_obs <- aperm(ano$ano_obs, perm = dimperm)
    }
        
    #Updating back permuted dimensions
    dim_exp <- dim(exp$data)
    dim_obs <- dim(obs$data)
  }
  
  # Adding dimensions names
  attr(ano$ano_exp, 'dimensions') <- dimnames_data
  attr(ano$ano_obs, 'dimensions') <- dimnames_data
  exp$data <- ano$ano_exp
  obs$data <- ano$ano_obs
  
  #  Outputs
  # ~~~~~~~~~
    if (case_obs == 1) {
      return(obs)
    } 
    else if (case_exp == 1) {
      return(exp)
    }
    else {
    return(list(exp = exp, obs = obs)) 
  }
}
.Loess <- function(clim, loess_span) {
  data <- data.frame(ensmean = clim, day = 1 : length(clim))
  loess_filt <- loess(ensmean ~ day, data, span = loess_span)
  output <- predict(loess_filt)
  return(output)
}

