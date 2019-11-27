#'@rdname CST_Analogs
#'@title Downscaling using Analogs based on large scale fields.
#' 
#'@author M. Carmen Alvarez-Castro, \email{carmen.alvarez-castro@cmcc.it}
#'@author Nuria Perez-Zanon \email{nuria.perez@bsc.es}
#'
#'@description This function perform a downscaling using Analogs. To compute 
#'the analogs, the function search for days with similar large scale conditions
#'to downscaled fields in the local scale. The large scale and the local scale 
#'regions are defined by the user. The large scale is usually given by 
#'atmospheric circulation as sea level pressure or geopotential height 
#'(Yiou et al, 2013) but the function gives the possibility to use another 
#'field. The local scale will be usually given by precipitation or temperature 
#'fields, but might be another variable.The analogs function will find the best
#'analogs based in three criterias:
#' (1) Minimal distance in the large scale pattern (i.e. SLP)
#' (2) Minimal distance in the large scale pattern (i.e. SLP) and minimal 
#' distance in the local scale pattern (i.e. SLP). 
#' (3) Minimal distance in the large scale pattern (i.e. SLP), minimal 
#' distance in the local scale pattern (i.e. SLP) and maxima correlation in the 
#' local variable to downscale (i.e Precipitation).
#'The search of analogs must be done in the longest dataset posible. This is 
#'important since it is necessary to have a good representation of the 
#'possible states of the field in the past, and therefore, to get better 
#'analogs. Once the search of the analogs is complete, and in order to used 
#'the three criterias the user can select a number of analogs, using parameter
#''nAnalogs' to restrict the selection of the best analogs in a short number 
#'of posibilities, the best ones.
#'This function has not constrains of specific regions, variables to downscale,
#'or data to be used (seasonal forecast data, climate projections data, 
#'reanalyses data). The regrid into a finner scale is done interpolating with 
#'CST_Load. Then, this interpolation is corrected selecting the analogs in the 
#'large and local scale in based of the observations. The function is an 
#'adapted version of the method of Yiou et al 2013.       
#'
#'@references Yiou, P., T. Salameh, P. Drobinski, L. Menut, R. Vautard,
#' and M. Vrac, 2013 : Ensemble reconstruction of the atmospheric column 
#' from surface pressure using analogues.  Clim. Dyn., 41, 1419-1437. 
#' \email{pascal.yiou@lsce.ipsl.fr}
#'
#'@param expL an 's2dv_cube' object containing the experimental field on the 
#'large scale for which the analog is aimed. This field is used to in all the
#'criterias. If parameter 'expVar' is not provided, the function will return 
#'the expL analog. The element 'data' in the 's2dv_cube' object must have, at
#'least, latitudinal and longitudinal dimensions. The object is expect to be 
#'already subset for the desired large scale region.
#'@param obsL an 's2dv_cube' object containing the observational field on the
#'large scale. The element 'data' in the 's2dv_cube' object must have the same 
#'latitudinal and longitudinal dimensions as parameter 'expL' and a temporal 
#'dimension with the maximum number of available observations.
#'@param time_obsL a character string indicating the date of the observations 
#'in the format "dd/mm/yyyy"
#'@param expVar an 's2dv_cube' object containing the experimental field on the
#'local scale, usually a different variable to the parameter 'expL'. If it is 
#'not NULL (by default, NULL), the returned field by this function will be the 
#'analog of parameter 'expVar'.
#'@param obsVar an 's2dv_cube' containing the field of the same variable as the 
#'passed in parameter 'expVar' for the same region.
#'@param region a vector of length four indicating the minimum longitude, the 
#'maximum longitude, the minimum latitude and the maximum latitude.
#'@param criteria a character string indicating the criteria to be used for the 
#'selection of analogs:
#'\itemize{
#'\item{Large_dist} minimal distance in the large scale pattern;
#'\item{Local_dist} minimal distance in the large scale pattern and minimal 
#' distance in the local scale pattern; and
#'\item{Local_cor} minimal distance in the large scale pattern, minimal 
#' distance in the local scale pattern and maxima correlation in the 
#' local variable to downscale.}
#' 
#'@import multiApply
#'@importFrom ClimProjDiags SelBox
#'@import abind
#'
#'@seealso code{\link{CST_Load}}, \code{\link[s2dverification]{Load}} and 
#'\code{\link[s2dverification]{CDORemap}}
#'
#'@return An 's2dv_cube' object containing the dowscaled values of the best 
#'analogs in the criteria selected. 
#'
#'@examples
#'res <- CST_Analogs(expL = lonlat_data$exp, obsL = lonlat_data$obs)
#'
#'@export
CST_Analogs <- function(expL, obsL, time_obsL, expVar = NULL, obsVar = NULL,
                        region = NULL, criteria = "Large_dist") {
  if (!inherits(expL, 's2dv_cube') || !inherits(obsL, 's2dv_cube')) {
    stop("Parameter 'expL' and 'obsL' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load.")
  }
  if (!is.null(expVar) || !is.null(obsVar)) {
    if (!inherits(expVar, 's2dv_cube') || !inherits(obsVar, 's2dv_cube')) {
      stop("Parameter 'expVar' and 'obsVar' must be of the class 's2dv_cube',
           ","as output by CSTools::CST_Load.")
    }
  }
  timevector <- obsL$Dates$start
  if (!is.null(expVar)) {
    region <- c(min(expVar$lon), max(expVar$lon), min(expVar$lat), 
                max(expVar$lon))
    lonVar <- expVar$lon
    latVar <- expVar$lat
  } else {
    region <- c(min(expL$lon), max(expL$lon), min(expL$lat), max(expL$lon))
    lonVar <- expL$lon
    latVar <- expL$lat
  }      
  result <- Analogs(expL$data, obsL$data, time_obsL = timevector, 
                    expVar = expVar$data, obsVar = obsVar$data,
                    criteria = criteria, 
                    lonVar = expVar$lon, latVar = expVar$lat, 
                    region = region, nAnalogs = 1, return_list = FALSE)
  if (!is.null(obsVar)) {
    obsVar$data <- result$AnalogsFields
    return(obsVar)
  } else {
    obsL$data <- result$AnalogsFields
    return(obsL)
  }
}
#'@rdname Analogs
#'@title Analogs based on large scale fields.
#' 
#'@author M. Carmen Alvarez-Castro, \email{carmen.alvarez-castro@cmcc.it}
#'@author Nuria Perez-Zanon \email{nuria.perez@bsc.es}
#'        
#'@description This function perform a downscaling using Analogs. To compute 
#'the analogs, the function search for days with similar large scale conditions
#'to downscaled fields in the local scale. The large scale and the local scale 
#'regions are defined by the user. The large scale is usually given by 
#'atmospheric circulation as sea level pressure or geopotential height (Yiou 
#'et al, 2013) but the function gives the possibility to use another field. The
#'local scale will be usually given by precipitation or temperature fields, but
#'might be another variable. 
#'The analogs function will find the best analogs based in three criterias:
#' (1) Minimum Euclidean distance in the large scale pattern (i.e. SLP)
#' (2) Minimum Euclidean distance in the large scale pattern (i.e. SLP) 
#' and minimum Euclidean distance in the local scale pattern (i.e. SLP). 
#' (3) Minimum Euclidean distance in the large scale pattern (i.e. SLP), minimum 
#' distance in the local scale pattern (i.e. SLP) and highest correlation in the
#' local variable to downscale (i.e Precipitation).
#'The search of analogs must be done in the longest dataset posible. This is 
#'important since it is necessary to have a good representation of the 
#'possible states of the field in the past, and therefore, to get better 
#'analogs. Once the search of the analogs is complete, and in order to used the 
#'three criterias the user can select a number of analogs , using parameter 
#''nAnalogs' to restrict 
#'the selection of the best analogs in a short number of posibilities, the best
#'ones. This function has not constrains of specific regions, variables to 
#'downscale, or data to be used (seasonal forecast data, climate projections 
#'data, reanalyses data). The regrid into a finner scale is done interpolating 
#'with CST_Load. Then, this interpolation is corrected selecting the analogs in
#'the large and local scale in based of the observations. The function is an 
#'adapted version of the method of Yiou et al 2013.       
#'
#'@references Yiou, P., T. Salameh, P. Drobinski, L. Menut, R. Vautard,
#'and M. Vrac, 2013 : Ensemble reconstruction of the atmospheric column 
#'from surface pressure using analogues.  Clim. Dyn., 41, 1419-1437. 
#'\email{pascal.yiou@lsce.ipsl.fr}
#'
#'@param expL an array of N named dimensions containing the experimental field
#'on the large scale for which the analog is aimed. This field is used to in 
#'all the criterias. If parameter 'expVar' is not provided, the function will
#'return the expL analog. The element 'data' in the 's2dv_cube' object must 
#'have, at least, latitudinal and longitudinal dimensions. The object is expect
#'to be already subset for the desired large scale region.
#'@param obsL an array of N named dimensions containing the observational field 
#'on the large scale. The element 'data' in the 's2dv_cube' object must have 
#'the same latitudinal and longitudinal dimensions as parameter 'expL' and a
#' temporal dimension with the maximum number of available observations.
#'@param time_obsL a character string indicating the date of the observations 
#'in the format "dd/mm/yyyy"
#'@param expVar an array of N named dimensions containing the experimental 
#'field on the local scale, usually a different variable to the parameter 
#''expL'. If it is not NULL (by default, NULL), the returned field by this 
#'function will be the analog of parameter 'expVar'.
#'@param obsVar an array of N named dimensions containing the field of the 
#'same variable as the passed in parameter 'expVar' for the same region.
#'@param criteria a character string indicating the criteria to be used for the
#'selection of analogs:
#'\itemize{
#'\item{Large_dist} minimum Euclidean distance in the large scale pattern; 
#'\item{Local_dist} minimum Euclidean distance in the large scale pattern 
#'and minimum Euclidean distance in the local scale pattern; and
#'\item{Local_cor} minimum Euclidean distance in the large scale pattern, 
#'minimum Euclidean distance in the local scale pattern and highest correlation
#' in the local variable to downscale.} 
#'@param lonVar a vector containing the longitude of parameter 'expVar'.
#'@param latVar a vector containing the latitude of parameter 'expVar'.
#'@param region a vector of length four indicating the minimum longitude, 
#'the maximum longitude, the minimum latitude and the maximum latitude.
#'@param return_list TRUE to get a list with the best analogs. FALSE
#'to get a single analog, the best analog. For Downscaling return_list must be
#'FALSE. 
#'@param nAnalogs number of Analogs to be selected to apply the criterias 
#''Local_dist' or 'Local_cor'. This is not the necessary the number of analogs 
#'that the user can get, but the number of events with minimum distance in 
#'which perform the search of the best Analog. The default value for the 
#''Large_dist' criteria is 1, for 'Local_dist' and 'Local_cor'criterias must
#' be selected by the user otherwise the default value will be set as 10. 
#'@import multiApply
#'@importFrom ClimProjDiags SelBox
#'@import abind
#'@return AnalogsFields, dowscaled values of the best analogs for the criteria 
#'selected.
#'@return AnalogsInfo, a dataframe with information about the number of the 
#'best analogs, the corresponding value of the metric used in the selected 
#'criteria (distance values for Large_dist and Local_dist,correlation values 
#'for Local_cor), date of the analog). The analogs are listed in decreasing 
#'order, the first one is the best analog (i.e if the selected criteria 
#'is Local_cor the best analog will be the one with highest correlation, while
#'for Large_dist criteria the best analog will be the day with minimum 
#'Euclidean distance)
#'@examples
#'require(zeallot)
#'
#'# Example 1:Downscaling using criteria 'Large_dist' and a single variable:
#'# The best analog is found using a single variable (i.e. Sea level pressure, 
#'# SLP). The downscaling will be done in the same variable used to study the 
#'# minimal distance (i.e. SLP). obsVar and expVar NULLS or equal to obsL and 
#'# expL respectively region, lonVar and latVar not necessary here. 
#'# return_list=FALSE 
#'
#'expSLP <- rnorm(1:20)
#'dim(expSLP) <- c(lat = 4, lon = 5)
#'obsSLP <- c(rnorm(1:180),expSLP*1.2)
#'dim(obsSLP) <- c(lat = 4, lon = 5, time = 10)
#'time_obsSLP <- paste(rep("01", 10), rep("01", 10), 1994 : 2003, sep = "-")
#'downscale_field <- Analogs(expL=expSLP, obsL=obsSLP, time_obsL=time_obsSLP)
#'str(downscale_field)
#'
#'# Example 2: Downscaling using criteria 'Large_dist' and 2 variables:
#'# The best analog is found using 2 variables (i.e. Sea Level Pressure, SLP 
#'# and precipitation, pr): one variable (i.e. sea level pressure, expL) to 
#'# find the best analog day (defined in criteria 'Large_dist' as the day, in 
#'# obsL, with the minimum Euclidean distance to the day of interest in expL) 
#'# The second variable will be the variable to donwscale (i.e. precipitation, 
#'# obsVar). Parameter obsVar must be different to obsL.The downscaled 
#'# precipitation will be the precipitation that belongs to the best analog day
#'# in SLP. Region not needed here since will be the same for both variables.
#'
#'expSLP <- rnorm(1:20)
#'dim(expSLP) <- c(lat = 4, lon = 5)
#'obsSLP <- c(rnorm(1:180),expSLP*2)
#'dim(obsSLP) <- c(lat = 4, lon = 5, time = 10)
#'time_obsSLP <- paste(rep("01", 10), rep("01", 10), 1994 : 2003, sep = "-")
#'obs.pr <- c(rnorm(1:200)*0.001)
#'dim(obs.pr)=dim(obsSLP)
#'downscale_field <- Analogs(expL=expSLP, obsL=obsSLP, 
#'                           obsVar=obs.pr,
#'                           time_obsL=time_obsSLP)
#'str(downscale_field)
#'
#'# Example 3:List of best Analogs using criteria 'Large_dist' and a single 
#'# variable: same as Example 1 but getting a list of best analogs (return_list
#'# =TRUE) with the corresponding downscaled values, instead of only 1 single 
#'# donwscaled value instead of 1 single downscaled value. Imposing nAnalogs 
#'# (number of analogs to do the search of the best Analogs). obsVar and expVar
#'# NULL or equal to obsL and expL,respectively.
#'
#'expSLP <- rnorm(1:20)
#'dim(expSLP) <- c(lat = 4, lon = 5)
#'obsSLP <- c(rnorm(1:1980),expSLP*1.5)
#'dim(obsSLP) <- c(lat = 4, lon = 5, time = 100)
#'time_obsSLP <- paste(rep("01", 100), rep("01", 100), 1920 : 2019, sep = "-")
#'downscale_field<- Analogs(expL=expSLP, obsL=obsSLP, time_obsSLP,
#'                          nAnalogs=5,return_list = TRUE)
#'str(downscale_field)
#'
#'# Example 4:List of best Analogs using criteria 'Large_dist' and 2 variables:
#'# same as example 2 but getting a list of best analogs (return_list =TRUE) 
#'# with the values instead of only a single downscaled value. Imposing 
#'# nAnalogs (number of analogs to do the search of the best Analogs). obsVar 
#'# and expVar must be different to obsL and expL.
#'
#'expSLP <- rnorm(1:20)
#'dim(expSLP) <- c(lat = 4, lon = 5)
#'obsSLP <- c(rnorm(1:180),expSLP*2)
#'dim(obsSLP) <- c(lat = 4, lon = 5, time = 10)
#'time_obsSLP <- paste(rep("01", 10), rep("01", 10), 1994 : 2003, sep = "-")
#'obs.pr <- c(rnorm(1:200)*0.001)
#'dim(obs.pr)=dim(obsSLP)
#'downscale_field <- Analogs(expL=expSLP, obsL=obsSLP, 
#'                           obsVar=obs.pr,
#'                           time_obsL=time_obsSLP,nAnalogs=5,
#'                           return_list = TRUE)
#'str(downscale_field)
#'
#'# Example 5: Downscaling using criteria 'Local_dist' and 2 variables:
#'# The best analog is found using 2 variables (i.e. Sea Level Pressure, 
#'# SLP and precipitation, pr). Parameter obsVar must be different to obsL.The 
#'# downscaled precipitation will be the precipitation that belongs to the best 
#'# analog day in SLP. lonVar, latVar and Region must be given here to select 
#'# the local scale
#'
#'expSLP <- rnorm(1:20)
#'dim(expSLP) <- c(lat = 4, lon = 5)
#'obsSLP <- c(rnorm(1:180),expSLP*2)
#'dim(obsSLP) <- c(lat = 4, lon = 5, time = 10)
#'time_obsSLP <- paste(rep("01", 10), rep("01", 10), 1994 : 2003, sep = "-")
#'obs.pr <- c(rnorm(1:200)*0.001)
#'dim(obs.pr)=dim(obsSLP)
#'# analogs of local scale using criteria 2
#'lonmin=-1
#'lonmax=2
#'latmin=30
#'latmax=33
#'region=c(lonmin,lonmax,latmin,latmax)
#'Local_scale <- Analogs(expL=expSLP,
#'                       obsL=obsSLP, time_obsL=time_obsSLP,
#'                       obsVar=obs.pr,
#'                       criteria="Local_dist",lonVar=seq(-1,5,1.5),
#'                       latVar=seq(30,35,1.5),region=region, 
#'                       nAnalogs = 10, return_list = FALSE)
#'str(Local_scale)
#'
#'# Example 6: list of best analogs using criteria 'Local_dist' and 2 
#'# variables:
#'# The best analog is found using 2 variables. Parameter obsVar must be 
#'# different to obsL in this case.The downscaled precipitation will be the 
#'# precipitation that belongs to the best analog day in SLP. lonVar, latVar 
#'# and Region needed. return_list=TRUE
#'
#'expSLP <- rnorm(1:20)
#'dim(expSLP) <- c(lat = 4, lon = 5)
#'obsSLP <- c(rnorm(1:180),expSLP*2)
#'dim(obsSLP) <- c(lat = 4, lon = 5, time = 10)
#'time_obsSLP <- paste(rep("01", 10), rep("01", 10), 1994 : 2003, sep = "-")
#'obs.pr <- c(rnorm(1:200)*0.001)
#'dim(obs.pr)=dim(obsSLP)
#'# analogs of local scale using criteria 2
#'lonmin=-1
#'lonmax=2
#'latmin=30
#'latmax=33
#'region=c(lonmin,lonmax,latmin,latmax)
#'Local_scale <- Analogs(expL=expSLP,
#'                       obsL=obsSLP, time_obsL=time_obsSLP,
#'                       obsVar=obs.pr,
#'                       criteria="Local_dist",lonVar=seq(-1,5,1.5),
#'                       latVar=seq(30,35,1.5),region=region, 
#'                       nAnalogs = 5, return_list = TRUE)
#'str(Local_scale)
#'
#'# Example 7: Downscaling using Local_dist criteria
#'# without parameters obsVar and expVar. return list =FALSE 
#'
#'expSLP <- rnorm(1:20)
#'dim(expSLP) <- c(lat = 4, lon = 5)
#'obsSLP <- c(rnorm(1:180),expSLP*2)
#'dim(obsSLP) <- c(lat = 4, lon = 5, time = 10)
#'time_obsSLP <- paste(rep("01", 10), rep("01", 10), 1994 : 2003, sep = "-")
#'# analogs of local scale using criteria 2
#'lonmin=-1
#'lonmax=2
#'latmin=30
#'latmax=33
#'region=c(lonmin,lonmax,latmin,latmax)
#'Local_scale <- Analogs(expL=expSLP,
#'                       obsL=obsSLP, time_obsL=time_obsSLP,
#'                       criteria="Local_dist",lonVar=seq(-1,5,1.5),
#'                       latVar=seq(30,35,1.5),region=region, 
#'                       nAnalogs = 10, return_list = TRUE)
#'str(Local_scale)
#'
#'# Example 8: Downscaling using criteria 'Local_cor' and 2 variables:
#'# The best analog is found using 2 variables. Parameter obsVar and expVar 
#'# are necessary and must be different to obsL and expL, respectively.
#'# The downscaled precipitation will be the precipitation that belongs to
#'# the best analog day in SLP large and local scales, and to the day with 
#'# the higher correlation in precipitation. return_list=FALSE. two options 
#'# for nAnalogs
#'
#'expSLP <- rnorm(1:20)
#'dim(expSLP) <- c(lat = 4, lon = 5)
#'obsSLP <- c(rnorm(1:180),expSLP*2)
#'dim(obsSLP) <- c(lat = 4, lon = 5, time = 10)
#'time_obsSLP <- paste(rep("01", 10), rep("01", 10), 1994 : 2003, sep = "-")
#'exp.pr <- c(rnorm(1:20)*0.001)
#'dim(exp.pr)=dim(expSLP)
#'obs.pr <- c(rnorm(1:200)*0.001)
#'dim(obs.pr)=dim(obsSLP)
#'# analogs of local scale using criteria 2
#'lonmin=-1
#'lonmax=2
#'latmin=30
#'latmax=33
#'region=c(lonmin,lonmax,latmin,latmax)
#'Local_scalecor <- Analogs(expL=expSLP,
#'                          obsL=obsSLP, time_obsL=time_obsSLP,
#'                          obsVar=obs.pr,expVar=exp.pr,
#'                          criteria="Local_cor",lonVar=seq(-1,5,1.5),
#'                          latVar=seq(30,35,1.5),nAnalogs=8,region=region, 
#'                          return_list = FALSE)
#'Local_scalecor$AnalogsInfo
#'Local_scalecor$DatesAnalogs
#'# same but without imposing nAnalogs, so nAnalogs will be set by default as 10 
#'Local_scalecor <- Analogs(expL=expSLP,
#'                          obsL=obsSLP, time_obsL=time_obsSLP,
#'                          obsVar=obs.pr,expVar=exp.pr,
#'                          criteria="Local_cor",lonVar=seq(-1,5,1.5),
#'                          latVar=seq(30,35,1.5),region=region, 
#'                          return_list = FALSE)
#'Local_scalecor$AnalogsInfo
#'Local_scalecor$DatesAnalogs
#'
#'# Example 9: List of best analogs in the three criterias Large_dist, 
#'# Local_dist, and Local_cor return list TRUE same variable
#'
#'expSLP <- rnorm(1:20)
#'dim(expSLP) <- c(lat = 4, lon = 5)
#'obsSLP <- c(rnorm(1:180),expSLP*2)
#'dim(obsSLP) <- c(lat = 4, lon = 5, time = 10)
#'time_obsSLP <- paste(rep("01", 10), rep("01", 10), 1994 : 2003, sep = "-")
#'# analogs of large scale using criteria 1
#'Large_scale <- Analogs(expL=expSLP,
#'                       obsL=obsSLP, time_obsL=time_obsSLP,
#'                       criteria="Large_dist", 
#'                       nAnalogs = 7, return_list = TRUE)
#'str(Large_scale)
#'Large_scale$AnalogsInfo
#'# analogs of local scale using criteria 2
#'lonmin=-1
#'lonmax=2
#'latmin=30
#'latmax=33
#'region=c(lonmin,lonmax,latmin,latmax)
#'Local_scale <- Analogs(expL=expSLP,
#'                       obsL=obsSLP, time_obsL=time_obsSLP,
#'                       criteria="Local_dist",lonVar=seq(-1,5,1.5),
#'                       latVar=seq(30,35,1.5),nAnalogs=7,region=region, 
#'                       return_list = TRUE)
#'str(Local_scale)
#'Local_scale$AnalogsInfo
#'# analogs of local scale using criteria 3
#'Local_scalecor <- Analogs(expL=expSLP,
#'                          obsL=obsSLP, time_obsL=time_obsSLP,
#'                          obsVar=obsSLP,expVar=expSLP,
#'                          criteria="Local_cor",lonVar=seq(-1,5,1.5),
#'                          latVar=seq(30,35,1.5),nAnalogs=7,region=region, 
#'                          return_list = TRUE)
#'str(Local_scalecor)
#'Local_scalecor$AnalogsInfo
#'
#'# Example 10: Downscaling in the three criteria Large_dist, Local_dist, and 
#'# Local_cor return list FALSE, different variable
#'
#'expSLP <- rnorm(1:20)
#'dim(expSLP) <- c(lat = 4, lon = 5)
#'obsSLP <- c(rnorm(1:180),expSLP*2)
#'dim(obsSLP) <- c(lat = 4, lon = 5, time = 10)
#'time_obsSLP <- paste(rep("01", 10), rep("01", 10), 1994 : 2003, sep = "-")
#'exp.pr <- c(rnorm(1:20)*0.001)
#'dim(exp.pr)=dim(expSLP)
#'obs.pr <- c(rnorm(1:200)*0.001)
#'dim(obs.pr)=dim(obsSLP)
#'# analogs of large scale using criteria 1
#'Large_scale <- Analogs(expL=expSLP,
#'                       obsL=obsSLP, time_obsL=time_obsSLP,
#'                       criteria="Large_dist", 
#'                       nAnalogs = 7, return_list = FALSE)
#'str(Large_scale)
#'Large_scale$AnalogsInfo
#'# analogs of local scale using criteria 2
#'lonmin=-1
#'lonmax=2
#'latmin=30
#'latmax=33
#'region=c(lonmin,lonmax,latmin,latmax)
#'Local_scale <- Analogs(expL=expSLP,
#'                       obsL=obsSLP, time_obsL=time_obsSLP,
#'                       obsVar=obs.pr,
#'                       criteria="Local_dist",lonVar=seq(-1,5,1.5),
#'                       latVar=seq(30,35,1.5),nAnalogs=7,region=region, 
#'                       return_list = FALSE)
#'str(Local_scale)
#'Local_scale$AnalogsInfo
#'# analogs of local scale using criteria 3
#'Local_scalecor <- Analogs(expL=expSLP,
#'                          obsL=obsSLP, time_obsL=time_obsSLP,
#'                          obsVar=obs.pr,expVar=exp.pr,
#'                          criteria="Local_cor",lonVar=seq(-1,5,1.5),
#'                          latVar=seq(30,35,1.5),nAnalogs=7,region=region, 
#'                          return_list = FALSE)
#'str(Local_scalecor)
#'Local_scalecor$AnalogsInfo
#'
#'@export 
Analogs <- function(expL, obsL, time_obsL, expVar = NULL, obsVar = NULL, 
                    criteria = "Large_dist",
                    lonVar = NULL, latVar = NULL, region = NULL, 
                    nAnalogs = NULL, return_list = FALSE) {
  # checks
  if (!all(c('lon', 'lat') %in% names(dim(expL)))) {
    stop("Parameter 'expL' must have the dimensions 'lat' and 'lon'.")
  }
  
  if (!all(c('lat', 'lon') %in% names(dim(obsL)))) {
    stop("Parameter 'obsL' must have the dimension 'lat' and 'lon'.")
  }
  
  if (any(is.na(expL)))  {
    warning("Parameter 'exp' contains NA values.")
  }
  
  if (any(is.na(obsL))) {
    warning("Parameter 'obs' contains NA values.")
  }
  if (!is.null(expVar) & is.null(obsVar)) {
    expVar <- NULL
    warning("Parameter 'expVar' is set to NULL as parameter 'obsVar', 
            large scale field will be returned.")
  }
  if (is.null(expVar) & is.null(obsVar)) {
    warning("Parameter 'expVar' and 'obsVar' are NULLs, downscaling/listing 
            same variable as obsL and expL'.")
  }
  if(!is.null(obsVar) & is.null(expVar) & criteria=="Local_cor"){
    stop("parameter 'expVar' cannot be NULL")
  }
  if(is.null(nAnalogs) & criteria!="Large_dist"){
    nAnalogs=10
    warning("Parameter 'nAnalogs' is NULL and is set to 10 by default")
  }
  if(is.null(nAnalogs) & criteria=="Large_dist"){
    nAnalogs=1
  }
  
  if (any(names(dim(obsL)) %in% 'ftime')) {
    if (any(names(dim(obsL)) %in% 'time')) {
      stop("Multiple temporal dimensions ('ftime' and 'time') found",
           "in parameter 'obsL'.")
    } else {
      time_pos_obsL <-  which(names(dim(obsL)) == 'ftime')
      names(dim(obsL))[time_pos_obsL] <- 'time'
      if (any(names(dim(expL)) %in% 'ftime')) {
        time_pos_expL <-  which(names(dim(expL)) == 'ftime')
        names(dim(expL))[time_pos_expL] <- 'time'
      }
    }
  }
  if (any(names(dim(obsVar)) %in% 'ftime')) {
    if (any(names(dim(obsVar)) %in% 'time')) {
      stop("Multiple temporal dimensions ('ftime' and 'time') found",
           "in parameter 'obsVar'.")
    } else {
      time_pos_obsVar <-  which(names(dim(obsVar)) == 'ftime')
      names(dim(obsVar))[time_pos_obsVar] <- 'time'
      if (any(names(dim(expVar)) %in% 'ftime')) {
        time_pos_expVar <-  which(names(dim(expVar)) == 'ftime')
        names(dim(expVar))[time_pos_expVar] <- 'time'
      }
    }
  }
  if (any(names(dim(obsL)) %in% 'sdate')) {
    if (any(names(dim(obsL)) %in% 'time')) {
      dims_obsL <- dim(obsL)
      pos_sdate <- which(names(dim(obsL)) == 'sdate')
      pos_time <- which(names(dim(obsL)) == 'time')
      pos <- 1 : length(dim(obsL))
      pos <- c(pos_time, pos_sdate, pos[-c(pos_sdate,pos_time)])
      obsL <- aperm(obsL, pos)
      dim(obsL) <- c(time = prod(dims_obsL[c(pos_time, pos_sdate)]),
                     dims_obsL[-c(pos_time, pos_sdate)])
    } else {
      stop("Parameter 'obsL' must have a temporal dimension.")
    }
  }
  if (any(names(dim(obsVar)) %in% 'sdate')) {
    if (any(names(dim(obsVar)) %in% 'time')) {
      dims_obsVar <- dim(obsVar)
      pos_sdate <- which(names(dim(obsVar)) == 'sdate')
      pos_time <- which(names(dim(obsVar)) == 'time')
      pos <- 1 : length(dim(obsVar))
      pos <- c(pos_time, pos_sdate, pos[-c(pos_sdate,pos_time)])
      obsVar <- aperm(obsVar, pos)
      dim(obsVar) <- c(time = prod(dims_obsVar[c(pos_time, pos_sdate)]),
                     dims_obsVar[-c(pos_time, pos_sdate)])
    } else {
      stop("Parameter 'obsVar' must have a temporal dimension.")
    }
  } 
  
  if (is.null(region) & criteria!="Large_dist") {
    if (!is.null(lonVar) & !is.null(latVar)) {
      region <- c(min(lonVar), max(lonVar), min(latVar), max(latVar))
    }else{
      stop("Parameters 'lonVar' and 'latVar' must be given in criteria 
           'Local_dist' and 'Local_cor'")
    }
  }
  position <- Select(expL = expL, obsL = obsL,  expVar = expVar, 
                     obsVar = obsVar, criteria = criteria, lonVar = lonVar, 
                     latVar = latVar, region = region)$position
  metrics<- Select(expL = expL, obsL = obsL,  expVar = expVar, 
                    obsVar = obsVar, criteria = criteria, lonVar = lonVar, 
                    latVar = latVar, region = region)$metric.original
  best <- Apply(list(position), target_dims = c('time', 'pos'), 
                fun = BestAnalog, criteria = criteria, 
                return_list = return_list, nAnalogs = nAnalogs)$output1
  Analogs_dates <- time_obsL[best]
  dim(Analogs_dates) <- dim(best)
  if (all(!is.null(region), !is.null(lonVar), !is.null(latVar))) {
    if (is.null(obsVar)) {
      obsVar <- SelBox(obsL, lon = lonVar, lat = latVar, region = region)$data
      expVar <- SelBox(expL, lon = lonVar, lat = latVar, region=region)$data
      Analogs_fields <- Subset(obsVar, 
                               along = which(names(dim(obsVar)) == 'time'),
                               indices = best)
      warning("obsVar is NULL, 
          obsVar will be computed from obsL (same variable)")
      
    } else {
      obslocal <- SelBox(obsVar, lon = lonVar, lat = latVar, 
                         region = region)$data
      Analogs_fields <- Subset(obslocal, 
                               along = which(names(dim(obslocal)) == 'time'),
                               indices = best)
    }
  } else {
    warning("One or more of the parameter 'region', 'lonVar' and 'latVar'",
            " are NULL and the large scale field will be returned.")
    if (is.null(obsVar)) {
      Analogs_fields <- Subset(obsL, along = which(names(dim(obsL)) == 'time'),
                               indices = best)
    } else {
      Analogs_fields <- Subset(obsVar,
                               along = which(names(dim(obsVar)) == 'time'),
                               indices = best)
    }
  }
  
  lon_dim <- which(names(dim(Analogs_fields)) == 'lon')
  lat_dim <- which(names(dim(Analogs_fields)) == 'lat')
  if (lon_dim < lat_dim) {
  dim(Analogs_fields) <- c(dim(Analogs_fields)[c(lon_dim, lat_dim)], dim(best))
  } else if (lon_dim > lat_dim) {
  dim(Analogs_fields) <- c(dim(Analogs_fields)[c(lat_dim, lon_dim)], dim(best))
  } else {
    stop("Dimensions 'lat' and 'lon' not found.")
  }
  Analogs_metrics <- Subset(metrics,
                            along = which(names(dim(metrics)) == 'time'),
                           indices = best)
  DistCorr <- data.frame(as.numeric(1:nrow(Analogs_metrics)),(Analogs_metrics),
                    Analogs_dates)
  if(dim(DistCorr)[2]==3){names(DistCorr) <- c("Analog","LargeDist","Dates")}
  if(dim(DistCorr)[2]==4){names(DistCorr) <- c("Analog","LargeDist",
                                                  "LocalDist","Dates")}
  if(dim(DistCorr)[2]==5){names(DistCorr) <- c("Analog","LargeDist",
                                                  "LocalDist","LocalCorr","Dates")}
  return(list(AnalogsFields = Analogs_fields,
              AnalogsInfo=DistCorr))
  }

BestAnalog <- function(position, criteria = 'Large_dist', return_list = FALSE, 
                       nAnalogs = nAnalogs) {
  pos_dim <- which(names(dim(position)) == 'pos')
  if (dim(position)[pos_dim] == 1) {
    pos1 <- Subset(position, along = pos_dim, indices = 1)
    if (criteria != 'Large_dist') {
      warning("Dimension 'pos' in parameter 'position' has length 1,",
              " criteria 'Large_dist' will be used.")
      criteria <- 'Large_dist'
    }
  } else if (dim(position)[pos_dim] == 2) {
    pos1 <- Subset(position, along = pos_dim, indices = 1)
    pos2 <- Subset(position, along = pos_dim, indices = 2)
    if (criteria == 'Local_cor') {
      warning("Dimension 'pos' in parameter 'position' has length 2,",
              " criteria 'Local_dist' will be used.")
      criteria <- 'Local_dist'
    }
  } else if (dim(position)[pos_dim] == 3) {
    pos1 <- Subset(position, along = pos_dim, indices = 1)
    pos2 <- Subset(position, along = pos_dim, indices = 2)
    pos3 <- Subset(position, along = pos_dim, indices = 3)
    if (criteria != 'Local_cor') {
      warning("Parameter 'criteria' is set to", criteria, ".")
    }
  } else {
    stop("Parameter 'position' has dimension 'pos' of different ",
         "length than expected (from 1 to 3).")
  }
  if (criteria == 'Large_dist') {
    if (return_list == FALSE) {
      pos <- pos1[1]
    } else {
      pos <- pos1[1 : nAnalogs]
    }
  } else if (criteria== 'Local_dist') {
    pos1 <- pos1[1 : nAnalogs]
    pos2 <- pos2[1 : nAnalogs]
    best <- match(pos1, pos2)
    if(length(best)==1){
      warning("Just 1 best analog matching Large_dist and ", 
              "Local_dist criteria")
    } 
    if(length(best)==1 & is.na(best[1])==TRUE){
      stop("no best analogs matching Large_dist and Local_dist criterias") 
    }
    pos <- pos2[as.logical(best)]
    pos <- pos[which(!is.na(pos))]
    if (return_list == FALSE) { 
      pos <- pos[1]
    }else {
        pos <- pos}
  } else if (criteria == 'Local_cor') {
    pos1 <- pos1[1 : nAnalogs]
    pos2 <- pos2[1 : nAnalogs]
    best <- match(pos1, pos2)
    pos <- pos1[as.logical(best)]
    pos <- pos[which(!is.na(pos))]
    pos3 <- pos3[1 : nAnalogs]
    best <- match(pos, pos3)
    pos <- pos[order(best, decreasing = F)]
    pos <- pos[which(!is.na(pos))]
    if (return_list == FALSE) { 
      pos <- pos[1]
    } else{
      pos <- pos
    }
    return(pos)
  }
}
Select <- function(expL, obsL,  expVar = NULL, obsVar = NULL, 
                   criteria = "Large_dist",
                   lonVar = NULL, latVar = NULL, region = NULL) {
names(dim(expL)) <- replace_repeat_dimnames(names(dim(expL)), names(dim(obsL)))
  metric1 <- Apply(list(obsL), target_dims = list(c('lat', 'lon')), 
                   fun = .select, expL, metric = "dist")$output1
  metric1.original=metric1
  if (length(dim(metric1)) > 1) {
    dim_time_obs <- which(names(dim(metric1)) == 'time' | 
                            names(dim(metric1)) == 'ftime')
    dim(metric1) <- c(dim(metric1), metric=1)
    margins <- c(1 : (length(dim(metric1))))[-dim_time_obs]
    pos1 <- apply(metric1, margins, order)      
    names(dim(pos1))[1] <- 'time'
    metric1.original=metric1
    metric1 <-  apply(metric1, margins, sort)
    names(dim(metric1))[1] <- 'time'
    names(dim(metric1.original))=names(dim(metric1))
  } else {
    pos1 <- order(metric1)
    dim(pos1) <- c(time = length(pos1))
    metric1 <- sort(metric1)
    dim(metric1) <- c(time = length(metric1))
    dim(metric1.original)=dim(metric1)
    dim_time_obs=1
  }
  if (criteria == "Large_dist") {
    dim(metric1) <- c(dim(metric1), metric = 1)
    dim(pos1) <- c(dim(pos1), pos = 1)
    dim(metric1.original)=dim(metric1)
    return(list(metric = metric1, metric.original=metric1.original,position = pos1))
  }
  if (criteria == "Local_dist" | criteria == "Local_cor") {
    obs <- SelBox(obsL, lon = lonVar, lat = latVar, region = region)$data
    exp <- SelBox(expL, lon = lonVar, lat = latVar, region = region)$data
    metric2 <- Apply(list(obs), target_dims = list(c('lat', 'lon')), 
                     fun = .select, exp, metric = "dist")$output1
    metric2.original=metric2
    dim(metric2) <- c(dim(metric2), metric=1)
    margins <- c(1 : (length(dim(metric2))))[-dim_time_obs]
    pos2 <- apply(metric2, margins, order)
    dim(pos2) <- dim(pos1)
    names(dim(pos2))[1] <- 'time'
    metric2 <-  apply(metric2, margins, sort)
    names(dim(metric2))[1] <- 'time'
    if (criteria == "Local_dist") {
      metric <- abind(metric1, metric2, along = length(dim(metric1))+1)
      metric.original <- abind(metric1.original,metric2.original,along=length(dim(metric1))+1)
      position <- abind(pos1, pos2, along = length(dim(pos1))+1) 
      names(dim(metric)) <- c(names(dim(pos1)), 'metric')
      names(dim(position)) <- c(names(dim(pos1)), 'pos')
      names(dim(metric.original)) = names(dim(metric)) 
      return(list(metric = metric, metric.original=metric.original,position = position))
    }   
  }
  if (criteria == "Local_cor") {
    obs <- SelBox(obsVar, lon = lonVar, lat = latVar, region = region)$data
    exp <- SelBox(expVar, lon = lonVar, lat = latVar, region = region)$data
    metric3 <- Apply(list(obs), target_dims = list(c('lat', 'lon')), 
                     fun = .select, exp, metric = "cor")$output1
    metric3.original=metric3
    dim(metric3) <- c(dim(metric3), metric=1)
    margins <- c(1 : (length(dim(metric3))))[-dim_time_obs]
    pos3 <- apply(abs(metric3), margins, order, decreasing = TRUE)
    names(dim(pos3))[1] <- 'time'
    #metric3 <-  apply(abs(metric3), margins, sort)
    metricsort <- metric3[pos3]
    dim(metricsort)=dim(metric3)
    names(dim(metricsort))[1] <- 'time'
    metric <- abind(metric1, metric2, metricsort, 
                    along = length(dim(metric1)) + 1)
    metric.original <- abind(metric1.original, metric2.original, metric3.original, 
                    along = length(dim(metric1)) + 1)
    position <- abind(pos1, pos2, pos3, along = length(dim(pos1)) + 1)   
    names(dim(metric)) <- c(names(dim(metric1)), 'metric')
    names(dim(position)) <- c(names(dim(pos1)), 'pos')
    names(dim(metric.original)) = names(dim(metric)) 
    return(list(metric = metric, metric.original=metric.original,position = position))
  }
  else {
    stop("Parameter 'criteria' must to be one of the: 'Large_dist', ",
         "'Local_dist','Local_cor'.")
  }
}
.select <- function(exp, obs, metric = "dist") {
  if (metric == "dist") {
    result <- Apply(list(obs), target_dims = list(c('lat', 'lon')), 
                    fun = function(x) {sum((x - exp) ^ 2)})$output1
  } else if (metric == "cor") {
    result <- Apply(list(obs), target_dims = list(c('lat', 'lon')), 
                    fun = function(x) {cor(as.vector(x), 
                                           as.vector(exp))})$output1
  } 
  result
}
replace_repeat_dimnames <- function(names_exp, names_obs, lat_name = 'lat', 
                                    lon_name = 'lon') {
  if (!is.character(names_exp)) {
    stop("Parameter 'names_exp' must be a vector of characters.")
  }
  if (!is.character(names_obs)) {
    stop("Parameter 'names_obs' must be a vector of characters.")
  }
  latlon_dim_exp <- which(names_exp == lat_name | names_exp == lon_name)
  latlon_dim_obs <- which(names_obs == lat_name | names_obs == lon_name)
  if (any(unlist(lapply(names_exp[-latlon_dim_exp],
                        function(x){x == names_obs[-latlon_dim_obs]})))) {
    original_pos <- lapply(names_exp, 
                           function(x) which(x == names_obs[-latlon_dim_obs]))
    original_pos <- lapply(original_pos, length) > 0
    names_exp[original_pos] <- paste0(names_exp[original_pos], "_exp")
  }
  return(names_exp)
  ## Improvements: other dimensions to avoid replacement for more flexibility.
}
