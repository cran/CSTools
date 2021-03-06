#' @rdname CST_WeatherRegimes
#' @title Function for Calculating the Cluster analysis
#'
#' @author Verónica Torralba - BSC, \email{veronica.torralba@bsc.es}
#'
#' @description This function computes the weather regimes from a cluster analysis.
#'It is applied on the array \code{data} in a 's2dv_cube' object. The dimensionality of this object can be also reduced
#'by using PCs obtained from the application of the #'EOFs analysis to filter the dataset.
#'The cluster analysis can be performed with the traditional k-means or those methods
#'included in the hclust (stats package).
#'
#'@references Cortesi, N., V., Torralba, N., González-Reviriego, A., Soret, and F.J., Doblas-Reyes (2019).
#' Characterization of European wind speed variability using weather regimes. Climate Dynamics,53,
#' 4961–4976, doi:10.1007/s00382-019-04839-5.
#'@references Torralba, V. (2019) Seasonal climate prediction for the wind energy sector: methods and tools
#' for the development of a climate service. Thesis. Available online: \url{https://eprints.ucm.es/56841/}
#'
#'@param data a 's2dv_cube' object
#'@param ncenters Number of clusters to be calculated with the clustering function.
#'@param EOFs Whether to compute the EOFs (default = 'TRUE') or not (FALSE) to filter the data.
#'@param neofs number of modes to be kept (default = 30).
#'@param varThreshold Value with the percentage of variance to be explained by the PCs.
#' Only sufficient PCs to explain this much variance will be used in the clustering.
#'@param method Different options to estimate the clusters. The most traditional approach is the k-means analysis (default=’kmeans’)
#'but the function also support the different methods included in the hclust . These methods are:
#'"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' For more details about these methods see the hclust function documentation included in the stats package.
#'@param iter.max Parameter to select the maximum number of iterations allowed (Only if method='kmeans' is selected).
#'@param nstart Parameter for the cluster analysis determining how many random sets to choose (Only if method='kmeans' is selected).
#'@param ncores The number of multicore threads to use for parallel computation. 
#'@return A list with two elements \code{$data} (a 's2dv_cube' object containing the composites cluster=1,..,K for case (*1)
# or only k=1 for any specific cluster, i.e., case (*2)) and \code{$statistics} that includes
#'         \code{$pvalue} (array with the same structure as \code{$data} containing the pvalue of the composites obtained through a t-test that accounts for the serial dependence.),
#'          \code{cluster} (A matrix or vector with integers (from 1:k) indicating the cluster to which each time step is allocated.),
#'         \code{persistence} (Percentage of days in a month/season before a cluster is replaced for a new one (only if method=’kmeans’ has been selected.)),
#'         \code{frequency} (Percentage of days in a month/season belonging to each cluster (only if method=’kmeans’ has been selected).),
#'@importFrom s2dverification EOF
#'@import multiApply
#'@examples
#'\dontrun{
#'res1 <- CST_WeatherRegimes(data = lonlat_data$obs, EOFs = FALSE, ncenters = 4)
#'res2 <- CST_WeatherRegimes(data = lonlat_data$obs, EOFs = TRUE, ncenters = 3)
#'}
#'@export
#'
CST_WeatherRegimes <- function(data, ncenters = NULL,
                              EOFs = TRUE, neofs = 30,
                              varThreshold = NULL,
                              method = "kmeans",
                              iter.max = 100, nstart = 30,
                              ncores = NULL)  {
  if (!inherits(data, 's2dv_cube')) {
    stop("Parameter 'data' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load.")
  }
  if ('lon' %in% names(data)){
    lon <- data$lon
  }else {
    lon <- NULL
  }
  result <- WeatherRegime(data$data,ncenters = ncenters,
                          EOFs = EOFs, neofs = neofs,
                          varThreshold = varThreshold, lon = lon,
                          lat = data$lat, method = method,
                          iter.max=iter.max, nstart = nstart,
                          ncores = ncores)
  data$data <- result$composite
  data$statistics <- result[-1]
  return(data)
}

#' @rdname WeatherRegimes
#' @title Function for Calculating the Cluster analysis
#'
#' @author Verónica Torralba - BSC, \email{veronica.torralba@bsc.es}
#'
#' @description This function computes the weather regimes from a cluster analysis.
#'It can be applied over the dataset with dimensions
#'c(year/month, month/day, lon, lat), or by using PCs obtained from the application of the
#'EOFs analysis to filter the dataset.
#'The cluster analysis can be performed with the traditional k-means or those methods
#'included in the hclust (stats package).
#'
#'@references Cortesi, N., V., Torralba, N., González-Reviriego, A., Soret, and F.J., Doblas-Reyes (2019).
#' Characterization of European wind speed variability using weather regimes. Climate Dynamics,53,
#' 4961–4976, doi:10.1007/s00382-019-04839-5.
#'@references Torralba, V. (2019) Seasonal climate prediction for the wind energy sector: methods and tools
#' for the development of a climate service. Thesis. Available online: \url{https://eprints.ucm.es/56841/}
#'
#'@param data an array containing anomalies with named dimensions with at least start date 'sdate', forecast time 'ftime', latitude 'lat' and longitude 'lon'.
#'@param ncenters Number of clusters to be calculated with the clustering function.
#'@param EOFs Whether to compute the EOFs (default = 'TRUE') or not (FALSE) to filter the data.
#'@param neofs number of modes to be kept only if EOFs = TRUE has been selected. (default = 30).
#'@param varThreshold Value with the percentage of variance to be explained by the PCs.
#' Only sufficient PCs to explain this much variance will be used in the clustering. 
#'@param lon Vector of longitudes.
#'@param lat Vector of latitudes.
#'@param method Different options to estimate the clusters. The most traditional approach is the k-means analysis (default=’kmeans’)
#'but the function also support the different methods included in the hclust . These methods are:
#'"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' For more details about these methods see the hclust function documentation included in the stats package.
#'@param iter.max Parameter to select the maximum number of iterations allowed (Only if method='kmeans' is selected).
#'@param nstart Parameter for the cluster analysis determining how many random sets to choose (Only if method='kmeans' is selected).
#'@param ncores The number of multicore threads to use for parallel computation. 
#'@return A list with elements \code{$composite} (array with at least 3-d ('lat', 'lon', 'cluster') containing the composites k=1,..,K for case (*1)
# or only k=1 for any specific cluster, i.e., case (*2)),
#'         \code{pvalue} (array with at least 3-d ('lat','lon','cluster') with the pvalue of the composites obtained through a t-test that accounts for the serial
# dependence of the data with the same structure as Composite.),
#'         \code{cluster} (A matrix or vector with integers (from 1:k) indicating the cluster to which each time step is allocated.),
#'         \code{persistence} (Percentage of days in a month/season before a cluster is replaced for a new one (only if method=’kmeans’ has been selected.)),
#'         \code{frequency} (Percentage of days in a month/season belonging to each cluster (only if method=’kmeans’ has been selected).),
#'@importFrom s2dverification EOF
#'@import multiApply
#'@examples
#'\dontrun{
#'res <- WeatherRegime(data = lonlat_data$obs$data, lat = lonlat_data$obs$lat,
#'                     EOFs = FALSE, ncenters = 4)
#'}
#'@export

WeatherRegime <- function(data, ncenters = NULL,
                          EOFs = TRUE,neofs = 30,
                          varThreshold = NULL, lon = NULL,
                          lat = NULL, method = "kmeans",
                          iter.max=100, nstart = 30,
                          ncores = NULL) {
  
  if (is.null(names(dim(data)))) {
    stop("Parameter 'data' must be an array with named dimensions.")
  }
  
  if (is.null(lat)) {
    stop("Parameter 'lat' must be specified.")
  }
  
  dimData <- names(dim(data))

  if ('sdate' %in% dimData && 'ftime' %in% dimData) {
    nsdates <- dim(data)['sdate']
    nftimes <- dim(data)['ftime']
    data <- MergeDims(data, 
                      merge_dims = c('ftime', 'sdate'), 
                      rename_dim = 'time')
  } else if ('sdate' %in% dimData | 'ftime' %in% dimData) {
    names(dim(data))[which(dimData == 'sdate' | dimData == 'ftime') ] = 'time' 
  } else {
      if (!('time' %in% dimData)) {
         stop("Parameter 'data' must have temporal dimensions.")
      }
  }
  

  output <- Apply(data = list(data),
                  target_dims = c('time','lat','lon'),
                  fun = .WeatherRegime,
                  EOFs = EOFs,  neofs = neofs,
                  varThreshold = varThreshold,
                  lon = lon, lat = lat,
                  ncenters = ncenters,
                  method = method,
                  ncores = ncores)
  
  if (method == 'kmeans' && 'sdate' %in% dimData && 'ftime' %in% dimData) {

    # The frequency and the persistency are computed as they are useful 
    # parameters in the cluster analysis
    extra_output <-  Apply(data = output$cluster,
                           target_dims = 'time',
                           fun = .freqPer,
                           nsdates = nsdates,
                           nftimes = nftimes ,
                           ncenters = ncenters)
  
     output$cluster <-  t(array(output$cluster, dim = c(nftimes, nsdates)))
     names(dim(output$cluster)) <- c('sdate', 'ftime')

     output <- list(composite = output$composite,
                   pvalue = output$pvalue,
                   cluster = output$cluster,
                   frequency = extra_output$frequency, 
                   persistence = extra_output$persistence)
  }
  return(output)
}

.WeatherRegime <- function(data, ncenters = NULL, EOFs = TRUE,neofs = 30,
                           varThreshold = NULL, lon = NULL,
                           lat = NULL, method = "kmeans",
                           iter.max=100, nstart = 30) {
  
  if (is.null(names(dim(data)))) {
    stop("Parameter 'data' must be an array with 'time', 'lat' and 'lon' dimensions.")
  }
  
  if (!is.null(lat) && dim(data)['lat'] != length(lat)) {
    stop("The length of the paramter 'lat' does not match with the ['lat'] dimension of 
         the parameter 'data'.")
  }
  if (is.null(ncenters)) {
    stop("Parameter 'ncenters' must be specified.")
  }
  if (EOFs  == TRUE && is.null(lon)) {
    stop("Parameter 'lon' must be specified.")
  }
  if (is.null(lat)) {
    stop("Parameter 'lat' must be specified.")
  }
  
  nlon <- dim(data)['lat']
  nlat <- dim(data)['lon']
  
  if (any(is.na(data))){
    nas_test <- MergeDims(data, merge_dims = c('lat','lon'),
                          rename_dim = 'space',na.rm = TRUE)
    if (dim(nas_test)['space']== c(nlat*nlon)){
      stop("Parameter 'data' contains NAs in the 'time' dimensions.")
    }
  }

  if (EOFs  == TRUE) {
    if (is.null(varThreshold)) {
      dataPC <- EOF(data,
                    lat = as.vector(lat),
                    lon = as.vector(lon),
                    neofs = neofs)
      cluster_input <- dataPC$PC
    } else {
      dataPC <- EOF(data,
                    lat = as.vector(lat),
                    lon = as.vector(lon),
                    neofs = neofs)
      minPC <- 
        head(as.numeric(which(cumsum(dataPC$var) > varThreshold)), 1)
      cluster_input <- dataPC$PC[, 1:minPC]
    }
  } else {
  
    dataW <- aperm(Apply(data, target_dims = 'lat', 
                         function (x, la) {
                           x * cos(la * pi / 180)},
                         la = lat)[[1]], c(2, 1, 3))
    
    cluster_input <- MergeDims(dataW, merge_dims = c('lat','lon'),
                               rename_dim = 'space',na.rm = TRUE)

  }
  
  if (method == "kmeans") {
   
    clust <- kmeans(
      cluster_input,
      centers = ncenters,
      iter.max = iter.max,
      nstart = nstart,
      trace = FALSE)
    
    result <- array(0, c(ncenters, nlat, nlon))
    # the order of the data dimensions is changed ('lat','lon','time')
    result <- Composite(aperm(data,c(2, 3, 1)), clust$cluster)
    
  } else {
    result <- hclust(dist(cluster_input), method = method)
    clusterCut <- cutree(result, ncenters)
    result <- Composite(aperm(data, c(2, 3, 1)), clusterCut)
  }
  result <- lapply(1:length(result),
                   function (n) {
                     names(dim(result[[n]])) <- c("lat", "lon", "cluster")
                     return (result[[n]])
                   })
  
  names(result) <- c('composite','pvalue')

  if (method == "kmeans") {
    clust <- as.array(clust$cluster)
    names(dim(clust)) <- 'time'
    return(list(
      composite = result$composite,
      pvalue = result$pvalue,
      cluster = clust))
  } else {
    clust <- as.array(clusterCut)
    names(dim(clust)) <- 'time'
    return(list(
      composite = result$composite,
      pvalue = result$pvalue,
      cluster = clust))
  }
}

.freqPer<- function (clust, nsdates, nftimes, ncenters){
  frequency <- persistence <- matrix(NA, nsdates, ncenters)
  x <- as.vector(clust)
  for (i in 1:nsdates) {
    occurences <-rle(x[((i * nftimes) + 1 - nftimes):(i * nftimes)])
    for (j in 1:ncenters) {
      frequency[i, j] <-(sum(occurences$lengths[occurences$values == j]) /  nftimes) * 100
      persistence[i, j] <- mean(occurences$lengths[occurences$values == j])
    }
  }
  return(list(frequency = frequency,
              persistence = persistence))
}
