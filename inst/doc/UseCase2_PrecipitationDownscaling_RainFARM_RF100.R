#!/usr/bin/env Rscript
rm(list=ls()); gc();

### Creation date: 9th June 2021
# Author: N. Pérez-Zanón
# Adapted from the original version Terzago et al. 2020
# https://drive.google.com/file/d/1qp2gbtKdBl4XmsyOeaEhFENwpeUuJwkf/view
# Refer CSTools package and manuscript when using it.
# ----------------------------------------
# This code is using divided in 4 steps plus visualization of the results
# Taking advantage of 
# STEP 1: Compute Slope values for RainFARM downscaling method
# STEP 2: Apply Quantile Mapping Correction to Seasonal Precipitation Forecast
# STEP 3: Compute Weights for RainFARM downscaling method
# STEP 4: Apply RainFARM downscaling method
# ----------------------------------------
# Note: STEP 3 requires a machine with internet connexion. 
#       In this file, the lines are commented since they have been run and the 
#       result saved on disk, then the result is loaded.
# ----------------------------------------
#
# Load required libraries and setup output directory:
library(CSTools)
library(ClimProjDiags)
library(zeallot)
library(ragg)
dir_output <- '/esarchive/scratch/nperez/CSTools_manuscript/v20210603/' #slash end

# --------------------------------------------
# STEP 1:
# --------------------------------------------
era5 <- list(name = 'era5',
             path = '/esarchive/recon/ecmwf/era5/$STORE_FREQ$_mean/$VAR_NAME$_f1h-r1440x721cds/$VAR_NAME$_$YEAR$$MONTH$.nc')
years <- unlist(lapply(1993:2018, function(x){
                       paste0(x, sprintf("%02d",1:12), '01')}))
era5 <- CST_Load(var = 'prlr',
                 exp = list(era5),
                 sdates = years, nmember = 1,
                 storefreq = "daily", sampleperiod = 1,
                 latmin = 37.5, latmax = 53.25, lonmin = 2.5, lonmax = 18.25,
                 output = 'lonlat', nprocs = 1)
era5$data <- era5$data * 24 * 3600 * 1000 # with or without this line -> same result
era5 <- CST_SplitDim(era5, split_dim = 'sdate', indices = rep(1:12, 26))
slope <- CST_RFSlope(era5, time_dim = c('sdate', 'ftime'), kmin = 5)
save(slope, file = paste0(dir_output, 'Slope.RDS'), version = 2)
slope_plot <- slope

# --------------------------------------------
# STEP 2:
# --------------------------------------------
StartDates <- paste0(1993:2018, '1101')
exp <- list(name = 'ecmwfS5',
            path = "/esarchive/exp/ecmwf/system5c3s/$STORE_FREQ$_mean/$VAR_NAME$_s0-24h/$VAR_NAME$_$START_DATE$.nc")
obs <- list(name = 'era5',
            path = '/esarchive/recon/ecmwf/era5/$STORE_FREQ$_mean/$VAR_NAME$_f1h-r1440x721cds/$VAR_NAME$_$YEAR$$MONTH$.nc')
#obs <- list(name = 'era5', path = '/esarchive/scratch/nperez/ERA5/daily_total_prlr_f1h-r1440x721cds/$VAR_NAME$_$YEAR$$MONTH$.nc')

c(exp, obs) %<-% CST_Load(var = 'prlr', exp = list(exp), obs = list(obs),
                          sdates = StartDates, nmember = 25,
                          storefreq = "daily", sampleperiod = 1,
                          latmin = 42, latmax = 49, lonmin = 4, lonmax = 11,
                          output = 'lonlat', nprocs = 1)
exp <- CST_SplitDim(exp, split_dim = c('ftime'))
obs <- CST_SplitDim(obs, split_dim = c('ftime'))
exp$data <- exp$data * 24 * 3600 * 1000
obs$data <- obs$data * 24 * 3600 * 1000
exp$data[which(exp$data < 0)] <- 0
exp.qm <- CST_QuantileMapping(exp, obs, method = "QUANT",
                              wet.day = FALSE,
                              sample_dims = c('member', 'sdate', 'ftime'),
                              ncores = 4)
save(exp.qm, file = paste0(dir_output, 'ExpQM.RDS'), version = 2)

# --------------------------------------------
# STEP 3:
# --------------------------------------------
#library(raster); library(s2dv); library(CSTools)
#worldclim <- getData("worldclim", var = "prec", res = 0.5, lon = 5, lat = 45)
#wc_month <- lapply(1:12, FUN = function(x) {
#                   res <- crop(worldclim[[x]],
#                               extent(3.5, 11.5, 41.5, 49.5))
#                   res <- as.array(res)
#                   names(dim(res)) <- c('lat', 'lon', 'month')
#                   return(res)
#                   })
#xy <- xyFromCell(crop(worldclim[[1]], extent(3.5, 11.5, 41.5, 49.5)),
#                 1:length(crop(worldclim[[1]], extent(3.5, 11.5, 41.5, 49.5))))
#lons <- unique(xy[,1])
#lats <- unique(xy[,2])
#wc_month <- unlist(wc_month)
#dim(wc_month) <- c(lat = length(lats), lon = length(lons), monthly = 12)
#wc_month <- Reorder(wc_month, c('lon', 'lat', 'monthly'))
#wc_month <- s2dv_cube(data = wc_month, lon = lons, lat = lats,
#                       Datasets = 'WorldClim')
#weight <- CST_RFWeights(wc_month, lon = 4:11, lat = 49:42, nf = 100)
#save(weight, file = paste0(dir_output, 'weightsRF100.RDS'))

load(paste0(dir_output, 'weightsRF100.RDS'))

# --------------------------------------------
# Visualization
# --------------------------------------------
agg_png(paste0(dir_output, "RF100_WeightsDec.png"),
        width = 1000, height = 1100, units = 'px',res = 144)
PlotEquiMap(weight$data[,,12], lon = weight$lon, lat = weight$lat,
            filled.continents = FALSE, title_scale = 1,
            intylat = 2, intxlon = 2,
            toptitle = 'December Weights RF 100')
dev.off()

# --------------------------------------------
# STEP 4:
# --------------------------------------------

weights <- Subset(weight$data, along = 'monthly', indices = c(11, 12, 1:6))
slope <- Subset(slope, along = 'monthly', indices = c(11, 12, 1:6), 
                drop = 'non-selected')
k = 1 # To create the member ID
for (realizations in 1:10) {
  for (member in 1:25) {
    result <- exp.qm # to store the data
    result$data <- NULL
    for (month in 1:8) {
      data <- exp.qm # to take the correct piece of data
      data$data <- data$data[1, member, , , , , month]
      fs <- CST_RainFARM(data, nf = 100,
                        weights = weights, slope = slope[month],
                        kmin = 1, nens = 1, verbose = TRUE,
                        nprocs = 8,
                        drop_realization = TRUE)
      result$data <- abind::abind(result$data, fs$data, along = 5)
      if (month == 2 & member == 1 & realization == 1) {
        # ----------------------------
        # Visualization:
        # ----------------------------
        agg_png(paste0(dir_output, "RF100_Down_11dec.png"),
                width = 1000, height = 1100, units = 'px',res = 144)
        PlotEquiMap(fs$data[1,11,,],lon = fs$lon, lat = fs$lat,
                    filled.continents = FALSE, bar_limits = c(0,40),
                    intylat = 2, intxlon = 2, title_scale = 1,
                    triangle_ends = c(TRUE, FALSE),
                    toptitle = 'Downsacaled RF 100', units = 'precipitation (mm)')
        dev.off()
    }
    result$lon <- fs$lon
    result$lat <- fs$lat
    result <- CST_MergeDims(result, merge_dims = c("ftime", "monthly"),
                         na.rm = TRUE)
    result$Dataset <- paste0('RF100_ECMWFC3S_QM_member_', member, '_real_', 
                            realizations)
    result$Dates[[1]] <- exp$Dates[[1]]
    CST_SaveExp(result, destination = dir_output, 
                extra_string = paste0('member', k))
    gc() 
    k = k + 1
    rm(list= list('fs', 'result', 'data'))
  }
} 


