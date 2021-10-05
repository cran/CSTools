#!/usr/bin/env Rscript
rm(list=ls()); gc();

### Creation date: 3rd June 2021
# Author: N. Pérez-Zanón
# Adapted from the original version Terzago et al. 2020
# https://drive.google.com/file/d/1qp2gbtKdBl4XmsyOeaEhFENwpeUuJwkf/view
# Refer CSTools package and manuscript when using it.
# ----------------------------------------
# This code is divided in 4 steps plus visualization of the results
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
save(exp, file = paste0(dir_output, 'Exp.RDS'), version = 2)
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
#weight <- CST_RFWeights(wc_month, lon = 4:11, lat = 49:42, nf = 4)
#save(weight, file = paste0(dir_output, 'weightsRF4.RDS'))

load(paste0(dir_output, 'weightsRF4.RDS'))

# --------------------------------------------
# STEP 4:
# --------------------------------------------
Rprof()
weights <- Subset(weight$data, along = 'monthly', indices = c(11, 12, 1:6))
slope <- Subset(slope, along = 'monthly', indices = c(11, 12, 1:6), 
                drop = 'non-selected')
fs <- CST_RainFARM(exp.qm, nf = 4,
                   weights = weights, slope = slope,
                   kmin = 1, nens = 10, verbose = TRUE,
                   time_dim = c("member", "sdate", "ftime"), nprocs = 8,
                   drop_realization = TRUE)
newfs <- CST_MergeDims(fs, merge_dims = c("ftime", "monthly"),
                         na.rm = TRUE)

newfs$Dates[[1]] <- exp$Dates[[1]]
CST_SaveExp(newfs, destination = paste0(dir_output, 'RF4/'))

Rprof(NULL)
profile.info <- summaryRprof(paste0(dir_output, "Rprof.out"))
# --------------------------------------------
# Visualization
# --------------------------------------------
library(s2dv)
agg_png(paste0(dir_output, "EXP_11dec.png"),
        width = 800, height = 900, units = 'px',res = 144)
PlotEquiMap(exp$data[1,1,1,11,,,2],lon = exp$lon, lat = exp$lat,
            filled.continents = FALSE, bar_limits = c(0,40),
           intylat = 2, intxlon = 2, title_scale = 0.8,
           bar_label_scale = 1.3, axes_label_scale = 1.2,
           triangle_ends = c(TRUE, FALSE), degree_sym = TRUE, units_scale = 1.5,
           toptitle = 'SEAS5', units = 'precipitation (mm)')
dev.off()
agg_png(paste0(dir_output, "EXPQM_11dec.png"),
        width = 800, height = 900, units = 'px',res = 144) 
PlotEquiMap(exp.qm$data[1,1,1,11,,,2],lon = exp$lon, lat = exp$lat,
            filled.continents = FALSE, bar_limits = c(0,40),
           intylat = 2, intxlon = 2, title_scale = 0.8,
           bar_label_scale = 1.3, axes_label_scale = 1.2,
           triangle_ends = c(TRUE, FALSE), degree_sym = TRUE, units_scale = 1.5,
           toptitle = 'Bias Adjusted', units = 'precipitation (mm)')
dev.off()
agg_png(paste0(dir_output, "RF4_Down_11dec.png"),
        width = 800, height = 900, units = 'px',res = 144)
PlotEquiMap(fs$data[1,1,1,11,,,2],lon = fs$lon, lat = fs$lat,
            filled.continents = FALSE, bar_limits = c(0,40), 
           intylat = 2, intxlon = 2, title_scale = 0.8,
           bar_label_scale = 1.3, axes_label_scale = 1.2,
           triangle_ends = c(TRUE, TRUE), degree_sym = TRUE, units_scale = 1.5,
           toptitle = 'Downscaled nf 4', units = 'precipitation (mm)')
dev.off()
agg_png(paste0(dir_output, "RF4_WeightsDec.png"),
        width = 800, height = 900, units = 'px',res = 144)
PlotEquiMap(weight$data[,,12], lon = weight$lon, lat = weight$lat,
            filled.continents = FALSE, title_scale = 0.8,
           bar_label_scale = 1.3, axes_label_scale = 1.2,
            intylat = 2, intxlon = 2, degree_sym = TRUE,
            toptitle = 'December Weights nf 4')
dev.off()
agg_png(paste0(dir_output, "Slope.png"),
        width = 700, height = 700, units = 'px',res = 144)
plot(1:12, slope_plot, type = 'b', main = 'Slope', pch = 16, xlab = 'month',
      ylab = 'Slope', bty = 'n', cex.main = 1.5, cex.lab = 1.3, cex = 1.5)
lines(12, slope_plot[12], type = 'p', col = 'red', pch = 16, cex = 1.3)
dev.off()

# Plot ForecastPDF
library(abind)
#obsts <- MeanDims(obs$data, c('lat', 'lon'), na.rm = T)
#print(quantile(obsts, c(0.1, 0.3, 0.6, 0.9), na.rm = T))
#expts <- MeanDims(exp$data, c('lat', 'lon'), na.rm = T)
#exp.qmts <- MeanDims(exp.qm$data,  c('lat', 'lon'), na.rm = T)
print("Quantiles gridpoint")
print(quantile(as.vector(exp.qm$data[1,,,11,5,4,2]), c(0.1,0.3,0.6,0.9)))
data <- rbind(exp$data[1, ,1,11,5,4,2], exp.qm$data[1,,1,11,5,4,2])
empty <- array(NA, c(dataset = 2, members = 225))
names(dim(data)) <- names(dim(empty))
data <- abind(data, empty, along = 2)
names(dim(data)) <- names(dim(empty))
data <- abind(data, fs$data[1,,1,11,15,18 ,2], along = 1)
names(dim(data)) <- names(dim(empty))
print(dim(data))
print(quantile(obs$data[1,,,11,5,4,2]))
agg_png(paste0(dir_output, "/FiguresPDF_RF4_11December_GridPoint.png"),
        width = 1000, height = 1000, units = 'px', res = 144)
PlotForecastPDF(data, tercile.limits = c(0.02, 1.17), 
                extreme.limits = c(0.0155555, 7.75), color.set = 'hydro', 
                var.name = "Precipitation (mm)", add.ensmemb = 'no',
               title = "Forecasts issued on Nov 1993 for 11th December 1993",
               fcst.names = c("SEAS5", "Bias Adjusted", "Downscaled nf 4"))
dev.off()

obsts <- MeanDims(obs$data, c('lat', 'lon'), na.rm = T)
print(quantile(obsts, c(0.1, 0.3, 0.6, 0.9), na.rm = T))
expts <- MeanDims(exp$data, c('lat', 'lon'), na.rm = T)
exp.qmts <- MeanDims(exp.qm$data,  c('lat', 'lon'), na.rm = T)
empty <- array(NA, c(dataset = 2, member = 225, sdate = 26, ftime = 31, monthly = 8))
data <- abind(expts, exp.qmts, along = 1)
names(dim(data)) <- names(dim(expts))
data <- abind(data, empty, along = 2)
names(dim(data)) <- names(dim(expts))
fsts <- MeanDims(fs$data, c('lat', 'lon'), na.rm = T)
data <- abind(data, fsts, along = 1)
names(dim(data)) <- c('dataset', 'members', 'sdate', 'ftime', 'monthly')
agg_png(paste0(dir_output, "/FiguresPDF_RF4_11December_Areave.png"),
        width = 1400, height = 800, units = 'px', res = 144)
PlotForecastPDF(data[,,1,11,2], tercile.limits = c(0.67, 2.5), 
                extreme.limits = c(0.09, 7.3), color.set = 'hydro', 
                var.name = "Precipitation (mm)", 
               title = "Forecasts issued on Nov 1993 for 11th December 1993",
               fcst.names = c("SEAS5", "Bias Adjusted", "Downscaled nf 4"))
dev.off()






