#!/usr/bin/env Rscript
rm(list=ls()); gc();

### Creation date: 3rd June 2021
# Author: N. Pérez-Zanón
# Refer CSTools package and manuscript when using it.
# ----------------------------------------
# Wind speed bias adjustment for assessment of an extreme event
# The System5-ECMWF downloaded from C3S seasonal forecasts initialized 
# in December 2017, January 2018 and February 2018 
# This code includes the bias adjustent and the results visualization 
# ----------------------------------------

#library(CSTools)
library(s2dv)
library(ragg)
library(multiApply)
output_dir <- "/esarchive/scratch/nperez/CSTools_manuscript/v20210603/"


exp_path <- list(name = "ECMWFS5", 
                 path = "/esarchive/exp/ecmwf/system5c3s/$STORE_FREQ$_mean/$VAR_NAME$_f6h/$VAR_NAME$_$START_DATE$.nc")
obs_path <- list(name = "ERA5",
                 path = "/esarchive/recon/ecmwf/era5/$STORE_FREQ$_mean/$VAR_NAME$_f1h/$VAR_NAME$_$YEAR$$MONTH$.nc")
#source("/esarchive/scratch/nperez/git/cstools/R/CST_BiasCorrection.R")
#library(multiApply)
  # Target months March (3)
  # Assess forecast from 1 to 3 months in advance
months_in_advance <- c('02', '01', '12')
wind_fsct_BC <- list()
wind_ref_terciles <- NULL
wind_ref_extremes <- NULL
wind_thres_latlon <- NULL
  # Observations March 2018
  wind_obs <- CSTools::CST_Load(var = 'windagl100', obs = list(obs_path),
                       sdates = '20180301', nmember = 1,
                        leadtimemin = 1, leadtimemax = 1,
                        storefreq = "monthly", sampleperiod = 1,
                        latmin = 36, latmax = 44, lonmin = -10, lonmax = 4,
#                       latmin = 42, latmax = 44, lonmin = -10, lonmax = 1,
                        output = 'lonlat', nprocs = 1, grid = 'r360x181')
  # For each month in advance:
for (mm in 1:3) {
  print("Three initializations:")
  print(mm)
  print(paste('Initialization', months_in_advance[mm]))
  # Generate the start dates of hindcast period
  year <- ifelse(mm == 3, 2017, 2018)
  print(paste('Hindcast period until', year - 1))
  hcst_sdates <- paste0(1993:(year - 1), months_in_advance[mm], '01') 
  wind_hcst <- CSTools::CST_Load(var = 'sfcWind', exp = list(exp_path),
                        sdates = hcst_sdates, nmember = 25, 
                        leadtimemin = mm + 1, leadtimemax = mm + 1,
                        storefreq = "monthly", sampleperiod = 1,
                        latmin = 36, latmax = 44, lonmin = -10, lonmax = 4,
                       # latmin = 42, latmax = 44, lonmin = -10, lonmax = 1,
                        output = 'lonlat', nprocs = 1)
  str(wind_hcst$Dates)
  dim(wind_hcst$data)
  fcst_sdates <- paste0(year, months_in_advance[mm], '01')
  wind_fcst <- CSTools::CST_Load(var = 'sfcWind', exp = list(exp_path),
                        sdates = fcst_sdates, nmember = 25,
                        leadtimemin = mm + 1, leadtimemax = mm + 1,
                        storefreq = "monthly", sampleperiod = 1,
                        latmin = 36, latmax = 44, lonmin = -10, lonmax = 4,
                      # latmin = 42, latmax = 44, lonmin = -10, lonmax = 1,

                        output = 'lonlat', nprocs = 1)
  str(wind_fcst$Dates)
  dim(wind_fcst$data)
  
  wind_ref <- CSTools::CST_Load(var = 'windagl100', obs = list(obs_path),
                       sdates = hcst_sdates, nmember = 1,
                        leadtimemin = mm + 1, leadtimemax = mm + 1,
                        storefreq = "monthly", sampleperiod = 1,
                       latmin = 36, latmax = 44, lonmin = -10, lonmax = 4,
#                       latmin = 42, latmax = 44, lonmin = -10, lonmax = 1,
                        output = 'lonlat', nprocs = 1,
                        grid = 'r360x181')
  str(wind_ref$Dates)
  dim(wind_ref$data)
  print(wind_ref$Dates$start)
                                    
  wind_ref_terciles <- rbind(wind_ref_terciles,
                             quantile(MeanDims(wind_ref$data, c('lat', 'lon')), c(0.3, 0.6)))
  wind_ref_extremes <- rbind(wind_ref_extremes,
                             quantile(MeanDims(wind_ref$data, c('lat', 'lon')), c(0.1, 0.9))) 
  source("/esarchive/scratch/nperez/git/cstools/R/CST_BiasCorrection.R")
  wind_fsct <- CST_BiasCorrection(exp = wind_hcst,
                                  obs = wind_ref, 
                                  exp_cor = wind_fcst)
  wind_fsct_BC[[mm]] <- wind_fsct
  # -----------------------------------------------
  #  PLOTTING: 
  # PlotMostlikely
  thres <- drop(Apply(list(wind_ref$data), target_dims = 'sdate', fun = function(x) {
                 quantile(x, c(1/3, 2/3))}, output_dims = 'probs')$output)
  PB <- Apply(list(wind_fsct$data, thres), target_dims = list('member', 'probs'),
              fun = function(x, z) {
                 result <- unlist(lapply(x, function(y) {
                   if (y <= z[1]) { 
                     res <- c(1, 0, 0)
                   } else if (y <= z[2]) {
                     res <- c(0, 1, 0)
                   } else {
                     res <- c(0, 0, 1)
                   }
                 return(res)}))
                 dim(result) <- c(bin = 3, member = 25)
                 return(result)
              })$output1
  Mean_PB <- drop(MeanDims(PB, 'member'))

  observed_tercile <- Apply(list(wind_obs$data, thres, Mean_PB),
             target_dims = list(NULL, 'probs', 'bin'),
             fun = function(x, y, z) {
                 val <- which.max(z)
                 if (val == 3) {
                   dot <- ifelse(x >= y[2], 1, 0)
                 } else if (val == 2) {
                   dot <- ifelse(x < y[2] && x >= y[1], 1, 0)
                 } else if (val == 1) {
                   dot <- ifelse(x < y[1], 1, 0)
                 } else {
                   stop('what')
                 }
                return(dot)
               })$output1

  filtered_obs_terciles <- Apply(list(Mean_PB, observed_tercile), 
             target_dims = list('bin', NULL), 
             fun = function(x,y) {
                if (sum(duplicated(x)) == 1) {
                  y <- 0
                } 
                return(y) })$output1

  wind_thres_latlon <- abind::abind(wind_thres_latlon, thres, along = 4)
  source("/esarchive/scratch/nperez/git/cstools/R/PlotCombinedMap.R")
  source("/esarchive/scratch/nperez/git/cstools/R/PlotMostLikelyQuantileMap.R")
  agg_png(paste0(output_dir, "Wind_MostLikely_", mm, "_obstercile.png"),
          width = 1050, height = 1000, units = 'px', res = 144) 
  PlotMostLikelyQuantileMap(probs = Mean_PB, lon = wind_fsct$lon,
                          lat = wind_fsct$lat, sizetit = 1.5,
                          intylat = 2, intxlon = 2, 
                          coast_width = 1.5, legend_scale = 0.8, 
                          cat_dim = 'bin', dot_size = 2.5,
                          axes_label_scale = 1.6, degree_sym = TRUE, 
                          dots = filtered_obs_terciles[,,1,1,1,1],
                          toptitle = c(paste0("Initialized on ", 
                              month.name[as.numeric(months_in_advance[mm])]))) 
  dev.off() 
}

visual <- data.frame(dec = as.vector(MeanDims(wind_fsct_BC[[3]]$data, c('lat', 'lon'))),
                     jan = as.vector(MeanDims(wind_fsct_BC[[2]]$data, c('lat', 'lon'))),
                     feb = as.vector(MeanDims(wind_fsct_BC[[1]]$data, c('lat', 'lon'))))

  wind_obs_areave <- CSTools::CST_Load(var = 'windagl100', obs = list(obs_path),
                       sdates = '20180301', nmember = 1,
                        leadtimemin = 1, leadtimemax = 1,
                        storefreq = "monthly", sampleperiod = 1,
                        latmin = 36, latmax = 44, lonmin = -10, lonmax = 4,
#                       latmin = 42, latmax = 44, lonmin = -10, lonmax = 1,
                        output = 'areave', nprocs = 1)

print("IS DATA LOADED")

print("Wait")
agg_png(paste0(output_dir, "Wind_PlotForecast_IP.png"),
    width = 1000, height = 1000, units = 'px',res = 150)
CSTools::PlotForecastPDF(visual, tercile.limits = wind_ref_terciles,
                extreme.limits = wind_ref_extremes,
                var.name = "Wind Speed 100 m (m/s)",
                title = "Bias Corrected forecasts at IP", 
                fcst.names = c("December", "January", "February"), 
                obs = as.vector(wind_obs_areave$data))
dev.off()

# Plotting observed terciles:
names(dim(wind_thres_latlon)) <- c('thres', 'lat', 'lon', 'sdate')
wind_thres_latlon <- ClimProjDiags::Subset(wind_thres_latlon, indices = 1, along = 'sdate')
wind_obs_obstercile <- Apply(list(wind_obs$data, wind_thres_latlon), 
                             target_dims = list(NULL, 'thres'),
                             fun = function(x, tercile) {
                               if (x <= tercile[1]) {
                                 res <- 1
                               } else if (x > tercile[2]) {
                                 res <- 3
                               } else {
                                 res <- 2
                               }
                               return(res)
                             })$output1
                              
agg_png(paste0(output_dir, "MostLikely_Observed_obstercile.png"),
        width = 1000, height = 1000, units = 'px', res = 144)

s2dv::PlotEquiMap(wind_obs_obstercile, 
                  lon = wind_obs$lon, lat = wind_obs$lat,
                  brks = c(0,1,2,3),
                  cols = c("#6BAED6FF", "#FFEDA0FF", "#FC4E2AFF"),
                          intylat = 2, intxlon = 2,
                          coast_width = 1.5, filled.continents = FALSE,
                          toptitle = c("Observed terciles March 2018"))
dev.off()

# All gridpoints are above normal observed tercile.

print("DONE")
 
