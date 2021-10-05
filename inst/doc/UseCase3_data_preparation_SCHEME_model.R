# Author: Bert Van Schaeybroeck
# Use Case 3: Seasonal forecasts for a river flow
# ----------------------------------------------- 
rm(list = ls())
library(CSTools)
library(s2dverification)
library(CSTools)
library(ClimProjDiags)

#SETUP PARAMETERS (TO FIX BEFORE RUNNING SCRIPT):
#------------------------------------------------
var.to.use <- "prlr" #which variable to correct (prlr=precipitation, tasmin, tasmax, tas)

init.yr <- 1993 #initial year (for ECMWF Sys5 = 1993)
end.yr <- 2019 #end year (for ECMWF Sys5 = 2019)
amt.ftime <- 214
n.cores.to.use <- 20 
use.chirps <- T
eval.method.to.use <- "leave-one-out"
domain.high.res <- "greece_high_res"
domain.low.res <- "greece_low_res"
sdate.mon.to.use <- 5 #Month of startdate for ECMWF Sys5 possibilities are 5 (May) or 11 (November)
sdate.day.to.use <- 1 #Day of startdate for ECMWF Sys5 only possibility is 1
make.plot.msk <- F #mask to indicate if figures need to be made based on the output of the Analog function.

#LOCAL PARAMETERS (to be adjusted for each working system)
#---------------------------------------------------------
dir.rdata <- "/mnt/HDS_URCLIM/URCLIM/bertvs/medscope/data/greece_rdata/"

#dir.scripts <- "/mnt/netapp/home/bertvs/ARCHIVE_bertvs/R/medscope/D3.1/"
dir.scratch <- "/scratch/bertvs/"
#Base path for C3S forecast (experiment) dataset:
dir.c3s <- "/mnt/HDS_URCLIM/URCLIM/bertvs/medscope/data/"
#Base path for ERA5 reference or obs dataset:
dir.era5 <- "/mnt/HDS_BREGILABEPOC/BREGILABEPOC/era5/europe/daily/per_mon/" 
#Base path for CHIRPS (reference or obs) rainfall dataset:
dir.chirps <- "/mnt/HDS_MEDYCLIM/MEDYCLIM/PREDANTAR/climate_data/obs/chirps/"
dir.chirps.low.res <- paste0(dir.chirps, "/", domain.low.res, "/per_mon/")
dir.chirps.high.res <- paste0(dir.chirps, "/", domain.high.res, "/per_mon/")
	

#AUXILIARY FUNCTIONS
#-------------------
set.msk <- function(x, msk, const){
	x[msk] = const
	return(x)
}

#FIXED PARAMETERS:
#-----------------
greece.coor.vec <- c(
	lonmin = 18.975, 
	lonmax = 24.025, 
	latmin = 37.975, 
	latmax = 43.025)
greece.coor.lst <- list(
	lon.min = 18.975, 
	lon.max = 24.025, 
	lat.min = 37.975, 
	lat.max = 43.025)
coor.to.use <- greece.coor.lst

europe.coor <- list(
	lon.min = -27, 
	lon.max = 45, 
	lat.min = 33, 
	lat.max = 73.5)

#Large-scale pressure field metadata (necessary for analogs)
var.msl <- "mean_sea_level_pressure"
nc.var.name.msl <- "msl"

#Depending on the variable loaded, different datasets and metadata are used
if(var.to.use == "prlr"){ #Precipitation
	var.era5 <- "total_precipitation"
	time.era5 <- "daily"
	nc.var.name.era5 <- "tp"
	var.chirps <- "precip"
	time.chirps <- "daily"
	nc.var.name.chirps <- "precip"
	cal.meth.to.use <- "bias" #method for bias calibration
	chirps.low.res.daily <- list(
		name = "chirps_low_res",
		path = paste0(dir.chirps.low.res, "chirps-v2.0.$YEAR$.days_greece_low_res-$MONTH$.nc"),
			nc_var_name = nc.var.name.chirps)
	chirps.high.res.daily <- list(
		name = "chirps_high_res",
		path = paste0(dir.chirps.high.res, 
			"chirps-v2.0.$YEAR$.days_greece_high_res-$MONTH$.nc"),
		nc_var_name = nc.var.name.chirps)
	#unit conversions
	mul.cor.era5 <- 1000 * 24
	add.cor.era5 <- 0
	mul.cor.chirps <- 1
	add.cor.chirps <- 0
	mul.cor.exp <- 1000 * 3600 * 24
	add.cor.exp <- 0
	if(use.chirps){
		add.cor.obs <- add.cor.chirps
		mul.cor.obs <- mul.cor.chirps
	} else {
		add.cor.obs <- add.cor.era5
		mul.cor.obs <- mul.cor.er5
	}
} else if(var.to.use == "tas"){
	var.era5 <- "2m_temperature"
	time.era5 <- "daily"
	nc.var.name.era5 <- "t2m"
	cal.meth.to.use <- "mse_min"
	
	#unit conversions
	mul.cor.era5 <- 0
	add.cor.era5 <- 0
	mul.cor.exp <- 0
	add.cor.exp <- 0
} else if(var.to.use == "tasmin"){
	var.era5 <- "2m_temperature"
	time.era5 <- "daily_min"
	nc_var_name.era5 <- "t2m"
	cal.meth.to.use <- "mse_min" #method for bias calibration
	
	#unit conversions
	mul.cor.era5 <- 0
	add.cor.era5 <- 0
	mul.cor.exp <- 0
	add.cor.exp <- 0
} else if(var.to.use == "tasmax"){
	var.era5 <- "2m_temperature"
	time.era5 <- "daily_max"
	nc_var_name.era5 <- "t2m"
	cal.meth.to.use <- "mse_min" #method for bias calibration
	
	#unit conversions
	mul.cor.era5 <- 0
	add.cor.era5 <- 0
	mul.cor.exp <- 0
	add.cor.exp <- 0
}


#Experiment path specification:
ecmwf.s5.daily <- list(
	name = "ecmwfS5",
	path = paste0(dir.c3s,
		"C3S/$EXP_NAME$/$STORE_FREQ$/$VAR_NAME$/",
		"$VAR_NAME$_$START_DATE$.nc"))
#Reference or obs path specifications (ERA5 data available over Europe):
era5.daily <- list(name = "era5",
	path = paste0(dir.era5, "era5-", time.era5, 
		"-europe-", var.era5, "-$YEAR$-$MONTH$.nc"),
	nc_var_name = nc.var.name.era5)
#Reference or obs path specifications for pressure field (ERA5 data available over Europe):
msl.era5.daily <- list(name = "msl",
	path = paste0(dir.era5, "era5-", time.era5, 
		"-europe-", var.msl, "-$YEAR$-$MONTH$.nc"),
	nc_var_name = nc.var.name.msl)

#UNIVERSAL PARAMETERS:
#---------------------
amt.mon.per.yr <- 12
amt.day.per.mon <- 31
sdate.day.str <- formatC(sdate.day.to.use, width = 2, flag = "0")
sdate.mon.str <- formatC(sdate.mon.to.use, width = 2, flag = "0")
day.lst <- formatC(seq(1, amt.day.per.mon), width = 2, flag = "0")
yr.lst <- seq(init.yr, end.yr)
amt.yr <- length(yr.lst)
sdate.lst <- paste0(yr.lst, sdate.mon.str, sdate.day.str ) 


#START
#-----

#1. LOAD THE DATA
#----------------

if(use.chirps){
	obs.set.to.use <- chirps.low.res.daily
} else {
	obs.set.to.use <- era5.daily
}

#Load mean sea level pressure field from ERA5 (no need to set the units)

file.to.load <- paste0(dir.rdata, "msl_all.RData")
if(file.exists(file.to.load)){
	load(file.to.load, verbose = T)
} else {
	msl.all <- CST_Load(
		var = "psl", #nc.var.name.msl,
		obs = list(msl.era5.daily),
		exp = list(ecmwf.s5.daily),
		nmember = NULL,
		sdates = sdate.lst,
		lonmin = europe.coor$lon.min, 
		lonmax = europe.coor$lon.max,
		latmin = europe.coor$lat.min, 
		latmax = europe.coor$lat.max,
		output = "lonlat",
		storefreq = "daily",
		nprocs = n.cores.to.use)
	save(file = file.to.load, msl.all)
}

#Data manipulation: first split lead times per month. Then merge all data per month and all sdates.
#This merged dataset will be used to calibrate (per month) and find analogs (per month).
obs.msl.eur.split <- CST_SplitDim(msl.all$obs, split_dim = c("ftime"))
exp.msl.eur.split <- CST_SplitDim(msl.all$exp, split_dim = c("ftime"))
obs.msl.eur.merge <- CST_MergeDims(
	obs.msl.eur.split,
	merge_dims = c("sdate", "ftime"),
	rename_dim = "sdate")
exp.msl.eur.merge <- CST_MergeDims(
	exp.msl.eur.split,
	merge_dims = c("sdate", "ftime"),
	rename_dim = "sdate")

obs.msl.eur.merge.an <- CST_Anomaly(exp = obs.msl.eur.merge, dim_anom = 3)
exp.msl.eur.merge.an <- CST_Anomaly(exp = exp.msl.eur.merge, dim_anom = 3)


#Load observational and forecast set of variable that needs to be calibrated and downscaled:
file.to.load <- paste0(dir.rdata, "data_all.RData")
if(file.exists(file.to.load)){
	load(file.to.load, verbose = T)
} else {
	data.all <- CST_Load(
		var = var.to.use,
		obs = list(obs.set.to.use),
		exp = list(ecmwf.s5.daily),
		nmember = NULL,
		sdates = sdate.lst,
		lonmin = coor.to.use$lon.min, 
		lonmax = coor.to.use$lon.max,
		latmin = coor.to.use$lat.min, 
		latmax = coor.to.use$lat.max,
		output = "lonlat",
		storefreq = "daily",
		nprocs = n.cores.to.use)
	save(file = file.to.load, data.all)
}
#Set the units:
data.all$obs$data <- data.all$obs$data * mul.cor.obs + add.cor.obs
data.all$exp$data <- data.all$exp$data * mul.cor.exp + add.cor.exp

#Data manipulation: first split lead times per month. Then merge all data per month and all sdates.
#This merged dataset will be used to calibrate (per month) and find analogs (per month).
obs.split <- CST_SplitDim(data.all$obs, split_dim = c("ftime"))
exp.split <- CST_SplitDim(data.all$exp, split_dim = c("ftime"))
obs.merge <- CST_MergeDims(
	obs.split,
	merge_dims = c("sdate", "ftime"),
	rename_dim = "sdate")
exp.merge <- CST_MergeDims(
	exp.split,
	merge_dims = c("sdate", "ftime"),
	rename_dim = "sdate")
	
#Calibrate the exp data (per month)
cal.merge <- CST_Calibration(
	exp = exp.merge, 
	obs = obs.merge, 
	cal.method = cal.meth.to.use, 
	eval.method = eval.method.to.use)
cal.merge$data[cal.merge$data < 0] <- 0

#LOAD HIGH RES CHIRPS DATA
file.to.load <- paste0(dir.rdata, "obs_high_res.RData")
if(file.exists(file.to.load)){
	load(file.to.load, verbose = T)
} else {
	obs.high.res <- CST_Load(var = var.to.use, 
		obs = list(chirps.high.res.daily),
		exp = NULL,
		sdates = sdate.lst,
		nmember = 1,
		leadtimemax = amt.ftime,
		sampleperiod = 1,
		lonmin = coor.to.use$lon.min, 
		lonmax = coor.to.use$lon.max,
		latmin = coor.to.use$lat.min, 
		latmax = coor.to.use$lat.max,
		output = "lonlat",
		storefreq = "daily",
		nprocs = n.cores.to.use)
	save(file = file.to.load, obs.high.res)
}
#set the units
obs.high.res$data <- obs.high.res$data * mul.cor.chirps + 
	add.cor.chirps
#split per month
obs.high.res.split <- CST_SplitDim(
	obs.high.res, 
	split_dim = c("ftime"))
#merge lead times and sdates
obs.high.res.merge <- CST_MergeDims(
	obs.high.res.split,
	merge_dims = c("sdate", "ftime"),
	rename_dim = "sdate")
	
#LOAD LOW RES CHIRPS DATA
file.to.load <- paste0(dir.rdata, "obs_low_res.RData")
if(file.exists(file.to.load)){
	load(file.to.load, verbose = T)
} else {
	obs.low.res <- CST_Load(var = var.to.use, 
	obs = list(chirps.low.res.daily),
	exp = NULL,
	sdates = sdate.lst,
	nmember = 1,
	leadtimemax = amt.ftime,
	sampleperiod = 1,
	lonmin = coor.to.use$lon.min, 
	lonmax = coor.to.use$lon.max,
	latmin = coor.to.use$lat.min, 
	latmax = coor.to.use$lat.max,
	output = "lonlat",
	storefreq = "daily",
	nprocs = n.cores.to.use)
	save(file = file.to.load, obs.low.res)
}
#set units
obs.low.res$data <- obs.low.res$data * mul.cor.chirps + 
	add.cor.chirps
#split per month
obs.low.res.split <- CST_SplitDim(
	obs.low.res, 
	split_dim = c("ftime"))
#merge lead times and sdates
obs.low.res.merge <- CST_MergeDims(
	obs.low.res.split,
	merge_dims = c("sdate", "ftime"),
	rename_dim = "sdate")


#2. PROCESS THE DATA
#-------------------

#amount of ensemble members from experiment. For ECMWF Sys5 it is 25:
amt.mbr <- as.numeric(dim(cal.merge$data)["member"])
lon.low.res <- as.vector(cal.merge$lon)
lat.low.res <- as.vector(cal.merge$lat)
lon.high.res <- as.vector(obs.high.res$lon)
lat.high.res <- as.vector(obs.high.res$lat)
lon.eur <- as.vector(obs.msl.eur.merge.an$lon)
lat.eur <- as.vector(obs.msl.eur.merge.an$lat)

#amount of lead times in months. For ECMWF Sys5 it is 7:
amt.lead.mon <- as.numeric(dim(cal.merge$data)["monthly"])
mon.seq.tmp <- seq(sdate.mon.to.use, sdate.mon.to.use + amt.lead.mon - 1)
mon.seq.tmp <- ((mon.seq.tmp - 1) %% amt.mon.per.yr) + 1
lead.mon.lst <- formatC(mon.seq.tmp, width = 2, flag = "0")
#amount of starting days from experiment. For ECMWF Sys5 it is 837:
amt.sdate <- as.numeric(dim(cal.merge$data)["sdate"])

sub.time <- outer( 
	as.vector(t(outer(yr.lst, day.lst, paste, sep="-"))), 
	lead.mon.lst,
	paste, sep = "-")	
#This step is necessary to set the non-existent dates to NA
sub.time <- format(as.Date(sub.time, format("%Y-%d-%m")), "%Y-%m-%d")
dim(sub.time) <- c(sdate = amt.yr * amt.day.per.mon, time = amt.lead.mon)

cal.high.res.merge <- obs.high.res.merge
cal.high.res.merge$data[] <- NA

#Determine spatial points with all obs.high.res.merge (CHIRPS) data equal to NA. These are the points over sea.
is.na.high.res.obs <- apply(
	obs.high.res.merge$data, 
	MARGIN = c(4, 5),
	FUN = function(x){all(is.na(x))})
#Determine spatial points with all obs.low.res.merge (CHIRPS) data equal to NA. These are the points over sea.
is.na.low.res.obs <- apply(
	obs.low.res.merge$data, 
	MARGIN = c(4, 5),
	FUN = function(x){all(is.na(x))})

#Set all calibrated exp data (cal.merge) equal to NA at the sea point.
cal.merge.tmp = Apply(
	data = list(x = cal.merge$data), 
	target_dims = list(x = c("lat", "lon")),
	fun = set.msk,
	msk = is.na.low.res.obs,
	const = 0,
	output_dims = list(c("lat", "lon"))
	)$output1
dex.match <- match(names(dim(cal.merge$data)), names(dim(cal.merge.tmp)))
cal.merge$data <- aperm(cal.merge.tmp, dex.match)
rm(cal.merge.tmp)

#2. PROCESS THE DATA
#-------------------

i.dataset <- 1
i.mbr.obs <- 1
for(i.mbr in seq(1, amt.mbr)){
	for(i.mon in seq(1, amt.lead.mon)){
		for(i.sdate in seq(1, amt.sdate)){
			#i.mbr <- 1
			#i.mon = 1
			#i.sdate = 24
			date.to.use <- sub.time[ i.sdate, i.mon]
			date.an.lst <- sub.time[ , i.mon]
			cat("i.mbr = ", i.mbr, ", i.mon =", i.mon, ", i.sdate = ", 
				i.sdate, "date: ", date.to.use,"\n")
			
			
			#Extract the (calibrated) forecast that you want to downscale:
			exp.low.res.tmp <- exp.merge$data[i.dataset, i.mbr, i.sdate, , , i.mon]
			cal.low.res.tmp <- cal.merge$data[i.dataset, i.mbr, i.sdate, , , i.mon]
			#Extract the large-scale pressure field of that day
			exp.msl.eur.tmp <- exp.msl.eur.merge.an$data[i.dataset, i.mbr, i.sdate, , , i.mon]
			
			#Extract all observations that will be used to find analogs 
			obs.msl.eur.tmp <- obs.msl.eur.merge.an$data[i.dataset, i.mbr.obs, , , , i.mon]#-i.sdate
			obs.low.res.tmp <- obs.low.res.merge$data[i.dataset, i.mbr.obs, , , , i.mon] #-i.sdate
			obs.high.res.tmp <- obs.high.res.merge$data[i.dataset, i.mbr.obs, , , , i.mon] #-i.sdate
			names(dim(obs.high.res.tmp)) <- c("time", "lat", "lon")
			names(dim(obs.low.res.tmp)) <- c("time", "lat", "lon")
			names(dim(obs.msl.eur.tmp)) <- c("time", "lat", "lon")
			if(!is.na(date.to.use) & !all(is.na(cal.low.res.tmp))){
				obs.low.res.tmp[is.na(obs.low.res.tmp)] <- 0
				
				res  <- Analogs(
					expL = exp.msl.eur.tmp, 
					obsL = obs.msl.eur.tmp, 
					time_obsL = date.an.lst,
					obsVar = obs.low.res.tmp,
					expVar = exp.low.res.tmp,
					lonVar = lon.low.res,
					latVar = lat.low.res,
					lonL = lon.eur, 
					latL = lat.eur, 
					region = greece.coor.vec,
					criteria = "Local_dist",
					time_expL = date.to.use,
					excludeTime = date.to.use,
					AnalogsInfo = T, 
					nAnalogs = 1000)
				
				
				if(make.plot.msk){
					corr.date <- as.character(res$dates[1]) #select the date of the most
					corr.dex <- which(date.an.lst == corr.date)
					
					#The following figure shows the uncalibrated raw model field (analogous to Fig. 9a)
					file.fig <- paste0("mbr_", i.mbr, "_mon_", i.mon, 
						"_sdate_", date.to.use, "_exp.low.res.pdf")
					pdf(file = file.fig)
					PlotEquiMap(
						exp.low.res.tmp[ , ], 
						lon = obs.low.res.merge$lon, 
						lat = obs.low.res.merge$lat,
						filled.continents = F, 
						intylat = 2, 
						intxlon = 2, 
						title_scale = 0.7, #bar_limits = c(0, 60),
						units = "precipitation (mm)") 
					dev.off()
					
					#The following figure includes the calibrated model field (analogous to Fig. 9b)
					file.fig <- paste0("mbr_", i.mbr, "_mon_", i.mon, 
						"_sdate_", date.to.use, "_cal.low.res.pdf")
					pdf(file = file.fig)
					PlotEquiMap(
						cal.low.res.tmp, 
						lon = obs.low.res.merge$lon, 
						lat = obs.low.res.merge$lat,
						filled.continents = F,
						intylat = 2, 
						intxlon = 2, 
						title_scale = 0.7, #bar_limits = c(0, 60),
						units = "precipitation (mm)")
					
					#The following figure includes the analog upscaled field (analogous to Fig. 9c)
					file.fig <- paste0("mbr_", i.mbr, "_mon_", i.mon, 
						"_sdate_", date.to.use, "_obs.low.res.pdf")
					pdf(file = file.fig)
					PlotEquiMap(
						obs.low.res.tmp[corr.dex, , ], 
						lon = obs.low.res.merge$lon, 
						lat = obs.low.res.merge$lat,
						filled.continents = F, 
						intylat = 2, 
						intxlon = 2, 
						title_scale = 0.7, #bar_limits = c(0, 60),
						units = "precipitation (mm)") 
					dev.off()
					
					#The following figure includes the analog field (analogous to Fig. 9d)
					file.fig <- paste0("mbr_", i.mbr, "_mon_", i.mon, 
						"_sdate_", date.to.use, "_obs.high.res.pdf")
					pdf(file = file.fig)
					PlotEquiMap(
						obs.high.res.tmp[corr.dex, , ], 
						lon = obs.high.res.merge$lon, 
						lat = obs.high.res.merge$lat,
						filled.continents = F,
						intylat = 2, 
						intxlon = 2, 
						title_scale = 0.7, #bar_limits = c(0, 60),
						units = "precipitation (mm)")
					dev.off()
					
					
				}
			}
		}
	}
}



