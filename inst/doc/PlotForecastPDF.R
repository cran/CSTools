## ----warning=FALSE,message=FALSE,error=FALSE-----------------------------
library(CSTools)

## ----fig.show = 'hide',warning=F-----------------------------------------
fcst <- data.frame(fcst1 = rnorm(mean = 25, sd = 3, n = 30), 
	fcst2 = rnorm(mean = 23, sd = 4.5, n = 30))
PlotForecastPDF(fcst, tercile.limits = c(20, 26))

## ----fig.show = 'hide',warning=F-----------------------------------------
fcst <- data.frame(fcst1 = rnorm(mean = 25, sd = 3, n = 30),
	fcst2 = rnorm(mean = 23, sd = 4.5, n = 30))
PlotForecastPDF(fcst, tercile.limits = c(20, 26), var.name = "Temperature (ºC)",
	title = "Forecasts valid on 2019-01-01 at Sunny Hills", 
	fcst.names = c("model a", "model b"))

## ----fig.show = 'hide',warning=F-----------------------------------------
fcst <- data.frame(fcst1 = rnorm(mean = 25, sd = 3, n = 30), 
	fcst2 = rnorm(mean = 28, sd = 4.5, n = 30), fcst3 = rnorm(mean = 17, sd = 3, n = 30))
PlotForecastPDF(fcst, tercile.limits = rbind(c(20, 26), c(22, 28), c(15, 22)), 
	var.name = "Temperature (ºC)", title = "Forecasts at Sunny Hills", 
	fcst.names = c("January", "February", "March"), obs = c(21, 24, 17), 
	extreme.limits = rbind(c(18, 28), c(20, 30), c(12, 24)))

## ----fig.show = 'hide',warning=F-----------------------------------------
fcst <- array(cbind(cbind(rnorm(mean = 25, sd = 3, n = 30),
	rnorm(mean = 23, sd = 4.5, n = 30)), rnorm(mean = 17, sd = 3, n = 30)), 
	dim = c(members = 30, 3))
PlotForecastPDF(fcst, tercile.limits = rbind(c(20, 26), c(22, 28), c(15, 22)), 
	var.name = "Temperature (ºC)", title = "Forecasts at Sunny Hills", 
	fcst.names = c("January", "February", "March"), obs = c(21, 24, 17), 
	extreme.limits = rbind(c(18, 28), c(20, 30), c(12, 24)))

## ----fig.show = 'hide',warning=F-----------------------------------------
fcst <- data.frame(fcst1 = lonlat_data$exp$data[1,,1,1,1,1] - 273.15,
                   fcst2 = lonlat_data$exp$data[1,,1,2,1,1] - 273.15)
PlotForecastPDF(fcst, tercile.limits = c(5, 7), extreme.limits = c(4, 8), 
  var.name = "Temperature (ºC)",
  title = "Forecasts valid on 2000-11 at sample mediterranean region", 
  fcst.names = c("November", "December"))

