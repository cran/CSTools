## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(eval = FALSE)

## ------------------------------------------------------------------------
#  install.packages('CSTools')
#  library(CSTools)

## ------------------------------------------------------------------------
#  library(CSTools)
#  exp <- lonlat_prec

## ------------------------------------------------------------------------
#  dim(exp$data)
#  #dataset  member   sdate   ftime     lat     lon
#  #      1       6       3      31       4      4

## ------------------------------------------------------------------------
#  ilon <- which ( exp$lon %in% 5:12 )
#  ilat <- which ( exp$lat %in% 40:47 )
#  exp$data <- exp$data[ , , , , ilon, ilat, drop=FALSE]
#  names(dim(exp$data)) <- names(dim(lonlat_prec$data))
#  exp$lon <- exp$lon[ilon]
#  exp$lat <- exp$lat[ilat]

## ------------------------------------------------------------------------
#  exp_down <- CST_RainFARM(exp, nf=20, kmin = 1, nens = 3,
#                           time_dim = c("member", "ftime"))
#  
#  dim(exp_down$data)
#  #    dataset      member realization       sdate       ftime         lat        lon
#  #          1           6           3           3          31          80         80
#  str(exp_down$lon)
#  # num [1:80] 5.53 5.58 5.62 5.67 5.72 ...
#  str(exp_down$lat)
#  # num [1:80] 47.5 47.4 47.4 47.3 47.3 ...

## ------------------------------------------------------------------------
#  downscaled <- RainFARM(exp$data, exp$lon, exp$lat,
#                         nf = 20, kmin = 1, nens = 3,
#                         time_dim = c("member", "ftime"))

## ------------------------------------------------------------------------
#  a <- exp$data[1, 1, 1, 17, , ] * 86400 * 1000
#  a[a > 60] <- 60
#  image(exp$lon, rev(exp$lat), t(apply(a, 2, rev)), xlab = "lon", ylab = "lat",
#        col = rev(terrain.colors(20)), zlim = c(0,60))
#  map("world", add = TRUE)
#  title(main = "pr 17/03/2010 original")
#  a <- exp_down$data[1, 1, 1, 1, 17, , ] * 86400 * 1000
#  a[a > 60] <- 60
#  image(exp_down$lon, rev(exp_down$lat), t(apply(a, 2, rev)), xlab = "lon", ylab = "lat",
#        col = rev(terrain.colors(20)), zlim = c(0, 60))
#  map("world", add = TRUE)
#  title(main = "pr 17/03/2010 downscaled")

## ------------------------------------------------------------------------
#  ww <- CST_RFWeights("./worldclim.nc", nf = 20, lon = exp$lon, lat = exp$lat)

## ------------------------------------------------------------------------
#  exp_down_weights <- CST_RainFARM(exp, nf = 20, kmin = 1, nens = 3,
#                                   weights = ww, time_dim = c("member", "ftime"))

## ------------------------------------------------------------------------
#  exp_down1 <- exp_down$data[1, , , 1, , , ]
#  exp_down_weights1 <- exp_down_weights$data[1, , , 1, , , ]
#  dim(exp_down1) <- c(member = 6 * 3 * 31, lat = 80, lon = 80)
#  dim(exp_down_weights1) <- c(member = 6 * 3 * 31, lat = 80, lon = 80)
#  ad <- apply(exp_down1, c(2, 3), mean)
#  adw <- apply(exp_down_weights1, c(2, 3), mean);
#  
#  png("Figures/RainFARM_fig2.png", width = 640, height = 243)
#  par(mfrow = c(1,3))
#  a <- exp_down_weights$data[1, 1, 1, 1, 17, , ] * 86400 * 1000
#  a[a > 60] <- 60
#  image(exp_down$lon, rev(exp_down$lat), t(apply(a, 2, rev)), xlab = "lon",
#        ylab = "lat", col = rev(terrain.colors(20)), zlim = c(0, 60))
#  map("world", add = TRUE)
#  title(main = "pr 17/03/2010 with weights")
#  a <- ad * 86400 * 1000
#  a[a > 5] <- 5
#  image(exp_down$lon, rev(exp_down$lat), t(apply(a, 2, rev)), xlab = "lon",
#        ylab="lat", col = rev(terrain.colors(20)), zlim = c(0, 5))
#  map("world", add = TRUE)
#  title(main = "climatology no weights")
#  a <- adw * 86400 * 1000
#  a[a > 5] <- 5
#  image(exp_down$lon, rev(exp_down$lat), t(apply(a, 2, rev)), xlab = "lon",
#        ylab = "lat", col = rev(terrain.colors(20)), zlim = c(0, 5))
#  map("world", add = TRUE)
#  title(main = "climatology with weights")
#  dev.off()

## ------------------------------------------------------------------------
#  slopes <- CST_RFSlope(exp, time_dim = c("member", "ftime"))
#  dim(slopes)
#  #    dataset   sdate
#  #          1       3
#  slopes
#  #         [,1]     [,2]     [,3]
#  #[1,] 1.532351 1.664028 1.459252

