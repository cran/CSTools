context("Generic tests")
test_that("Sanity checks", {
  expect_error(
    CST_BiasCorrection(exp = 1),
    paste0("Parameter 'exp' and 'obs' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load."))
  
  mod1 <- 1 : (1 * 3 * 4 * 5 * 6 * 7)
  obs1 <- 1 : (1 * 1 * 4 * 5 * 6 * 7)
  dim(mod1) <- c(dataset = 1, member = 3, sdate = 4, ftime = 5, 
                 lat = 6, lon = 7)
  dim(obs1) <- c(dataset = 1, member = 1, sdate = 4, ftime = 5, 
                 lat = 6, lon = 7)
  lon <- seq(0, 30, 5)
  lat <- seq(0, 25, 5)
  exp <- list(data = mod1, lat = lat, lon = lon)
  obs <- list(data = obs1, lat = lat, lon = lon)
  attr(exp, 'class') <- 's2dv_cube'
  attr(obs, 'class') <- 's2dv_cube'

  bc <- CST_BiasCorrection(exp = exp, obs = obs)  
  expect_equal(length(bc), 3)
  expect_equal(dim(bc$data), 
               c(dataset = 1, member = 3, sdate = 4, ftime = 5, 
                 lat = 6, lon = 7))
  expect_equal(bc$lat, lat)
  expect_equal(bc$lon, lon)
  
  expect_error(CST_BiasCorrection(exp = exp, obs = exp),
        paste0("The length of the dimension 'member' in the component 'data' ",
         "of the parameter 'obs' must be equal to 1."))
  
  exp2 <- exp
  exp2$data[1, 2, 1, 1, 1, 1] <- NA
  expect_warning(CST_BiasCorrection(exp = exp2, obs = obs),
    "Parameter 'exp' contains NA values.")
  
  obs2 <- obs
  obs2$data[1, 1, 2, 1, 1, 1] <- NA
  expect_warning(CST_BiasCorrection(exp = exp, obs = obs2),
    "Parameter 'obs' contains NA values.")
  
  expect_warning(CST_BiasCorrection(exp = exp2, obs = obs2),
    "Parameter 'obs' contains NA values", "Parameter 'exp' contains NA values.")
})
