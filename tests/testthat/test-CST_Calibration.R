context("Generic tests")
test_that("Sanity checks", {
  expect_error(
    CST_Calibration(exp = 1),
    "Parameter 'exp' and 'obs' must be of the class 's2dv_cube', "
  )
  expect_error(
    CST_Calibration(obs = 1),
    c("argument \"exp\" is missing, with no default")
  )
  library(zeallot)  
  c(exp, obs) %<-% lonlat_data
  cal <- CST_Calibration(exp = exp, obs = obs)
  expect_equal(length(cal), 9)
  expect_equal(dim(cal$data), dim(exp$data))
  expect_equal(cal$lat, exp$lat)
  expect_equal(cal$lat, obs$lat)
  expect_equal(cal$lon, exp$lon)
  expect_equal(cal$lon, obs$lon)
  expect_error(
    CST_Calibration(exp = exp, obs = exp),
    "The length of the dimension 'member' in the component 'data' "
  )

  exp2 <- exp
  exp2$data[1, 2, 1, 1, 1, 1] <- NA
  expect_warning(
    CST_Calibration(exp = exp2, obs = obs),
    "Parameter 'exp' contains NA values."
  )

  obs2 <- obs
  obs2$data[1, 1, 2, 1, 1, 1] <- NA
  expect_warning(
    CST_Calibration(exp = exp, obs = obs2),
    "Parameter 'obs' contains NA values."
  )

  expect_warning(
    CST_Calibration(exp = exp2, obs = obs2),
    "Parameter 'obs' contains NA values", "Parameter 'exp' contains NA values."
  )
})  
