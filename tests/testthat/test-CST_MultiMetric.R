context("Generic tests")
test_that("basic use case", {
  mod <- 1 : (2 * 3 * 4 * 5 * 6 * 8)
  dim(mod) <- c(dataset = 2, member = 3, sdate = 4, ftime = 5, lat = 6, lon = 8)
  obs <- 1 : (1 * 1 * 4 * 5 * 6 * 8)
  dim(obs) <- c(dataset = 1, member = 1, sdate = 4, ftime = 5, lat = 6, lon = 8)
  lon <- seq(0, 30, 5)
  lat <- seq(0, 30, 5)
  exp <- list(data = mod, lat = lat, lon = lon)
  obs <- list(data = obs, lat = lat, lon = lon)
  attr(exp, 'class') <- 's2dv_cube'
  attr(obs, 'class') <- 's2dv_cube'
  
  result <- list(data = array(rep(c(rep(1, 9), 
rep(0.89999999999999991118215802998747676610950, 3)), 48), 
                dim = c(dataset = 3, dataset = 1, statistics = 4, 
                       lat = 6, lon = 8)),
                lat = lat, lon = lon) 
  attr(result, 'class') <- 's2dv_cube'
  expect_equal(CST_MultiMetric(exp = exp, obs = obs), result)
  
  exp2 <- exp
  exp2$data[1, 1, 1, 2, 1, 1] = NA
  CST_MultiMetric(exp = exp2, obs = obs) 
  res <- CST_MultiMetric(exp = exp, obs = obs, metric = 'rms')
  expect_equal(length(res), 3)
  expect_equal(dim(res$data), 
        c(dataset = 3, dataset = 1, statistics = 3, lat = 6, lon = 8))
  res <- CST_MultiMetric(exp = exp, obs = obs, metric = 'rms', 
                         multimodel = FALSE)
  expect_equal(dim(res$data), 
        c(dataset = 2, dataset = 1, statistics = 3, lat = 6, lon = 8))
  res <- CST_MultiMetric(exp = exp, obs = obs, metric = 'rmsss')
  expect_equal(dim(res$data), 
        c(dataset = 3, dataset = 1, statistics = 2, lat = 6, lon = 8))
  res <- CST_MultiMetric(exp = exp, obs = obs, metric = 'rmsss', multimodel = FALSE)
  expect_equal(dim(res$data), 
        c(dataset = 2, dataset = 1, statistics = 2, lat = 6, lon =8))
  })


test_that("Sanity checks", {
  expect_error(
    CST_MultiMetric(exp = 1),
    paste0("Parameter 'exp' and 'obs' must be of the class 's2dv_cube', ",
         "as output by CSTools::CST_Load."))
  mod <- 1 : (2 * 3 * 4 * 5 * 6 * 8)
  dim(mod) <- c(dataset = 2, member = 3, sdate = 4, ftime = 5, lat = 6, lon = 8)
  obs <- 1 : (1 * 1 * 4 * 5 * 6 * 8)
  dim(obs) <- c(dataset = 1, member = 1, sdate = 4, ftime = 5, lat = 6, lon = 8)
  lon <- seq(0, 30, 5)
  lat <- seq(0, 30, 5)
  exp <- list(data = mod, lat = lat, lon = lon)
  obs <- list(data = obs, lat = lat, lon = lon)
  attr(exp, 'class') <- 's2dv_cube'
  attr(obs, 'class') <- 's2dv_cube'

  
  expect_error(
    CST_MultiMetric(exp = exp, obs = obs, metric = 1),
    paste0("Parameter 'metric' must be a character string indicating one ",
           "of the options: 'correlation', 'rms' or 'rmse'"))
  expect_error(
    CST_MultiMetric(exp = exp, obs = obs, metric = NA),
    "missing value where TRUE/FALSE needed")
  expect_error(
    CST_MultiMetric(exp = exp, obs = obs, metric = NULL),
    "argument is of length zero")
  expect_error(
    CST_MultiMetric(exp = exp, obs = obs, metric = "correlation", 
                    multimodel = NULL),
    "Parameter 'multimodel' must be a logical value.")
})  
