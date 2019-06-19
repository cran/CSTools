context("Generic tests")
test_that("Sanity checks", {
  expect_error(
    PlotForecastPDF(fcst, tercile.limits),
    "object 'tercile.limits' not found")
  expect_error(
    PlotForecastPDF(fcst, tercile.limits = c(0.25, 0.55)),
    "object 'fcst' not found")
  expect_error(
    PlotForecastPDF(fcst, tercile.limits = 10),
    "Parameter 'tercile.limits' should be an array with two limits for delimiting tercile categories")
  expect_error(
    PlotForecastPDF(fcst, tercile.limits = c(10, 20)),
    "object 'fcst' not found")
  fcsts2 <- array(rnorm(100),dim = c(members = 20, fcst = 5))
  expect_error(
    PlotForecastPDF(fcst = fcsts2, tercile.limits),
    "object 'tercile.limits' not found")
  expect_error( 
    PlotForecastPDF(fcst = fcsts2, tercile.limits = c(-0.5, 0.5), extreme.limits = NA),
    "Parameter 'extreme.limits' should be an array with two limits for delimiting extreme categories")
})

