context("Generic tests")
test_that("some checks", {
  fcsts2 <- array(1:100, dim = c(members = 20, fcst = 5))
  expect_equal(
    length(PlotForecastPDF(fcsts2, c(-0.66, 0.66), extreme.limits = c(-1.2, 1.2), 
                    fcst.names = paste0("random fcst ", 1 : 5), obs = 0.7)),
    10)
  expect_equal(
    names(PlotForecastPDF(fcsts2, c(-0.66, 0.66), extreme.limits = c(-1.2, 1.2), 
                           fcst.names = paste0("random fcst ", 1 : 5), obs = 0.7)),
    c("data", "layers", "scales", "mapping", "theme", "coordinates", "facet", 
      "plot_env", "labels", "guides"))
  expect_equal(
    length(PlotForecastPDF(fcsts2, c(-0.66, 0.66), extreme.limits = c(-1.2, 1.2), 
                    fcst.names = paste0("random fcst ", 1 : 5), obs = 0.7)[[5]]),
    61)
})

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

