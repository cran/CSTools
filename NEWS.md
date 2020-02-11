### CSTools 3.0.0
**Submission date to CRAN: 10-02-2020**

- New features:
    + CST_MergeDims and MergeDims
    + Version working with R 3.4.2
    + PlotForecastPDF handles independent terciles, extremes and observations for each panel
- Fixes
    + CST_Calibration handles missing values
    + BEI functions handle missing values
    
### CSTools 2.0.0 
**Submission date to CRAN: 25-11-2019**

- New features: 
    + CST_Analogs Analogs downscaling method, 
    + CST_MultiEOFS for multiple variables, 
    + Ensemble Clustering, 
    + Categorical Ensemble Combination,
    + new Calibration methods included in CST_Calibration, 
    + Best Estimated Index method, 
    + CST_QuantileMapping,
    + CST_SplitDim to split dimension, if it is a temporal dimension, it can be split by days, months and years or other inidices,
    + creation and transformation to class 's2dv_cube', 
    + CST_SaveExp function for saving experiments to be loadable with CST_Load, 
- Parallelization of RainFARM downscaling
- Adding unit tests using testthat for BEI and RainFarm functions
- New vignette Best Estimate Index
- Minor fix in CST_BiasCorrection when checking parameter 'obs'
- Addapting CST_Load to use 'as.s2dv_cube' function
- Minor fix in data lonlat_prec to be of class 's2dv_cube'
- Minor fix in RainFARM vignette
- Adding reference to S2S4E H2020 project into the DESCRIPTION file
- Adding NEWS.md file

### CSTools 1.0.1 
**Release date on CRAN: 19-06-2019**

- Correcting test of PlotForecastPDF for compatibility with ggplot2 release
- New function PlotCombinedMap 
- Adding reference to MEDSCOPE ERA4CS Project into the DESCRIPTION file
- Documentation minor fix in CST_RFWeights
- Minor fix in PlotMostLikelyQuantileMap for bar_titles
- MultiModelSkill vignette updated to use PlotCombinedMap



### CSTools 1.0.0 
**Release date on CRAN: 24-04-2019**

- Features included: Load, Anomaly, MultiMetric, MultivarRMSE, Calibration, BiasCorrection, RainFARM Downscaling, PlotForecastPDF, PlotMostLikelyQuantileMap
- Three sample data: lonlat_data, lonlat_prec, areave_data
- Unit tests using testthat: BiasCorrection, Calibration, MultiMetric, PlotForecast
- Vignettes: MultiMetric, Multivar and RainFARM
