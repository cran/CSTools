### CSTools 4.0.0 
**Submission date to CRAN: XX-12-2020**

- New features:
    + ADAMONT downscaling method: requires CST_AdamontAnalogs and CST_AdamontQQCor functions
    + Analogs method using Predictors: requires training_analogs and  CST_AnalogsPredictors
    + PlotPDFsOLE includes parameters to modify legend style
    + CST_RFSlope handless missing values in the temporal dimension and new 'ncores' parameter allows parallel computation
    + CST_RFWeights accepts s2dv_cube objects as input and new 'ncores' paramenter allows parallel computation
    + RFWeights is exposed to users
    + CST_RainFARM accepts multi-dimensional slopes and weights and handless missing values in sample dimensions.
    + QuantileMapping is exposed to users
    + CST_MultiMetric includes 'rpss' metric and it is addapted to s2dv.
    + PlotMostLikelyQuantileMap vignette
    + PlotTriangles4Categories includes two parameters to adjust axis and margins
    + CategoricalEnsCombination is exposed to users
    + CST_SplitDims includes parameter 'insert_ftime'
    + Analogs vignette
    + Data Storage and retrieval vignette

- Fixes:
    + PlotForecastPDF correctly displays terciles labels 
    + CST_SaveExp correctly save time units
    + CST_SplitDims returns ordered output following ascending order provided in indices when it is numeric
    + qmap library moved from Imports to Depends
    + CST_QuantileMapping correctly handles exp_cor
    + Figures resize option from vignettes has been removed
    + Fix Analogs to work with three diferent criteria
    + Vignette PlotForecastPDF updated plots
    + Decrease package size compresing vignettes figures and removing areave_data sample

### CSTools 3.1.0
**Submission date to CRAN: 02-07-2020**

- New features:
    + EnsClustering vignette
    + EnsClustering has a new parameter 'time_dim'
    + CST_BiasCorrection has na.rm paramter
    + CST_Anomaly allows to smooth the climatology with filter.span parameter
    + PlotTriangles4Categories new plotting function to convert any 3-d numerical array to a grid of coloured triangles.
    + CST_WeatherRegimes/WeatherRegimes and CST_RegimeAssign/RegimeAssign
    + PlotPDFsOLE plots two probability density gaussian functions and the optimal linear estimation
    + CST_RFTemp/RF_Temp functions available for downscaling temperature
    + Weather Regimes vignette

- Fixes
    + CST_Anomaly handles exp, obs or both
    + PlotForecastPDF vignette displays figures correctly
    + Calibration function is exposed to users
    + MultiMetric vignette fixed typo text description
    + RainFARM checks 'slope' is not a vector
    + DESCRIPTION specifies the minimum multiApply version required
    + EnsClustering has a fixed 'closest_member' output
    + PlotCombinedMap handles masks correctly
    + CST_SaveExp uses multiApply and save time dimension correctly

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
