CSTools
=======

#### Assessing Skill of Climate Forecasts on Seasonal-to-Decadal Timescales

The Climate Services Tools, CSTools, is an easy-to-use R package designed and built to assess and improve the quality of climate forecasts for seasonal to multi–annual scales. The package contains process-based state-of-the-art methods for forecast calibration, bias correction, statistical and stochastic downscaling, optimal forecast combination and multivariate verification, as well as basic and advanced tools to obtain tailored products.

This package was developed in the context of the ERA4CS project MEDSCOPE and the H2020 S2S4E project and includes contributions from ArticXchange project founded by EU-PolarNet 2. This GitLab project allows you to monitor its progress and to interact with other developers via the Issues section.

A scientific publication including use cases was published in the Geoscientific Model Development Journal, and it can be cited as follows:

> Pérez-Zanón, N., Caron, L.-P., Terzago, S., Van Schaeybroeck, B., Lledó, L., Manubens, N., Roulin, E., Alvarez-Castro, M. C., Batté, L., Bretonnière, P.-A., Corti, S., Delgado-Torres, C., Domínguez, M., Fabiano, F., Giuntoli, I., von Hardenberg, J., Sánchez-García, E., Torralba, V., and Verfaillie, D.: Climate Services Toolbox (CSTools) v4.0: from climate forecasts to climate forecast information, Geosci. Model Dev., 15, 6115–6142, https://doi.org/10.5194/gmd-15-6115-2022, 2022.

On-line resources
-----------------

A part from this GitLab project, that allows you to monitor CSTools progress, to interact with other developers via the Issues section and to contribute, you can find:

- The CRAN repository [https://CRAN.R-project.org/package=CSTools](https://CRAN.R-project.org/package=CSTools) which includes the user manual and vignettes.
- Video tutorials [https://www.medscope-project.eu/products/tool-box/cstools-video-tutorials/](https://www.medscope-project.eu/products/tool-box/cstools-video-tutorials/).
- Other resources are under-development such [training material](https://earth.bsc.es/gitlab/external/cstools/-/tree/MEDCOF2022/inst/doc/MEDCOF2022) and a [full reproducible use case for forecast calibration](https://earth.bsc.es/gitlab/external/cstools/-/tree/develop-CalibrationVignette/FOCUS_7_2).

Installation
------------

CSTools has a system dependency, the CDO libraries, for interpolation of grid data
and retrieval of metadata. Make sure you have these libraries installed in the
system or download and install from
[https://code.zmaw.de/projects/cdo](https://code.zmaw.de/projects/cdo).

You can then install the public released version of CSTools from CRAN:

```r
install.packages("CSTools")
```

Or the development version from the GitLab repository:

```r
# install.packages("devtools")
devtools::install_git("https://earth.bsc.es/gitlab/external/cstools.git")
```

Overview
--------

The CSTools package functions can be distributed in the following methods:

- **Data retrieval and formatting:** CST_Load, CST_Anomaly, CST_MergeDims, CST_SplitDims, CST_Subset, as.s2dv_cube, s2dv_cube, CST_SaveExp.
- **Classification:** CST_MultiEOF, CST_WeatherRegimes, CST_RegimsAssign, CST_CategoricalEnsCombination, CST_EnsClustering.
- **Downscaling:** CST_Analogs, CST_RainFARM, CST_RFTemp, CST_AdamontAnalog, CST_AnalogsPredictors.
- **Correction:** CST_BEI_Weighting, CST_BiasCorrection, CST_Calibration, CST_QuantileMapping, CST_DynBiasCorrection.
- **Assessment:** CST_MultiMetric, CST_MultivarRMSE
- **Visualization:** PlotCombinedMap, PlotForecastPDF, PlotMostLikelyQuantileMap, PlotPDFsOLE, PlotTriangles4Categories, PlotWeeklyClim.

This package is designed to be compatible with other R packages such as [s2dv](https://CRAN.R-project.org/package=s2dv), [startR](https://CRAN.R-project.org/package=startR), [CSIndicators](https://CRAN.R-project.org/package=CSIndicators), [CSDownscale](https://earth.bsc.es/gitlab/es/csdownscale). Functions with the prefix **CST_** deal with a common object called `s2dv_cube` as inputs. Also, this object can be created from Load (s2dv) and from Start (startR) directly. Multiple functions from different packages can operate on this common data structure to easily define a complete post-processing workflow.

The class `s2dv_cube` is mainly a list of named elements to keep data and metadata in a single object. Basic structure of the object:

```r
$ data: [data array]
$ dims: [dimensions vector]
$ coords: [List of coordinates vectors]
  $ sdate
  $ time
  $ lon
  [...]
$ attrs: [List of the attributes]
  $ Variable:
    $ varName
    $ metadata 
  $ Datasets
  $ Dates
  $ source_files
  $ when
  $ load_parameters
```

More information about the `s2dv_cube` object class can be found here: [description of the s2dv_cube object structure document](https://docs.google.com/document/d/1ko37JFl_h6mOjDKM5QSQGikfLBKZq1naL11RkJIwtMM/edit?usp=sharing).

The current `s2dv_cube` object (CSTools 5.0.0) differs from the original object used in the previous versions of the packages. If you have **questions** on this change you can follow some of the points below:

- [New s2dv_cube object discussion](https://earth.bsc.es/gitlab/external/cstools/-/issues/94)
- [How to deal with the compatibility break](https://earth.bsc.es/gitlab/external/cstools/-/issues/112)
- [Testing issue and specifications](https://earth.bsc.es/gitlab/external/cstools/-/issues/110)

Contribute
----------

Before adding a development, we suggest to contact the package mantainer. Details on the procedure and development guidelines can be found in [this issue](https://earth.bsc.es/gitlab/external/cstools/-/issues/3).

If you plan on contributing, you should rather clone the project on your workstation and modify it using the basic Git commands (clone, branch, add, commit, push, merge, ...).

The code of each function should live in a separate file with the .R extension under the R folder, and the documentation of each function should live in a separate file with the .Rd extension under the man folder.

For an introductory video on Git, you can have a look at [https://vimeo.com/41027679](https://vimeo.com/41027679).

You can also find all the necessary documentation on git here: [https://git-scm.com/book/en/v2](https://git-scm.com/book/en/v2). A lot of it may be a bit complicated for beginners (and not necessary for us), but the "Getting started" and "Git basics" sections are a good resources.
