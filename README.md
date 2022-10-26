Welcome to the CSTools GitLab website
======================================

The Climate Services Tools, CSTools, is an easy-to-use R package designed and built to assess and improve the quality of climate forecasts for seasonal to multi–annual scales. The package contains process-based state-of-the-art methods for forecast calibration, bias correction, statistical and stochastic downscaling, optimal forecast combination and multivariate verification, as well as basic and advanced tools to obtain tailored products. 

This package was developed in the context of the ERA4CS project MEDSCOPE and the H2020 S2S4E project. This GitLab project allows you to monitor its progress and to interact with other developers via the Issues section.

A scientific publication including use cases was published in the Geoscientific Model Development Journal, and it can be cited as follows:  
> Pérez-Zanón, N., Caron, L.-P., Terzago, S., Van Schaeybroeck, B., Lledó, L., Manubens, N., Roulin, E., Alvarez-Castro, M. C., Batté, L., Bretonnière, P.-A., Corti, S., Delgado-Torres, C., Domínguez, M., Fabiano, F., Giuntoli, I., von Hardenberg, J., Sánchez-García, E., Torralba, V., and Verfaillie, D.: Climate Services Toolbox (CSTools) v4.0: from climate forecasts to climate forecast information, Geosci. Model Dev., 15, 6115–6142, https://doi.org/10.5194/gmd-15-6115-2022, 2022.


On-line resources
-----------------

A part from this GitLab project, that allows you to monitor CSTools progress, to interact with other developers via the Issues section and to contribute, you can find:

- The CRAN repository <https://CRAN.R-project.org/package=CSTools> which includes the user manual and vignettes

- Video tutorials <https://www.medscope-project.eu/products/tool-box/cstools-video-tutorials/>

- Other resources are under-development such [training material](https://earth.bsc.es/gitlab/external/cstools/-/tree/MEDCOF2022/inst/doc/MEDCOF2022) and a [full reproducible use case for forecast calibration](https://earth.bsc.es/gitlab/external/cstools/-/tree/develop-CalibrationVignette/FOCUS_7_2)

Installation
------------

CSTools has a system dependency, the CDO libraries, for interpolation of grid data
and retrieval of metadata. Make sure you have these libraries installed in the
system or download and install from
<https://code.zmaw.de/projects/cdo>.

You can then install the public released version of CSTools from CRAN:
```r
install.packages("CSTools")
```
Or the development version from the GitLab repository:
```r
# install.packages("devtools")
devtools::install_git("https://earth.bsc.es/gitlab/external/cstools.git")
```

How to contribute
-----------------

Before adding a development, we suggest to contact the package mantainer. Details on the procedure and development guidelines can be found in [this issue](https://earth.bsc.es/gitlab/external/cstools/-/issues/3)

If you plan on contributing, you should rather clone the project on your workstation and modify it using the basic Git commands (clone, branch, add, commit, push, merge, ...).

The code of each function should live in a separate file with the .R extension under the R folder, and the documentation of each function should live in a separate file with the .Rd extension under the man folder.

For an introductory video on Git, you can have a look at https://vimeo.com/41027679.

You can also find all the necessary documentation on git here: https://git-scm.com/book/en/v2
A lot of it may be a bit complicated for beginners (and not necessary for us), but the "Getting started" and "Git basics" sections are a good resources.

