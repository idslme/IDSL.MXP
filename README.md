# IDSL.MXP<img src='MXP_educational_files/Figures/IDSL.MXP-logo.png' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.MXP)](https://cran.r-project.org/package=IDSL.MXP)
![](http://cranlogs.r-pkg.org/badges/IDSL.MXP?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.MXP?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.MXP)](https://cran.r-project.org/package=IDSL.MXP)
<!-- badges: end -->

A tiny parser to extract mass spectra data and metadata table of MS acquisition properties from mzML, mzXML and netCDF mass spectrometry files.

	install.packages("IDSL.MXP")

Visit https://ipa.idsl.me/mxp for the detailed documentation and tutorial.

To use this package, follow below commands:

	path <- “location address of the mass spectrometry file”
	msfile <- "name of the mass spectrometry file with its extension"
	msobject <- IDSL.MXP::peak2list(path, msfile)

msobject is a list with two objects - 1) scanTable, a data.frame of different scan properties, and 2) spectraList, a list of m/z and intensity values for each scan