# IDSL.MXP <img src='MXP_educational_files/Figures/IDSL.MXP-logo.png' width="250px" align="right" />

<!-- badges: start -->
[![Developed-by](https://img.shields.io/badge/Developed_by-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.MXP)](https://cran.r-project.org/package=IDSL.MXP)
![](http://cranlogs.r-pkg.org/badges/IDSL.MXP?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.MXP?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.MXP)](https://cran.r-project.org/package=IDSL.MXP)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1gXwwuI1zzDHykKfodLSQQt5rwTuFEMpD)
<!-- badges: end -->

[**Mass Spectrometry Parser (MXP)**](https://ipa.idsl.me/mxp) by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me/) is a tiny parser to extract mass spectra data and metadata table of instrumentational acquisition properties from mzML, mzXML and netCDF mass spectrometry files.

## Installation

	install.packages("IDSL.MXP")
	
**Note:** In case you want to process **netCDF/CDF** mass spectrometry data by IDSL.MXP, you should also install the [**RnetCDF**](https://CRAN.R-project.org/package=RNetCDF) package separately using the below command.

	install.packages("RNetCDF")

## Workflow
To use this package, follow below commands:

	path <- “location address of the mass spectrometry file”
	MSfileName <- "name of the mass spectrometry file with its extension"
	mxpObject <- IDSL.MXP::peak2list(path, MSfileName)

**mxpObject** is a list with two objects - **1) scanTable**, a metadata of different scan properties, and **2) spectraList**, a list of m/z and intensity values for each scan.

Visit [**wiki**](https://github.com/idslme/IDSL.MXP/wiki/Example-for-IDSL.MXP) and [Google colab](https://colab.research.google.com/drive/1gXwwuI1zzDHykKfodLSQQt5rwTuFEMpD) to illustrate performance of IDSL.MXP.

## Citation

[1] Fakouri Baygi, S., Kumar, Y. Barupal, D.K. [IDSL. IPA characterizes the organic chemical space in untargeted LC/HRMS datasets](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00120). *Journal of proteome research*, **2022**, *21(6)*, 1485-1494.