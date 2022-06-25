\name{getNetCDF}
\alias{getNetCDF}
\title{
getNetCDF
}
\description{
This function returns a list of two data objects needed for the mass spectrometry data processing.
}
\usage{
getNetCDF(MSfile)
}
\arguments{
  \item{MSfile}{
name of the mass spectrometry file with .cdf extension
}
}
\value{
  \item{scanTable}{
a dataframe of different scan properties including `seqNum`, `msLevel`, `polarity`, `peaksCount`, `totIonCurrent`, `retentionTime`, `basePeakMZ`, `basePeakIntensity`, `collisionEnergy`, `lowMZ`, `highMZ`, `precursorScanNum`, `precursorMZ`, `precursorCharge`, `precursorIntensity`, `injectionTime`, `filterString`, `scanType`, `centroided`, `isolationWindowTargetMZ`, `isolationWindowLowerOffset`, `isolationWindowUpperOffset`, `scanWindowLowerLimit`, and `scanWindowUpperLimit`.
}
  \item{spectraList}{
a list of matrices of m/z and intensity values for each chromatogram scan
}
}
\note{
`retentionTime` column in the `scanTable' object is presented in minute.
}