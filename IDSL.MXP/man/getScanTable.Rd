\name{getScanTable}
\alias{getScanTable}
\title{
getScanTable
}
\description{
This function creates a scanTable for the chromatogram scans of the mass spectrometry data.
}
\usage{
getScanTable(xmlData, msFormat)
}
\arguments{
  \item{xmlData}{
A structured data of the mass spectrometry data created by the `read_xml` function.
}
  \item{msFormat}{
format extension of the mass spectrometry file c("mzML", "mzXML")
}
}
\value{
a dataframe of different scan properties including `seqNum`, `msLevel`, `polarity`, `peaksCount`, `totIonCurrent`, `retentionTime`, `basePeakMZ`, `basePeakIntensity`, `collisionEnergy`, `lowMZ`, `highMZ`, `precursorScanNum`, `precursorMZ`, `precursorCharge`, `precursorIntensity`, `injectionTime`, `filterString`, `scanType`, `centroided`, `isolationWindowTargetMZ`, `isolationWindowLowerOffset`, `isolationWindowUpperOffset`, `scanWindowLowerLimit`, and `scanWindowUpperLimit`. `scanType` is only provided for the mzXML data format.
}
\note{
`retentionTime` column is presented in minute.
}
\examples{
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd,"/testfiles.zip")
download.file(
  "https://github.com/idslme/IDSL.IPA/raw/main/idsl_ipa_test_files.zip",
  destfile = temp_wd_zip)
unzip(temp_wd_zip, exdir = temp_wd)
xmlData <- xml2::read_xml(paste0(path = temp_wd, "/", MSfile = "003.mzML"))
scanTable <- getScanTable(xmlData, msFormat = "mzML")
}