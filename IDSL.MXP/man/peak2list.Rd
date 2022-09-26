\name{peak2list}
\alias{peak2list}
\title{
Peak to List (The main function)
}
\description{
This function returns a list of two data objects required for the mass spectrometry data processing.
}
\usage{
peak2list(path = getwd(), MSfileName = "")
}
\arguments{
  \item{path}{
address of the mass spectrometry file
}
  \item{MSfileName}{
name of the mass spectrometry file with .mzML or .mzXML extensions
}
}
\details{
\url{https://ipa.idsl.me/mxp}
}
\value{
  \item{scanTable}{
a dataframe of different scan properties including `seqNum`, `msLevel`, `polarity`, `peaksCount`, `totIonCurrent`, `retentionTime`, `basePeakMZ`, `basePeakIntensity`, `collisionEnergy`, `lowMZ`, `highMZ`, `precursorScanNum`, `precursorMZ`, `precursorCharge`, `precursorIntensity`, `injectionTime`, `filterString`, `scanType`, `centroided`, `isolationWindowTargetMZ`, `isolationWindowLowerOffset`, `isolationWindowUpperOffset`, `scanWindowLowerLimit`, and `scanWindowUpperLimit`. `scanType` is only provided for the mzXML data format.
}
  \item{spectraList}{
a list of matrices of m/z and intensity values for each chromatogram scan
}
}
\note{
`retentionTime` column in the `scanTable' object is presented in minute.
}
\seealso{
\url{https://colab.research.google.com/drive/1gXwwuI1zzDHykKfodLSQQt5rwTuFEMpD}
}
\examples{
\donttest{
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd,"/idsl_ipa_test_files.zip")
download.file(paste0("https://github.com/idslme/IDSL.IPA/blob/main/",
"IPA_educational_files/idsl_ipa_test_files.zip?raw=true"),
destfile = temp_wd_zip, mode = "wb")
unzip(temp_wd_zip, exdir = temp_wd)
p2l <- peak2list(path = temp_wd, MSfileName = "003.mzML")
}
}