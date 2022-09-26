\name{getSpectra}
\alias{getSpectra}
\title{
getSpectra
}
\description{
This function creates a spectraList for the chromatogram scans of the mass spectrometry data.
}
\usage{
getSpectra(xmlData, msFormat)
}
\arguments{
  \item{xmlData}{
a structured data of the mass spectrometry data created by the `read_xml` function.
}
  \item{msFormat}{
format extension of the mass spectrometry file c("mzML", "mzXML")
}
}
\value{
a list of matrices of m/z and intensity values for each chromatogram scan
}
\examples{
\donttest{
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd,"/idsl_ipa_test_files.zip")
download.file(paste0("https://github.com/idslme/IDSL.IPA/blob/main/",
"IPA_educational_files/idsl_ipa_test_files.zip?raw=true"),
destfile = temp_wd_zip, mode = "wb")
unzip(temp_wd_zip, exdir = temp_wd)
xmlData <- xml2::read_xml(paste0(path = temp_wd, "/", MSfile = "003.mzML"))
spectraList <- getSpectra(xmlData, msFormat = "mzML")
}
}