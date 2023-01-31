peak2list <- function(path, MSfileName = "") {
  ##
  MSfileLocation <- paste0(path, "/", MSfileName)
  MSfileLocation <- gsub("\\", "/", MSfileLocation, fixed = TRUE)
  strMSfileLocation <- strsplit(MSfileLocation, "/")[[1]]
  MSfileName <- strMSfileLocation[length(strMSfileLocation)]
  MSfileLocation <- paste0(strMSfileLocation, collapse = "/")
  ##
  msFormat <- strsplit(MSfileName, "[.]")[[1]]
  msFormat <- tolower(msFormat[length(msFormat)])
  ##
  if ((msFormat == "mzml") | (msFormat == "mzxml")) {
    ##
    xmlData <- read_xml(MSfileLocation)
    ##
    scanTable <- getScanTable(xmlData, msFormat)
    ##
    spectraList <- getSpectra(xmlData, msFormat)
    ##
    p2l <- list(scanTable, spectraList)
    ##
    names(p2l) <- c("scanTable", "spectraList")
    ##
  } else if ((msFormat == "cdf")) {
    ##
    p2l <- getNetCDF(MSfileLocation)
    ##
  } else {
    stop(paste0(MSfileName, " is not consistent with the IDSL.MXP package!"))
  }
  return(p2l)
}