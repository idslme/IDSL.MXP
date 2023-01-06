getNetCDF <- function(MSfile) {
  ##
  RNetCDFpackageCheck <- tryCatch(requireNamespace('RNetCDF', quietly = TRUE), error = function(e) {FALSE})
  if (!RNetCDFpackageCheck) {
    message("IDSL.MXP requires the 'RNetCDF' package to deconvolute `.netCDF` mass spectrometry data!")
    stop(" <<< install.packages('RNetCDF') >>> ")
  }
  netCDFmassSpectrometryfile <- RNetCDF::open.nc(MSfile)
  netCDFList <- RNetCDF::read.nc(netCDFmassSpectrometryfile)
  nameNetCDFList <- names(netCDFList)
  nameNetCDFList <- setdiff(nameNetCDFList, c("mass_values", "intensity_values"))
  ##############################################################################
  scanTable <- do.call(cbind, lapply(nameNetCDFList, function(n) {
    netCDFList[[n]]
  }))
  ##
  scanTable <- data.frame(scanTable)
  totalNumberScans <- dim(scanTable)[1]
  ##
  x_peaksCount <- grep("point_count", nameNetCDFList, ignore.case = TRUE)
  if (length(x_peaksCount) > 0) {
    nameNetCDFList[x_peaksCount] <- "peaksCount"
    scanTable[, x_peaksCount] <- as.numeric(scanTable[, x_peaksCount])
  } else {
    stop("`point_count` column is not available!")
  }
  ##############################################################################
  netCDFdataFrame <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
  colnames(netCDFdataFrame) <- "toBeRemoved"
  ##############################################################################
  x_filterString <- grep("error_log", nameNetCDFList, ignore.case = TRUE)
  if (length(x_filterString) > 0) {
    nameNetCDFList[x_filterString] <- "filterString"
    scanTable[, x_filterString] <- scanTable[, x_filterString]
  } else {
    filterStringColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(filterStringColumn) <- "filterString"
    netCDFdataFrame <- cbind(netCDFdataFrame, filterStringColumn)
  }
  ##############################################################################
  x_seqNum <- grep("seqNum", nameNetCDFList, ignore.case = TRUE)
  if (length(x_seqNum) > 0) {
    nameNetCDFList[x_seqNum] <- "seqNum"
    scanTable[, x_seqNum] <- as.numeric(scanTable[, x_seqNum])
    ##
    if (scanTable[1, x_seqNum] < 0) {
      scanTable[, x_seqNum] <- matrix(seq(1, totalNumberScans, 1), ncol = 1)
    }
    ##
  } else {
    seqNumColumn <- data.frame(matrix(seq(1, totalNumberScans, 1), ncol = 1))
    colnames(seqNumColumn) <- "seqNum"
    netCDFdataFrame <- cbind(netCDFdataFrame, seqNumColumn)
  }
  ##############################################################################
  x_acquisitionNum <- grep("actual_scan_number", nameNetCDFList, ignore.case = TRUE)
  if (length(x_acquisitionNum) > 0) {
    nameNetCDFList[x_acquisitionNum] <- "acquisitionNum"
    scanTable[, x_acquisitionNum] <- as.numeric(scanTable[, x_acquisitionNum])
    ##
    if (scanTable[1, x_acquisitionNum] < 0) {
      scanTable[, x_acquisitionNum] <- matrix(seq(1, totalNumberScans, 1), ncol = 1)
    }
    ##
  } else {
    acquisitionNumColumn <- data.frame(matrix(seq(1, totalNumberScans, 1), ncol = 1))
    colnames(acquisitionNumColumn) <- "acquisitionNum"
    netCDFdataFrame <- cbind(netCDFdataFrame, acquisitionNumColumn)
  }
  ##############################################################################
  x_totIonCurrent <- grep("total_intensity", nameNetCDFList, ignore.case = TRUE)
  if (length(x_totIonCurrent) > 0) {
    nameNetCDFList[x_totIonCurrent] <- "totIonCurrent"
    scanTable[, x_totIonCurrent] <- as.numeric(scanTable[, x_totIonCurrent])
  } else {
    totIonCurrentColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(totIonCurrentColumn) <- "totIonCurrent"
    netCDFdataFrame <- cbind(netCDFdataFrame, totIonCurrentColumn)
  }
  ##############################################################################
  x_retentionTime <- grep("scan_acquisition_time", nameNetCDFList, ignore.case = TRUE)
  if (length(x_retentionTime) > 0) {
    nameNetCDFList[x_retentionTime] <- "retentionTime"
    scanTable[, x_retentionTime] <- as.numeric(scanTable[, x_retentionTime])/60
  } else {
    retentionTimeColumn <- data.frame(matrix(rep(0, totalNumberScans), ncol = 1))
    colnames(retentionTimeColumn) <- "retentionTime"
    netCDFdataFrame <- cbind(netCDFdataFrame, retentionTimeColumn)
  }
  ##############################################################################
  x_lowMZ <- grep("mass_range_min", nameNetCDFList, ignore.case = TRUE)
  if (length(x_lowMZ) > 0) {
    nameNetCDFList[x_lowMZ] <- "lowMZ"
    scanTable[, x_lowMZ] <- as.numeric(scanTable[, x_lowMZ])
  } else {
    lowMZColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(lowMZColumn) <- "lowMZ"
    netCDFdataFrame <- cbind(netCDFdataFrame, lowMZColumn)
  }
  ##############################################################################
  x_highMZ <- grep("mass_range_max", nameNetCDFList, ignore.case = TRUE)
  if (length(x_highMZ) > 0) {
    nameNetCDFList[x_highMZ] <- "highMZ"
    scanTable[, x_highMZ] <- as.numeric(scanTable[, x_highMZ])
  } else {
    highMZColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(highMZColumn) <- "highMZ"
    netCDFdataFrame <- cbind(netCDFdataFrame, highMZColumn)
  }
  ##############################################################################
  x_scanWindowLowerLimit <- grep("time_range_min", nameNetCDFList, ignore.case = TRUE)
  if (length(x_scanWindowLowerLimit) > 0) {
    nameNetCDFList[x_scanWindowLowerLimit] <- "scanWindowLowerLimit"
    scanTable[, x_scanWindowLowerLimit] <- as.numeric(scanTable[, x_scanWindowLowerLimit])
  } else {
    scanWindowLowerLimitColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(scanWindowLowerLimitColumn) <- "scanWindowLowerLimit"
    netCDFdataFrame <- cbind(netCDFdataFrame, scanWindowLowerLimitColumn)
  }
  ##############################################################################
  x_scanWindowUpperLimit <- grep("time_range_max", nameNetCDFList, ignore.case = TRUE)
  if (length(x_scanWindowUpperLimit) > 0) {
    nameNetCDFList[x_scanWindowUpperLimit] <- "scanWindowUpperLimit"
    scanTable[, x_scanWindowUpperLimit] <- as.numeric(scanTable[, x_scanWindowUpperLimit])
  } else {
    scanWindowUpperLimitColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(scanWindowUpperLimitColumn) <- "scanWindowUpperLimit"
    netCDFdataFrame <- cbind(netCDFdataFrame, scanWindowUpperLimitColumn)
  }
  ##############################################################################
  x_msLevel <- grep("ms_level", nameNetCDFList, ignore.case = TRUE)
  if (length(x_msLevel) > 0) {
    nameNetCDFList[x_msLevel] <- "msLevel"
    scanTable[, x_msLevel] <- as.numeric(scanTable[, x_msLevel])
  } else {
    msLevelColumn <- data.frame(matrix(rep(1, totalNumberScans), ncol = 1))
    colnames(msLevelColumn) <- "msLevel"
    netCDFdataFrame <- cbind(netCDFdataFrame, msLevelColumn)
  }
  ##############################################################################
  x_polarity <- grep("polarity", nameNetCDFList, ignore.case = TRUE)
  if (length(x_polarity) > 0) {
    nameNetCDFList[x_polarity] <- "polarity"
    scanTable[, x_polarity] <- as.numeric(scanTable[, x_polarity])
  } else {
    polarityColumn <- data.frame(matrix(rep(1, totalNumberScans), ncol = 1))
    colnames(polarityColumn) <- "polarity"
    netCDFdataFrame <- cbind(netCDFdataFrame, polarityColumn)
  }
  ##############################################################################
  x_basePeakMZ <- grep("basePeakMZ", nameNetCDFList, ignore.case = TRUE)
  if (length(x_basePeakMZ) > 0) {
    nameNetCDFList[x_basePeakMZ] <- "basePeakMZ"
    scanTable[, x_basePeakMZ] <- as.numeric(scanTable[, x_basePeakMZ])
  } else {
    basePeakMZColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(basePeakMZColumn) <- "basePeakMZ"
    netCDFdataFrame <- cbind(netCDFdataFrame, basePeakMZColumn)
  }
  ##############################################################################
  x_basePeakIntensity <- grep("basePeakIntensity", nameNetCDFList, ignore.case = TRUE)
  if (length(x_basePeakIntensity) > 0) {
    nameNetCDFList[x_basePeakIntensity] <- "basePeakIntensity"
    scanTable[, x_basePeakIntensity] <- as.numeric(scanTable[, x_basePeakIntensity])
  } else {
    basePeakIntensityColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(basePeakIntensityColumn) <- "basePeakIntensity"
    netCDFdataFrame <- cbind(netCDFdataFrame, basePeakIntensityColumn)
  }
  ##############################################################################
  x_precursorScanNum <- grep("precursorScanNum", nameNetCDFList, ignore.case = TRUE)
  if (length(x_precursorScanNum) > 0) {
    nameNetCDFList[x_precursorScanNum] <- "precursorScanNum"
    scanTable[, x_precursorScanNum] <- as.numeric(scanTable[, x_precursorScanNum])
  } else {
    precursorScanNumColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(precursorScanNumColumn) <- "precursorScanNum"
    netCDFdataFrame <- cbind(netCDFdataFrame, precursorScanNumColumn)
  }
  ##############################################################################
  x_precursorMZ <- grep("precursorMZ", nameNetCDFList, ignore.case = TRUE)
  if (length(x_precursorMZ) > 0) {
    nameNetCDFList[x_precursorMZ] <- "precursorMZ"
    scanTable[, x_precursorMZ] <- as.numeric(scanTable[, x_precursorMZ])
  } else {
    precursorMZColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(precursorMZColumn) <- "precursorMZ"
    netCDFdataFrame <- cbind(netCDFdataFrame, precursorMZColumn)
  }
  ##############################################################################
  x_precursorIntensity <- grep("precursorIntensity", nameNetCDFList, ignore.case = TRUE)
  if (length(x_precursorIntensity) > 0) {
    nameNetCDFList[x_precursorIntensity] <- "precursorIntensity"
    scanTable[, x_precursorIntensity] <- as.numeric(scanTable[, x_precursorIntensity])
  } else {
    precursorIntensityColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(precursorIntensityColumn) <- "precursorIntensity"
    netCDFdataFrame <- cbind(netCDFdataFrame, precursorIntensityColumn)
  }
  ##############################################################################
  x_precursorCharge <- grep("precursorCharge", nameNetCDFList, ignore.case = TRUE)
  if (length(x_precursorCharge) > 0) {
    nameNetCDFList[x_precursorCharge] <- "precursorCharge"
    scanTable[, x_precursorCharge] <- as.numeric(scanTable[, x_precursorCharge])
  } else {
    precursorChargeColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(precursorChargeColumn) <- "precursorCharge"
    netCDFdataFrame <- cbind(netCDFdataFrame, precursorChargeColumn)
  }
  ##############################################################################
  x_injectionTime <- grep("injectionTime", nameNetCDFList, ignore.case = TRUE)
  if (length(x_injectionTime) > 0) {
    nameNetCDFList[x_injectionTime] <- "injectionTime"
    scanTable[, x_injectionTime] <- as.numeric(scanTable[, x_injectionTime])
  } else {
    injectionTimeColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(injectionTimeColumn) <- "injectionTime"
    netCDFdataFrame <- cbind(netCDFdataFrame, injectionTimeColumn)
  }
  ##############################################################################
  x_isolationWindowTargetMZ <- grep("isolationWindowTargetMZ", nameNetCDFList, ignore.case = TRUE)
  if (length(x_isolationWindowTargetMZ) > 0) {
    nameNetCDFList[x_isolationWindowTargetMZ] <- "isolationWindowTargetMZ"
    scanTable[, x_isolationWindowTargetMZ] <- as.numeric(scanTable[, x_isolationWindowTargetMZ])
  } else {
    isolationWindowTargetMZColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(isolationWindowTargetMZColumn) <- "isolationWindowTargetMZ"
    netCDFdataFrame <- cbind(netCDFdataFrame, isolationWindowTargetMZColumn)
  }
  ##############################################################################
  x_isolationWindowLowerOffset <- grep("isolationWindowLowerOffset", nameNetCDFList, ignore.case = TRUE)
  if (length(x_isolationWindowLowerOffset) > 0) {
    nameNetCDFList[x_isolationWindowLowerOffset] <- "isolationWindowLowerOffset"
    scanTable[, x_isolationWindowLowerOffset] <- as.numeric(scanTable[, x_isolationWindowLowerOffset])
  } else {
    isolationWindowLowerOffsetColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(isolationWindowLowerOffsetColumn) <- "isolationWindowLowerOffset"
    netCDFdataFrame <- cbind(netCDFdataFrame, isolationWindowLowerOffsetColumn)
  }
  ##############################################################################
  x_isolationWindowUpperOffset <- grep("isolationWindowUpperOffset", nameNetCDFList, ignore.case = TRUE)
  if (length(x_isolationWindowUpperOffset) > 0) {
    nameNetCDFList[x_isolationWindowUpperOffset] <- "isolationWindowUpperOffset"
    scanTable[, x_isolationWindowUpperOffset] <- as.numeric(scanTable[, x_isolationWindowUpperOffset])
  } else {
    isolationWindowUpperOffsetColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(isolationWindowUpperOffsetColumn) <- "isolationWindowUpperOffset"
    netCDFdataFrame <- cbind(netCDFdataFrame, isolationWindowUpperOffsetColumn)
  }
  ##############################################################################
  x_collisionEnergy <- grep("collisionEnergy", nameNetCDFList, ignore.case = TRUE)
  if (length(x_collisionEnergy) > 0) {
    nameNetCDFList[x_collisionEnergy] <- "collisionEnergy"
    scanTable[, x_collisionEnergy] <- as.numeric(scanTable[, x_collisionEnergy])
  } else {
    collisionEnergyColumn <- data.frame(matrix(rep(NA, totalNumberScans), ncol = 1))
    colnames(collisionEnergyColumn) <- "collisionEnergy"
    netCDFdataFrame <- cbind(netCDFdataFrame, collisionEnergyColumn)
  }
  ##############################################################################
  x_centroided <- grep("centroided", nameNetCDFList, ignore.case = TRUE)
  if (length(x_centroided) > 0) {
    nameNetCDFList[x_centroided] <- "centroided"
    scanTable[, x_centroided] <- as.numeric(scanTable[, x_centroided])
  } else {
    centroidedColumn <- data.frame(matrix(rep(TRUE, totalNumberScans), ncol = 1))
    colnames(centroidedColumn) <- "centroided"
    netCDFdataFrame <- cbind(netCDFdataFrame, centroidedColumn)
  }
  ##############################################################################
  colnames(scanTable) <- nameNetCDFList
  ##
  numberProperties <- dim(netCDFdataFrame)[2]
  if (numberProperties > 1) {
    netCDFdataFrame <- netCDFdataFrame[, 2:numberProperties]
    scanTable <- cbind(scanTable, netCDFdataFrame)
  }
  ##############################################################################
  ##############################################################################
  mzPeakCountsIndex <- rep(0, (totalNumberScans + 1))
  for (i in 1:totalNumberScans) {
    mzPeakCountsIndex[i + 1] <- mzPeakCountsIndex[i] + scanTable$peaksCount[i]
  }
  ##
  spectraList <- lapply(1:totalNumberScans, function(i) {
    ##
    MZ <- netCDFList[["mass_values"]][(mzPeakCountsIndex[i] + 1):mzPeakCountsIndex[i + 1]]
    INT <- netCDFList[["intensity_values"]][(mzPeakCountsIndex[i] + 1):mzPeakCountsIndex[i + 1]]
    ##
    cbind(MZ, INT)
  })
  ##############################################################################
  ##############################################################################
  p2l <- list(scanTable, spectraList)
  ##
  names(p2l) <- c("scanTable", "spectraList")
  ##
  return(p2l)
}