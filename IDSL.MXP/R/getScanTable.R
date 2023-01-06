getScanTable <- function(xmlData, msFormat) {
  scanTable <- NA
  if (tolower(msFormat) == "mzml") {
    spectrumNodes <- xml_find_all(xmlData, '//d1:spectrum')
    xmlNSspectrumNodes <- xml_ns(spectrumNodes)
    ############################################################################
    positiveScanNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="positive scan"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    negativeScanNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="negative scan"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    centroidSpectrumNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="centroid spectrum"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    profileSpectrumNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="profile spectrum"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    basePeakMZNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="base peak m/z"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    basePeakIntensityNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="base peak intensity"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    totIonCurrentNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="total ion current"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    msLevelNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="ms level"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    lowMZNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="lowest observed m/z"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    highMZNodes <- xml_find_all(spectrumNodes, 'd1:cvParam[@name="highest observed m/z"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    retentionTimeNodes <- xml_find_all(spectrumNodes, 'd1:scanList/d1:scan/d1:cvParam[@name="scan start time"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    filterStringNodes <- xml_find_all(spectrumNodes, 'd1:scanList/d1:scan/d1:cvParam[@name="filter string"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    injectionTimeNodes <- xml_find_all(spectrumNodes, 'd1:scanList/d1:scan/d1:cvParam[@name="ion injection time"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    scanWindowLowerLimitNodes <- xml_find_all(spectrumNodes, 'd1:scanList/d1:scan/d1:scanWindowList/d1:scanWindow/d1:cvParam[@name="scan window lower limit"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    scanWindowUpperLimitNodes <- xml_find_all(spectrumNodes, 'd1:scanList/d1:scan/d1:scanWindowList/d1:scanWindow/d1:cvParam[@name="scan window upper limit"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    isolationWindowTargetMZNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:isolationWindow/d1:cvParam[@name="isolation window target m/z"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    isolationWindowLowerOffsetNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:isolationWindow/d1:cvParam[@name="isolation window lower offset"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    isolationWindowUpperOffsetNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:isolationWindow/d1:cvParam[@name="isolation window upper offset"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    precursorScanNumNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor', ns = xmlNSspectrumNodes, flatten = FALSE)
    precursorMZNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:selectedIonList/d1:selectedIon/d1:cvParam[@name="selected ion m/z"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    precursorChargeNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:selectedIonList/d1:selectedIon/d1:cvParam[@name="charge state"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    precursorIntensityNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:selectedIonList/d1:selectedIon/d1:cvParam[@name="peak intensity"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    collisionEnergyNodes <- xml_find_all(spectrumNodes, 'd1:precursorList/d1:precursor/d1:activation/d1:cvParam[@name="collision energy"]', ns = xmlNSspectrumNodes, flatten = FALSE)
    ############################################################################
    scanTable <- do.call(rbind, lapply(1:length(spectrumNodes), function(i) {
      S <- xml_attrs(spectrumNodes[[i]])
      ##
      seqNum <- as.numeric(S[1]) + 1
      ##
      spectrumId <- S[2]
      if (length(spectrumId) == 0) {
        spectrumId <- NA
      }
      ##
      x_acquisitionNum <- MXP_locate_regex(spectrumId, "=")
      if (!is.null(x_acquisitionNum)) {
        x_acquisitionNum <- x_acquisitionNum[nrow(x_acquisitionNum), 2]
        acquisitionNum <- substr(spectrumId, (x_acquisitionNum + 1), nchar(spectrumId))
      } else {
        acquisitionNum <- NULL
      }
      if (length(acquisitionNum) == 0) {
        acquisitionNum <- NA
      }
      ##
      peaksCount <- S[3]
      if (length(peaksCount) == 0) {
        peaksCount <- 0
      }
      ##########################################################################
      positiveScan <- xml_has_attr(positiveScanNodes[[i]], "value")
      ##
      negativeScan <- xml_has_attr(negativeScanNodes[[i]], "value")
      ##
      if (length(negativeScan) > 0) {
        if (negativeScan) {
          polarity <- 0
        }
      } else if (length(positiveScan) > 0) {
        if (positiveScan) {
          polarity <- 1
        }
      } else {
        polarity <- -1
      }
      ##########################################################################
      centroidSpectrum <- xml_has_attr(centroidSpectrumNodes[[i]], "value")
      ##
      profileSpectrum <- xml_has_attr(profileSpectrumNodes[[i]], "value")
      ##
      if (length(centroidSpectrum) > 0) {
        centroidSpectrum <- centroidSpectrum[1]
        if (centroidSpectrum) {
          centroided <- TRUE
        } else {
          centroided <- FALSE
        }
      } else if (length(profileSpectrum) > 0) {
        centroided <- FALSE
      } else {
        centroided <- FALSE
      }
      ##########################################################################
      basePeakMZ <- xml_attr(basePeakMZNodes[[i]], "value")
      if (length(basePeakMZ) == 0) {
        basePeakMZ <- NA
      } else {
        basePeakMZ <- basePeakMZ[1]
      }
      ##
      basePeakIntensity <- xml_attr(basePeakIntensityNodes[[i]], "value")
      if (length(basePeakIntensity) == 0) {
        basePeakIntensity <- NA
      } else {
        basePeakIntensity <- basePeakIntensity[1]
      }
      ##
      totIonCurrent <- xml_attr(totIonCurrentNodes[[i]], "value")
      if (length(totIonCurrent) == 0) {
        totIonCurrent <- NA
      } else {
        totIonCurrent <- totIonCurrent[1]
      }
      ##
      msLevel <- xml_attr(msLevelNodes[[i]], "value")
      if (length(msLevel) == 0) {
        msLevel <- 0
      } else {
        msLevel <- msLevel[1]
      }
      ##
      lowMZ <- xml_attr(lowMZNodes[[i]], "value")
      if (length(lowMZ) == 0) {
        lowMZ <- NA
      } else {
        lowMZ <- lowMZ[1]
      }
      ##
      highMZ <- xml_attr(highMZNodes[[i]], "value")
      if (length(highMZ) == 0) {
        highMZ <- NA
      } else {
        highMZ <- highMZ[1]
      }
      ##########################################################################
      retentionTime <- xml_attr(retentionTimeNodes[[i]], "value")
      if (length(retentionTime) == 0) {
        retentionTime <- NA
      } else {
        retentionTime <- retentionTime[1]
      }
      ##
      filterString <- xml_attr(filterStringNodes[[i]], "value")
      if (length(filterString) == 0) {
        filterString <- NA
      } else {
        filterString <- filterString[1]
      }
      ##
      injectionTime <- xml_attr(injectionTimeNodes[[i]], "value")
      if (length(injectionTime) == 0) {
        injectionTime <- NA
      } else {
        injectionTime <- injectionTime[1]
      }
      ##########################################################################
      scanWindowLowerLimit <- xml_attr(scanWindowLowerLimitNodes[[i]], "value")
      if (length(scanWindowLowerLimit) == 0) {
        scanWindowLowerLimit <- NA
      } else {
        scanWindowLowerLimit <- scanWindowLowerLimit[1]
      }
      ##
      scanWindowUpperLimit <- xml_attr(scanWindowUpperLimitNodes[[i]], "value")
      if (length(scanWindowUpperLimit) == 0) {
        scanWindowUpperLimit <- NA
      } else {
        scanWindowUpperLimit <- scanWindowUpperLimit[1]
      }
      ##########################################################################
      isolationWindowTargetMZ <- xml_attr(isolationWindowTargetMZNodes[[i]], "value")
      if (length(isolationWindowTargetMZ) == 0) {
        isolationWindowTargetMZ <- NA
      } else {
        isolationWindowTargetMZ <- isolationWindowTargetMZ[1]
      }
      ##
      isolationWindowLowerOffset <- xml_attr(isolationWindowLowerOffsetNodes[[i]], "value")
      if (length(isolationWindowLowerOffset) == 0) {
        isolationWindowLowerOffset <- NA
      } else {
        isolationWindowLowerOffset <- isolationWindowLowerOffset[1]
      }
      ##
      isolationWindowUpperOffset <- xml_attr(isolationWindowUpperOffsetNodes[[i]], "value")
      if (length(isolationWindowUpperOffset) == 0) {
        isolationWindowUpperOffset <- NA
      } else {
        isolationWindowUpperOffset <- isolationWindowUpperOffset[1]
      }
      ##########################################################################
      precursorScanNum <- xml_attr(precursorScanNumNodes[[i]], "spectrumRef")
      if (length(precursorScanNum) == 0) {
        precursorScanNum <- NA
      } else {
        if (!is.na(precursorScanNum)) {
          precursorScanNum <- precursorScanNum[1]
          x_scan <- MXP_locate_regex(precursorScanNum, "scan=")
          if (!is.null(x_scan)) {
            precursorScanNum <- substr(precursorScanNum, (x_scan[2] + 1), nchar(precursorScanNum))
          } else {
            precursorScanNum <- NA
          }
        }
      }
      ##
      precursorMZ <- xml_attr(precursorMZNodes[[i]], "value")
      if (length(precursorMZ) == 0) {
        precursorMZ <- NA
      } else {
        precursorMZ <- precursorMZ[1]
      }
      ##
      precursorCharge <- xml_attr(precursorChargeNodes[[i]], "value")
      if (length(precursorCharge) == 0) {
        precursorCharge <- NA
      } else {
        precursorCharge <- precursorCharge[1]
      }
      ##
      precursorIntensity <- xml_attr(precursorIntensityNodes[[i]], "value")
      if (length(precursorIntensity) == 0) {
        if (is.na(precursorMZ)) {
          precursorIntensity <- NA
        } else {
          precursorIntensity <- 0
        }
      } else {
        precursorIntensity <- precursorIntensity[1]
      }
      ##
      collisionEnergy <- xml_attr(collisionEnergyNodes[[i]], "value")
      if (length(collisionEnergy) == 0) {
        collisionEnergy <- NA
      } else {
        collisionEnergy <- collisionEnergy[1]
      }
      ##########################################################################
      c(seqNum, acquisitionNum, msLevel, polarity, peaksCount, totIonCurrent, retentionTime,
        basePeakMZ, basePeakIntensity, collisionEnergy, lowMZ, highMZ, precursorScanNum, precursorMZ,
        precursorCharge, precursorIntensity, injectionTime, filterString, spectrumId, centroided,
        isolationWindowTargetMZ, isolationWindowLowerOffset, isolationWindowUpperOffset,
        scanWindowLowerLimit, scanWindowUpperLimit)
    }))
    ############################################################################
    scanTable <- data.frame(scanTable)
    colnames(scanTable) <- c("seqNum", "acquisitionNum", "msLevel", "polarity", "peaksCount", "totIonCurrent", "retentionTime",
                             "basePeakMZ", "basePeakIntensity", "collisionEnergy", "lowMZ", "highMZ", "precursorScanNum", "precursorMZ",
                             "precursorCharge", "precursorIntensity", "injectionTime", "filterString", "spectrumId", "centroided",
                             "isolationWindowTargetMZ", "isolationWindowLowerOffset", "isolationWindowUpperOffset",
                             "scanWindowLowerLimit", "scanWindowUpperLimit")
    ##
    scanTable$acquisitionNum <- as.numeric(scanTable$acquisitionNum)
    ############################################################################
  } else if (tolower(msFormat) == "mzxml") {
    ##
    scanNodes <- xml_find_all(xmlData, '//d1:scan')
    ##
    scanNodeTable <- do.call(rbind, lapply(1:length(scanNodes), function(i) {
      S <- xml_attrs(scanNodes[[i]])
      ##########################################################################
      nameS <- tolower(names(S))
      ##########################################################################
      x_num <- which(nameS == "num")
      if (length(x_num) > 0) {
        seqNum <- S[x_num]
      } else {
        seqNum <- 0
      }
      ##########################################################################
      x_scanType <- which(nameS == "scantype")
      if (length(x_scanType) > 0) {
        scanType <- S[x_scanType]
      } else {
        scanType <- NA
      }
      ##########################################################################
      x_peaksCount <- which(nameS == "peakscount")
      if (length(x_peaksCount) > 0) {
        peaksCount <- S[x_peaksCount]
      } else {
        peaksCount <- 0
      }
      ##########################################################################
      x_polarity <- which(nameS == "polarity")
      if (length(x_polarity) > 0) {
        pol <- S[x_polarity]
        if (pol == "-" | grepl("N", pol, ignore.case = TRUE)) {
          polarity <- 0
        } else if (pol == "+" | grepl("P", pol, ignore.case = TRUE)) {
          polarity <- 1
        }
      } else {
        polarity <- -1
      }
      ##########################################################################
      x_centroided <- which(nameS == "centroided")
      if (length(x_centroided) > 0) {
        cent <- as.numeric(S[x_centroided])
        if (cent == 1) {
          centroided <- TRUE
        } else if (cent == 0) {
          centroided <- FALSE
        }
      } else {
        centroided <- FALSE
      }
      ##########################################################################
      x_basePeakMZ <- which(nameS == "basepeakmz")
      if (length(x_basePeakMZ) > 0) {
        basePeakMZ <- S[x_basePeakMZ]
      } else {
        basePeakMZ <- NA
      }
      ##########################################################################
      x_basePeakIntensity <- which(nameS == "basepeakintensity")
      if (length(x_basePeakIntensity) > 0) {
        basePeakIntensity <- S[x_basePeakIntensity]
      } else {
        basePeakIntensity <- NA
      }
      ##########################################################################
      x_totIonCurrent <- which(nameS == "totioncurrent")
      if (length(x_totIonCurrent) > 0) {
        totIonCurrent <- S[x_totIonCurrent]
      } else {
        totIonCurrent <- 0
      }
      ##########################################################################
      x_msLevel <- which(nameS == "mslevel")
      if (length(x_msLevel) > 0) {
        msLevel <- S[x_msLevel]
      } else {
        msLevel <- 0
      }
      ##########################################################################
      x_lowMZ <- which(nameS == "lowmz")
      if (length(x_lowMZ) > 0) {
        lowMZ <- S[x_lowMZ]
      } else {
        lowMZ <- NA
      }
      ##########################################################################
      x_highMZ <- which(nameS == "highmz")
      if (length(x_highMZ) > 0) {
        highMZ <- S[x_highMZ]
      } else {
        highMZ <- NA
      }
      ##########################################################################
      x_retentionTime <- which(nameS == "retentiontime")
      if (length(x_retentionTime) > 0) {
        retentionTime <- S[x_retentionTime]
      } else {
        retentionTime <- NA
      }
      ##########################################################################
      x_filterString <- which(nameS == "filterstring")
      if (length(x_filterString) > 0) {
        filterString <- S[x_filterString]
      } else {
        filterString <- NA
      }
      ##########################################################################
      x_injectionTime <- which(nameS == "injectiontime")
      if (length(x_injectionTime) > 0) {
        injectionTime <- S[x_injectionTime]
      } else {
        injectionTime <- NA
      }
      ##########################################################################
      x_scanWindowLowerLimit <- which(nameS == "scanwindowLowerlimit")
      if (length(x_scanWindowLowerLimit) > 0) {
        scanWindowLowerLimit <- S[x_scanWindowLowerLimit]
      } else {
        scanWindowLowerLimit <- NA
      }
      ##########################################################################
      x_scanWindowUpperLimit <- which(nameS == "scanwindowUpperlimit")
      if (length(x_scanWindowUpperLimit) > 0) {
        scanWindowUpperLimit <- S[x_scanWindowUpperLimit]
      } else {
        scanWindowUpperLimit <- NA
      }
      ##########################################################################
      x_isolationWindowTargetMZ <- which(nameS == "isolationwindowtargetmz")
      if (length(x_isolationWindowTargetMZ) > 0) {
        isolationWindowTargetMZ <- S[x_isolationWindowTargetMZ]
      } else {
        isolationWindowTargetMZ <- NA
      }
      ##########################################################################
      x_collisionEnergy <- which(nameS == "collisionenergy")
      if (length(x_collisionEnergy) > 0) {
        collisionEnergy <- S[x_collisionEnergy]
      } else {
        collisionEnergy <- NA
      }
      ##########################################################################
      c(seqNum, msLevel, polarity, peaksCount, totIonCurrent, retentionTime, basePeakMZ,
        basePeakIntensity, collisionEnergy, lowMZ, highMZ, injectionTime, filterString,
        scanType, centroided, isolationWindowTargetMZ, scanWindowLowerLimit, scanWindowUpperLimit)
    }))
    ##
    scanNodeTable <- data.frame(scanNodeTable)
    colnames(scanNodeTable) <- c("seqNum", "msLevel", "polarity", "peaksCount", "totIonCurrent", "retentionTime", "basePeakMZ",
                                 "basePeakIntensity", "collisionEnergy", "lowMZ", "highMZ", "injectionTime", "filterString",
                                 "scanType", "centroided", "isolationWindowTargetMZ", "scanWindowLowerLimit", "scanWindowUpperLimit")
    ############################################################################
    mslevel <- as.numeric(scanNodeTable$msLevel)
    ##
    precursorMatrix <- matrix(rep(NA, 7*length(mslevel)), ncol = 7)
    ##
    precursorNodes <- xml_find_all(xmlData, '//d1:precursorMz', ns = xml_ns(xmlData), flatten = FALSE)
    L_precursorNodes <- length(xml_attrs(precursorNodes))
    if (L_precursorNodes == 0) {
      precursorNodes <- xml_find_all(xmlData, '//d1:precursorMZ', ns = xml_ns(xmlData), flatten = FALSE)
      L_precursorNodes <- length(xml_attrs(precursorNodes))
      if (L_precursorNodes == 0) {
        precursorNodes <- xml_find_all(xmlData, '//d1:precursormz', ns = xml_ns(xmlData), flatten = FALSE)
        L_precursorNodes <- length(xml_attrs(precursorNodes))
      }
    }
    ############################################################################
    if (L_precursorNodes > 0) {
      counterNodes <- 1
      counter <- 0
      for (i in mslevel) {
        counter <- counter + 1
        if (i > 1) {
          ##
          precursorMatrix[counter, 1] <- xml_text(precursorNodes[[counterNodes]]) 
          ######################################################################
          precursorVec <- xml_attrs(precursorNodes[[counterNodes]])
          ######################################################################
          nameP <- tolower(names(precursorVec))
          ######################################################################
          x_precursorScanNum <- which(nameP == "precursorscannum")
          if (length(x_precursorScanNum) > 0) {
            precursorMatrix[counter, 2] <- precursorVec[x_precursorScanNum]
          }
          ######################################################################
          x_precursorIntensity <- which(nameP == "precursorintensity")
          if (length(x_precursorIntensity) > 0) {
            precursorMatrix[counter, 3] <- precursorVec[x_precursorIntensity]
          }
          ######################################################################    
          x_precursorCharge <- which(nameP == "precursorcharge")
          if (length(x_precursorCharge) > 0) {
            precursorMatrix[counter, 4] <- precursorVec[x_precursorCharge]
          } 
          ######################################################################
          x_activationMethod <- which(nameP == "activationmethod")
          if (length(x_activationMethod) > 0) {
            precursorMatrix[counter, 5] <- precursorVec[x_activationMethod]
          }
          ######################################################################
          x_windowWideness <- which(nameP == "windowwideness")
          if (length(x_windowWideness) > 0) {
            isolationWindowLowerOffset <- as.numeric(precursorVec[x_windowWideness])/2
            ##
            precursorMatrix[counter, 6] <- isolationWindowLowerOffset
            precursorMatrix[counter, 7] <- isolationWindowLowerOffset
          }
          ######################################################################
          counterNodes <- counterNodes + i - 1
        }
      }
    }
    ##
    precursorMatrix <- data.frame(precursorMatrix)
    colnames(precursorMatrix) <- c("precursorMZ", "precursorScanNum", "precursorIntensity", "precursorCharge", "activationMethod", "isolationWindowLowerOffset", "isolationWindowUpperOffset")
    ##
    scanTable <- cbind(scanNodeTable, precursorMatrix)
    ##
    scanTable$retentionTime <- as.numeric(gsub("[a-zA-Z]", "", scanTable$retentionTime, ignore.case = TRUE))/60
    ##
  } else {
    stop("The MSfile is not consistent with the IDSL.MXP package!")
  }
  ##############################################################################
  scanTable$seqNum <- as.numeric(scanTable$seqNum)
  scanTable$msLevel <- as.numeric(scanTable$msLevel)
  scanTable$polarity <- as.numeric(scanTable$polarity)
  scanTable$peaksCount <- as.numeric(scanTable$peaksCount)
  scanTable$totIonCurrent <- as.numeric(scanTable$totIonCurrent)
  scanTable$retentionTime <- as.numeric(scanTable$retentionTime)
  scanTable$basePeakMZ <- as.numeric(scanTable$basePeakMZ)
  scanTable$basePeakIntensity <- as.numeric(scanTable$basePeakIntensity)
  scanTable$lowMZ <- as.numeric(scanTable$lowMZ)
  scanTable$highMZ <- as.numeric(scanTable$highMZ)
  scanTable$precursorScanNum <- as.numeric(scanTable$precursorScanNum)
  scanTable$precursorMZ <- as.numeric(scanTable$precursorMZ)
  scanTable$precursorIntensity <- as.numeric(scanTable$precursorIntensity)
  scanTable$injectionTime <- as.numeric(scanTable$injectionTime)
  scanTable$isolationWindowTargetMZ <- as.numeric(scanTable$isolationWindowTargetMZ)
  scanTable$isolationWindowLowerOffset <- as.numeric(scanTable$isolationWindowLowerOffset)
  scanTable$isolationWindowUpperOffset <- as.numeric(scanTable$isolationWindowUpperOffset)
  scanTable$scanWindowLowerLimit <- as.numeric(scanTable$scanWindowLowerLimit)
  scanTable$scanWindowUpperLimit <- as.numeric(scanTable$scanWindowUpperLimit)
  ##############################################################################
  return(scanTable)
}