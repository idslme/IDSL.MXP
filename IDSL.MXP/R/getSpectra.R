getSpectra <- function(xmlData, msFormat) {
  spectraList <- NA
  if (tolower(msFormat) == "mzml") {
    spectrumNodes1 <- xml_find_first(xmlData, '//d1:spectrum')
    ############################################################################
    firstNode <- xml_find_first(spectrumNodes1, "//d1:spectrum")
    ##
    BitNode1000523 <- xml_find_first(firstNode, '//d1:cvParam[@accession="MS:1000523"]')
    BitType1000523 <- xml_attr(BitNode1000523, "name")
    precision1000523 <- as.numeric(gsub("-bit float", "", BitType1000523))/8
    ##
    BitNode1000521 <- xml_find_first(firstNode, '//d1:cvParam[@accession="MS:1000521"]')
    BitType1000521 <- xml_attr(BitNode1000521, "name")
    precision1000521 <- as.numeric(gsub("-bit float", "", BitType1000521))/8
    #######
    if (!is.na(precision1000521) & is.na(precision1000523)) {
      mzPrecision <- precision1000521
      intPrecision <- precision1000521
    } else if (is.na(precision1000521) & !is.na(precision1000523)) {
      mzPrecision <- precision1000523
      intPrecision <- precision1000523
    } else if (!is.na(precision1000521) & !is.na(precision1000523)) {
      mzPrecision <- precision1000523
      intPrecision <- precision1000521
    }
    ############################################################################
    compressionTypeFirst <- xml_find_first(firstNode, '//d1:cvParam[@accession="MS:1000574"]|//d1:cvParam[@accession="MS:1000576"]')
    ##
    compressionType <- "none"
    #
    if (length(compressionTypeFirst) > 0) {
      compressionTypeFirst <- xml_attr(compressionTypeFirst, "name")
      if (length(compressionTypeFirst) > 0) {
        if ((compressionTypeFirst == "zlib") | (compressionTypeFirst == "zlib compression")) {
          compressionType <- "gzip"
        }
      }
    }
    ############################################################################
    ## Thanks to the RaMS package
    deCompresser_mzML <- function(binary, compressionType, precision) {
      decodedBinary <- base64decode(binary)
      rawBinary <- as.raw(decodedBinary)
      decompressedBinary <- memDecompress(rawBinary, type = compressionType)
      deComp <- readBin(decompressedBinary, what = "double", n = length(decompressedBinary)/precision, size = precision)
      return(deComp)
    }
    ############################################################################
    binaryNodes <- xml_find_all(xmlData, '//d1:spectrum/d1:binaryDataArrayList/d1:binaryDataArray/d1:binary', ns = xml_ns(xmlData), flatten = FALSE)
    spectraList <- lapply(seq(1, (length(binaryNodes) - 1), by = 2), function(i) {
      binary_mz <- xml_text(binaryNodes[[i]])
      if (binary_mz != "") {
        MZ <- deCompresser_mzML(binary_mz, compressionType, mzPrecision)
      } else {
        MZ <- matrix(0, nrow = 1)
      }
      ##
      binary_int <- xml_text(binaryNodes[[(i + 1)]])
      if (binary_int != "") {
        INT <- deCompresser_mzML(binary_int, compressionType, intPrecision)
      } else {
        INT <- matrix(0, nrow = 1)
      }
      ##
      cbind(MZ, INT)
    })
    ##
    ############################################################################
    ##
  } else if (tolower(msFormat) == "mzxml") {
    peaksNodes <- xml_find_all(xmlData, '//d1:scan/d1:peaks')
    compressionOrder <- xml_attrs(peaksNodes[[1]])
    c_Order <- names(compressionOrder)
    compressionOrder <- data.frame(matrix(compressionOrder, nrow = 1))
    colnames(compressionOrder) <- c_Order
    ##
    compressionTypeFirst <- compressionOrder$compressionType
    #
    compressionType <- "none"
    #
    if (length(compressionTypeFirst) > 0) {
      if ((compressionTypeFirst == "zlib") | (compressionTypeFirst == "zlib compression")) {
        compressionType <- "gzip"
      }
    }
    ##
    precision <- as.numeric(compressionOrder$precision)/8
    ##
    endianTypeFirst <- compressionOrder$byteOrder
    #
    endianType <- "swap"
    #
    if (length(endianTypeFirst) > 0) {
      if ((endianTypeFirst == "network")) {
        endianType <- "big"
      }
    }
    ############################################################################
    ## Thanks to the RaMS package
    deCompresser_mzXML <- function(binary, compressionType, precision, endianType) {
      decodedBinary <- base64decode(binary)
      rawBinary <- as.raw(decodedBinary)
      decompressedBinary <- memDecompress(rawBinary, type = compressionType)
      deComp <- readBin(decompressedBinary, what = "numeric", n = length(decompressedBinary)/precision, size = precision, endian = endianType)
      return(deComp)
    }
    ##
    spectraList <- lapply(1:length(peaksNodes), function(i) {
      binary <- xml_text(peaksNodes[[i]], '[@contentType="m/z-int"]')
      if (binary != "") {
        MZ_INT <- deCompresser_mzXML(binary, compressionType, precision, endianType)
        matrix(MZ_INT, ncol = 2, byrow = TRUE)
      } else {
        matrix(c(0, 0), ncol = 2)
      }
    })
    ##
  } else {
    stop("The MSfile is not consistent with the IDSL.MXP package!")
  }
  return(spectraList)
}