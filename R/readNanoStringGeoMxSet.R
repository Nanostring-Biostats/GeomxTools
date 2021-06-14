readNanoStringGeoMxSet <-
function(dccFiles,
         pkcFiles,
         phenoDataFile,
         phenoDataSheet,
         phenoDataDccColName = "Sample_ID",
         phenoDataColPrefix = "",
         protocolDataColNames = c("slide name"),
         experimentDataColNames = c("panel"))
{
  # Read data rccFiles
  data <- structure(lapply(dccFiles, readDccFile), names = basename(dccFiles))

  # Create assayData
  assay <- lapply(data, function(x)
    structure(x[["Code_Summary"]][["Count"]],
              names = rownames(x[["Code_Summary"]])))

  # Create phenoData
  if (is.null(phenoDataFile)) {
    stop("Please specify an input for phenoDataFile.")
  } else {
    pheno <- readxl::read_xlsx(phenoDataFile, col_names = TRUE, sheet = phenoDataSheet)
    pheno <- data.frame(pheno, stringsAsFactors = FALSE, check.names = FALSE)
    j <- grep(phenoDataDccColName, colnames(pheno), ignore.case = TRUE)
    if (length(j) == 0L){
      stop("Column `phenoDataDccColName` not found in `phenoDataFile`")
    } else if (length(j) > 1L){
      stop("Multiple columns in `phenoDataFile` match `phenoDataDccColName`")
    }
    # check protocolDataColNames
    if (!(all(protocolDataColNames %in% colnames(pheno)))){
      stop("Columns specified in `protocolDataColNames` are not found in `phenoDataFile`")
    }
    pheno[[j]] <- paste0(pheno[[j]], ".dcc")
    if ("slide name" %in% colnames(pheno)) {
        ntcs <- which(tolower(pheno[["slide name"]]) == "no template control")
        if (length(ntcs) > 0) {
            ntcData <- lapply(seq_along(ntcs), function(x) {
                ntcID <- pheno[ntcs[x], j]
                if(!is.na(ntcs[x + 1L])) {
                    ntcNames <- rep(ntcID, ntcs[x + 1L] - ntcs[x])
                    ntcCounts <-
                        rep(sum(assay[[ntcID]]), ntcs[x + 1L] - ntcs[x])
                    ntcDF <- data.frame("NTC_ID"=ntcNames, "NTC"=ntcCounts)
                } else {
                    ntcNames <- rep(ntcID, dim(pheno)[1L] - ntcs[x] + 1L)
                    ntcCounts <-
                        rep(sum(assay[[ntcID]]), dim(pheno)[1L] - ntcs[x] + 1L)
                    ntcDF <- data.frame("NTC_ID"=ntcNames, "NTC"=ntcCounts)
                }
                return(ntcDF)
            })
            if (length(ntcs) > 1L) {
                ntcData <- do.call(rbind, ntcData)
            } else {
                ntcData <- ntcData[[1L]]
            }
            pheno <- cbind(pheno, ntcData)
            pheno <- pheno[!rownames(pheno) %in% ntcs, ]
            assay <- assay[!names(assay) %in% unique(pheno[["NTC_ID"]])]
            data <- data[!names(data) %in% unique(pheno[["NTC_ID"]])]
            protocolDataColNames <- c(protocolDataColNames, "NTC_ID", "NTC")
        }
    }
    rownames(pheno) <- pheno[[j]]
    zeroReads <- names(which(lapply(assay, length) == 0L))
    if (length(zeroReads) > 0L) {
        warning("The following DCC files had no counts: ",
                paste0(zeroReads, sep=", "),
                "These will be excluded from the GeoMxSet object.")
        pheno <- pheno[!rownames(pheno) %in% zeroReads, ]
        assay <- assay[!names(assay) %in% zeroReads]
        data <- data[!names(data) %in% zeroReads]
    }
    missingDCCFiles <- pheno[[j]][!pheno[[j]] %in% names(assay)]
    missingPhenoData <- names(assay)[!names(assay) %in% pheno[[j]]]
    assay <- assay[names(assay) %in% pheno[[j]]]
    data <- data[names(data) %in% pheno[[j]]]
    pheno <- pheno[names(assay), , drop = FALSE]
    if (length(missingDCCFiles) > 0L) {
      warning("DCC files missing for the following: ",
              paste0(missingDCCFiles, sep=", "),
              "These will be excluded from the GeoMxSet object.")
    }
    if (length(missingPhenoData) > 0L) {
      warning("Annotations missing for the following: ",
              paste0(missingPhenoData, sep=", "),
              "These will be excluded from the GeoMxSet object.")
    }
    pheno[[j]] <- NULL
    if (phenoDataColPrefix != "") {
      colnames(pheno) <- paste0(phenoDataColPrefix, colnames(pheno))
      protocolDataColNames <- paste0(phenoDataColPrefix, protocolDataColNames)
    }
    pheno <- Biobase::AnnotatedDataFrame(pheno,
                                dimLabels = c("sampleNames", "sampleColumns"))
  }

  #stopifnot(all(sapply(feature, function(x) identical(feature[[1L]], x))))
  if (is.null(pkcFiles)) {
    stop("Please specify an input for pkcFiles")
  } else if (!is.null(pkcFiles)) {
    pkcData <- readPKCFile(pkcFiles)

    pkcHeader <- S4Vectors::metadata(pkcData)
    pkcHeader[["PKCFileDate"]] <- as.character(pkcHeader[["PKCFileDate"]])

    pkcData$RTS_ID <- gsub("RNA", "RTS00", pkcData$RTS_ID)

    pkcData <- as.data.frame(pkcData)
    rownames(pkcData) <- pkcData[["RTS_ID"]]
  }

  probeAssay <- lapply(names(data), function(x)
    data.frame(data[[x]][["Code_Summary"]],
               Sample_ID = x))
  probeAssay <- do.call(rbind, probeAssay)
  zeroProbes <- setdiff(rownames(pkcData), unique(probeAssay[["RTS_ID"]]))
  zeroProbeAssay <- data.frame(RTS_ID=pkcData[zeroProbes, "RTS_ID"],
    Count=rep(0, length(zeroProbes)),
    Sample_ID=rep(probeAssay[1, "Sample_ID"], length(zeroProbes)))
  probeAssay <- rbind(probeAssay, zeroProbeAssay)
  probeAssay[["Module"]] <- pkcData[probeAssay[["RTS_ID"]], "Module"]
  probeAssay <- reshape2::dcast(probeAssay, RTS_ID + Module ~ Sample_ID,
      value.var="Count", fill=0)
  rownames(probeAssay) <- probeAssay[, "RTS_ID"]
  assay <- as.matrix(probeAssay[, names(data)])

  # Create featureData
  feature <- pkcData[rownames(assay), , drop = FALSE]
  # change the colnames of feature data to match with dimLables
  colnames(feature)[which(colnames(feature)=="Target")] <- "TargetName"
  feature <- AnnotatedDataFrame(feature,
                                dimLabels = c("featureNames", "featureColumns"))

  # Create experimentData
  experimentList<- lapply(experimentDataColNames,
                            function(experimentDataColName)
                              unique(S4Vectors::na.omit(pheno@data[[experimentDataColName]])))
  names(experimentList) <- experimentDataColNames

  experiment <- Biobase::MIAME(name = "",
                      other = c(experimentList, pkcHeader, list(shiftedByOne=FALSE)))

  # Create annotation
  annotation <- sort(sapply(strsplit(pkcFiles, "/"), function(x) x[length(x)]))
  if(!identical(annotation, paste0(sort(unique(probeAssay[["Module"]])), ".pkc"))) {
    stop("Name mismatch between pool and PKC files")
  }

  # Create protocolData
  protocol <-
    do.call(rbind,
            lapply(names(data), function(i) {
              cbind(data[[i]][["Header"]], data[[i]][["Scan_Attributes"]],
                    data[[i]][["NGS_Processing_Attributes"]])
            }))

  protocol <- data.frame(protocol,
                         pheno@data[, which(colnames(pheno@data) %in% protocolDataColNames)],
                         check.names = FALSE)

  pheno <- pheno[, setdiff(colnames(pheno@data),
                           c(protocolDataColNames, experimentDataColNames))]

  annot_labelDescription <-
    data.frame(
      labelDescription=rep(NA_character_, length(protocolDataColNames) + 1L),
      row.names = c(protocolDataColNames, "DeduplicatedReads"),
      stringsAsFactors = FALSE)
  protocol <-
    protocol[, names(protocol) %in% c(rownames(.dccMetadata[["protocolData"]]),
                                      rownames(annot_labelDescription))]
  protocol <- AnnotatedDataFrame(protocol,
                                 rbind(.dccMetadata[["protocolData"]],
                                       annot_labelDescription),
                                 dimLabels = c("sampleNames", "sampleColumns"))

  # Create NanoStringGeoMxSet
  return( NanoStringGeoMxSet(assayData = assay,
                             phenoData = pheno,
                             featureData = feature,
                             experimentData = experiment,
                             annotation = annotation,
                             protocolData = protocol,
                             check = FALSE) )
}
