readNanoStringGeoMxSet <-
function(dccFiles,
         pkcFiles,
         phenoDataFile,
         phenoDataSheet,
         phenoDataDccColName = "Sample_ID",
         phenoDataColPrefix = "",
         protocolDataColNames = NULL,
         experimentDataColNames = NULL,
         configFile = NULL,
         analyte = "RNA",
         defaultPKCVersions = NULL,
         ...)
{
  # check inputs
  if (!(sum(grepl("\\.dcc$",dccFiles)) == length(dccFiles) && length(dccFiles) > 0L)){
    stop("Specify valid dcc files." )
  }
  if (!(sum(grepl("\\.pkc$",pkcFiles)) == length(pkcFiles) && length(pkcFiles) > 0L)){
    stop( "Specify valid PKC files." )
  }
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
    pheno <- readxl::read_xlsx(phenoDataFile, col_names = TRUE, sheet = phenoDataSheet, ...)
    pheno <- data.frame(pheno, stringsAsFactors = FALSE, check.names = FALSE)
    j <- colnames(pheno)[colnames(pheno) == phenoDataDccColName]
    if (length(j) == 0L){
      stop("Column `phenoDataDccColName` not found in `phenoDataFile`")
    } else if (length(j) > 1L){
      stop("Multiple columns in `phenoDataFile` match `phenoDataDccColName`")
    }
    # check protocolDataColNames
    if (!(all(protocolDataColNames %in% colnames(pheno))) &
          !(is.null(protocolDataColNames))) {
      stop("Columns specified in `protocolDataColNames` are not found in `phenoDataFile`")
    }
    # check experimentDataColNames
    if (!(all(experimentDataColNames %in% colnames(pheno))) &
        !(is.null(experimentDataColNames))) {
      stop("Columns specified in `experimentDataColNames` are not found in `phenoDataFile`")
    }
    # add ".dcc" to the filenames if there is none
    pheno[[j]] <- ifelse(grepl(".dcc", pheno[[j]]), paste0(pheno[[j]]),
                         paste0(pheno[[j]], ".dcc"))
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
    
    if("area" %in% colnames(pheno)){
        pheno[["area"]] <- as.numeric(pheno[["area"]])
    }
    if("nuclei" %in% colnames(pheno)){
        pheno[["nuclei"]] <- as.numeric(pheno[["nuclei"]])
    }
    
    pheno <- Biobase::AnnotatedDataFrame(pheno,
                                dimLabels = c("sampleNames", "sampleColumns"))
  }

  #stopifnot(all(sapply(feature, function(x) identical(feature[[1L]], x))))
  if (is.null(pkcFiles)) {
    stop("Please specify an input for pkcFiles")
  } else if (!is.null(pkcFiles)) {
    pkcData <- readPKCFile(pkcFiles, default_pkc_vers=defaultPKCVersions)

    pkcHeader <- S4Vectors::metadata(pkcData)
    # pkcHeader[["PKCFileDate"]] <- as.character(pkcHeader[["PKCFileDate"]])

    pkcData$RTS_ID <- gsub("RNA", "RTS00", pkcData$RTS_ID)

    pkcData <- as.data.frame(pkcData)
    rownames(pkcData) <- pkcData[["RTS_ID"]]
  }

  probeAssay <- lapply(names(data), function(x)
    data.frame(data[[x]][["Code_Summary"]],
               Sample_ID = x))
  probeAssay <- do.call(rbind, probeAssay)
  # check for missing probes in PKC that are in assay
  missingProbes <- setdiff(unique(probeAssay[["RTS_ID"]]), rownames(pkcData))
  if (length(missingProbes) > 0L){
    warning("Not all probes are found within PKC probe metadata.",
    " The following probes are ignored from analysis",
    " and were most likely removed from metadata while",
    " resolving multiple module PKC version conflicts.\n",
    paste(missingProbes, sep=", "))
  }
  
  if(!is.null(configFile)){
    pkcProbes <- compareToConfig(config = configFile, pkcProbes = pkcData, pkcHeader = pkcHeader)
    
    pkcData <- pkcProbes$pkcData
    pkcHeader <- pkcProbes$pkcHeader
    
    pkcFiles <- paste0(unique(pkcData$Module), ".pkc")
  }
  
  if(tolower(analyte) == "protein"){
    analyte <- "Protein"
  }else if(tolower(analyte) == "rna"){
    analyte <- "RNA"
  }
  
  pkcs <- names(which(tolower(pkcHeader$AnalyteType) == tolower(analyte)))
  
  if(length(pkcs) == 0){
    stop(paste("Given analyte is not valid: options =", paste(unique(pkcHeader$AnalyteType), collapse = ", ")))
  }
  
  pkcFiles <- paste0(pkcs, ".pkc")
  
  pkcData <- pkcData[pkcData$Module %in% pkcs,]
  for(i in names(pkcHeader)){
    pkcHeader[[i]] <- pkcHeader[[i]][names(pkcHeader[[i]]) %in% pkcs]
  }
  
  probeAssay <- probeAssay[probeAssay[["RTS_ID"]] %in% pkcData[["RTS_ID"]],]
  
  zeroProbes <- setdiff(rownames(pkcData), unique(probeAssay[["RTS_ID"]]))
  zeroProbeAssay <- data.frame(RTS_ID=pkcData[zeroProbes, "RTS_ID"],
    Count=rep(0, length(zeroProbes)),
    Sample_ID=rep(probeAssay[1, "Sample_ID"], length(zeroProbes)))
  probeAssay <- rbind(probeAssay, zeroProbeAssay)
  probeAssay[["Module"]] <- pkcData[probeAssay[["RTS_ID"]], "Module"]
  probeAssay <- reshape2::dcast(probeAssay, RTS_ID + Module ~ Sample_ID,
      value.var="Count", fill=0)
  rownames(probeAssay) <- probeAssay[, "RTS_ID"]
  
  if(!all(names(data) %in% names(probeAssay))){
    missingAnalyteDCC <- names(data)[which(!names(data) %in% names(probeAssay))]
    warning("The following DCC files had no counts and will be excluded from the GeoMxSet object: ",
            paste0(missingAnalyteDCC, sep=", "))
    
    data <- data[which(names(data) %in% names(probeAssay))]
    
    pheno <- pheno[which(names(data) %in% 
                           sampleNames(pheno)),]
  }
  
  assay <- as.matrix(probeAssay[, names(data)])

  # Create featureData
  feature <- pkcData[rownames(assay), , drop = FALSE]
  # change the colnames of feature data to match with dimLabels
  colnames(feature)[which(colnames(feature)=="Target")] <- "TargetName"
  feature <- AnnotatedDataFrame(feature,
                                dimLabels = c("featureNames", "featureColumns"))

    # Create experimentData
    if (!(is.null(experimentDataColNames))) {
        experimentList <- 
            lapply(experimentDataColNames,
                   function(experimentDataColName) {
                       unique(S4Vectors::na.omit(pheno@data[[experimentDataColName]]))})
        names(experimentList) <- experimentDataColNames
        experiment <- 
            Biobase::MIAME(name = "", 
                           other = c(experimentList, 
                                     pkcHeader, 
                                     list(shiftedByOne=FALSE)))
    } else {
        experiment <- 
            Biobase::MIAME(name = "",
                           other = c(pkcHeader, 
                                     list(shiftedByOne=FALSE)))
    }

    # Create annotation
    annotation <- sort(sapply(strsplit(pkcFiles, "/"), function(x) x[length(x)]))
    if(!identical(annotation, paste0(sort(unique(probeAssay[["Module"]])), ".pkc"))) {
        stop("Name mismatch between pool and PKC files")
    }

  # Create protocolData
  protocol <-
    do.call(dplyr::bind_rows,
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
  GxT <- NanoStringGeoMxSet(assayData = assay,
                            phenoData = pheno,
                            featureData = feature,
                            experimentData = experiment,
                            annotation = annotation,
                            protocolData = protocol,
                            check = FALSE, 
                            dimLabels = c("RTS_ID", "SampleID"),
                            analyte = analyte)
  
  if(analyte(GxT) == "Protein"){
    GxT <- suppressWarnings(aggregateCounts(GxT))
  }
  
  return(GxT)
}

#' Compare given PKC probes to probes in config file
#' 
#' Check if extra PKCs are given based on probes in config file
#' 
#' @param config file path to config file
#' @param pkcProbes probe information from readPKCFile
#' @param pkcHeader pkc metadata from readPKCFile
#' 

compareToConfig <- function(config, pkcProbes, pkcHeader){
  if(!file.exists(config)){
    stop("Given config file does not exist")
  }
  
  config <- readLines(config)
  targets <- config[(which(config == "[Targets]")+1):length(config)]
  targets <- str_split(string = targets, pattern = " = ", simplify = T)[,1]
  
  extraProbes <- which(!pkcProbes$RTS_ID %in% targets)
  
  if(length(extraProbes) > 0){
    extraPKCs <- unique(pkcProbes$Module[extraProbes])
    
    for(pkc in extraPKCs){
      if(all(which(pkcProbes$Module == pkc) %in% extraProbes) == FALSE){
        extraPKCs <- extraPKCs[-which(extraPKCs == pkc)]
      }
    }
    
    if(length(extraPKCs) > 0){
      warning(paste0("Extra PKC files were provided and will be removed from further analysis: ", paste(extraPKCs, collapse = ", ")),
              immediate. = TRUE)
      
      pkcs <- unique(pkcProbes$Module[!pkcProbes$Module %in% extraPKCs])
      
      pkcProbes <- pkcProbes[pkcProbes$Module %in% pkcs,]
      
      for(i in names(pkcHeader)){
        pkcHeader[[i]] <- pkcHeader[[i]][names(pkcHeader[[i]]) %in% pkcs]
      }
    }
  }
  
  return(list(pkcData=pkcProbes, pkcHeader=pkcHeader))
}


# analyteSubset <- function(object, analyte){
#   if(length(unique(fData(object)$AnalyteType)) > 1){
#     if(!tolower(analyte) %in% tolower(unique(fData(object)$AnalyteType))){
#       stop(paste("Given analyte is not in dataset; options:", paste(unique(fData(object)$AnalyteType), collapse = ", ")))
#     }
#     
#     objects <- list()
#     pkcs <- object@annotation
#     
#     object <- subset(object, subset = tolower(AnalyteType) == tolower(analyte))
#     object@annotation <- pkcs[pkcs %in% paste0(unique(fData(object)$Module), ".pkc")]
#     
#     return(object)
#   }else{
#     warning("Only one analyte present, no subsetting neccesary")
#   }
# }
