readNanoStringGeomxSet <-
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

  # remove any zero reads in Dcc files
  zeroRead <- which(sapply(seq_len(length(data)), 
                           function(x) nrow(data[[x]]$Code_Summary))==0)
  if(length(zeroRead) > 0){
    warning("The following DCC files are removed: ", names(data)[zeroRead])
    data <- data[-zeroRead]
    dccFiles <- dccFiles[-zeroRead]
  }
  
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
    missingPhenoCount <- sum(!(colnames(assay) %in% pheno[[j]]))
    pheno[[j]] <- paste0(pheno[[j]], ".dcc")
    rownames(pheno) <- pheno[[j]]
    pheno[[j]] <- NULL
    pheno <- pheno[names(assay), , drop = FALSE]
    if (missingPhenoCount != 0L) {
      rownames(pheno) <- colnames(assay)
      warning(sprintf("Column `phenoDataDccColName` in `phenoDataFile` is missing %d of %d Samples",
                      missingPhenoCount, ncol(assay)))
    }
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
    
    pkcData <- as.data.frame(pkcData)
    rownames(pkcData) <- pkcData[["RTS_ID"]]
  }
  
  for (index in seq_len(length(data))) {
    countMat <- data[[index]]$Code_Summary
    countMat <- countMat[which(rownames(countMat) %in% rownames(pkcData)), , drop = FALSE]
    countMat$GeneName <- pkcData$Gene[match(countMat$RNAID, 
                                        pkcData$RTS_ID)]
    
    countMat$Pool <- pkcData$Module[match(countMat$RNAID, 
                                          pkcData$RTS_ID)]
    data[[index]]$Code_Summary <- countMat
  }

  probe_assay <- lapply(seq_len(length(data)), function(x)
    data.frame(data[[x]][["Code_Summary"]],
               Sample_ID = names(data)[x]))
  probe_assay <- do.call(rbind, probe_assay)
  gene_assay <- reshape2::dcast(probe_assay, GeneName + Pool ~ Sample_ID, 
                      value.var = 'Count', fun.aggregate = ngeoMean, fill = 1)
  
  if ( length(unique(gene_assay$GeneName)) != nrow(gene_assay) ) {
    duplicatedGenes <- gene_assay$GeneName[which(duplicated(gene_assay$GeneName))]
    warning(sprintf('Some genes are listed in multiple pools including %s', 
            paste0(duplicatedGenes, collapse = ",")))
    for (duplicatedGene in duplicatedGenes ){
      gene_assay$GeneName[which(gene_assay$GeneName == duplicatedGene)] <- 
        sapply(which(gene_assay$GeneName == duplicatedGene), 
                      function(index) paste0(gene_assay[index, c("GeneName", "Pool")], 
                                             collapse = "_"))
    }
  }
  
  rownames(gene_assay) <- gene_assay[, "GeneName"]
  assay <- as.matrix(gene_assay[, -seq_len(2)])
  
  # Create featureData
  feature <- gene_assay[, "GeneName", drop = FALSE]
  rownames(feature) <- feature[["GeneName"]]
  
  feature <- AnnotatedDataFrame(feature,
                                dimLabels = c("featureNames", "featureColumns"))

  # Create experimentData
  experimentList<- lapply(experimentDataColNames, 
                            function(experimentDataColName) 
                              unique(S4Vectors::na.omit(pheno@data[[experimentDataColName]])))
  names(experimentList) <- experimentDataColNames
  
  experiment <- Biobase::MIAME(name = "", 
                      other = c(experimentList, pkcHeader))
  
  # Create annotation
  annotation <- sort(sapply(strsplit(pkcFiles, "/"), function(x) x[length(x)]))
  if(!identical(annotation, paste0(sort(unique(probe_assay[['Pool']])), ".pkc"))) {
    stop("Name mismatch between pool and PKC files")
  }

  # Create protocolData
  protocol <-
    do.call(rbind,
            lapply(seq_along(dccFiles), function(i) {
              cbind(data[[i]][["Header"]], data[[i]][["Scan_Attributes"]],
                    data[[i]][["NGS_Processing_Attributes"]])
            }))
  
  protocol <- data.frame(protocol, 
                         pheno@data[, which(colnames(pheno@data) %in% protocolDataColNames)])
  
  pheno <- pheno[, setdiff(colnames(pheno@data), 
                           c(protocolDataColNames, experimentDataColNames))]
  
  annot_labelDescription <-  data.frame(labelDescription =
                                           rep(NA_character_, length(protocolDataColNames)),
                                         row.names = protocolDataColNames,
                                         stringsAsFactors = FALSE)
  protocol <- AnnotatedDataFrame(protocol,
                                 rbind(.dccMetadata[["protocolData"]], 
                                       annot_labelDescription),
                                 dimLabels = c("sampleNames", "sampleColumns"))

  # Create NanoStringGeomxSet
  return( NanoStringGeomxSet(assayData = assay,
                             phenoData = pheno,
                             featureData = feature,
                             experimentData = experiment,
                             annotation = annotation,
                             protocolData = protocol,
                             check = FALSE) )
}
