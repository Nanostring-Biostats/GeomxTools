HOUSEKEEPERS <- c(
  "C1orf43", "GPI", "OAZ1", "POLR2A", "PSMB2", "RAB7A",
  "SDHA", "SNRPD3", "TBC1D10B", "TPM4", "TUBB", "UBB"
)
#' normalize
#' @description normalize GeoMxSet using different normalization methods
#' @param object name of the object class to perform normalization on
#' @param norm_method the normalization method to be applied on the object
#' @param fromElt name of the assayDataElement to normalize
#' @param toElt name of the assayDataElement to store normalized values
#' @param housekeepers optional vector of housekeeper target names
#' @param ... optional arguments
#' @return a NanoStringGeoMxSet object with normalized counts and normalized factors
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'   package = "GeomxTools"
#' )
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' norm_object <- normalize(demoData[1:1000,1:10])
#' @export

setMethod(
  "normalize", "NanoStringGeoMxSet",
  function(object, norm_method = c("quant", "neg", "hk", "subtractBackground"),
           fromElt = "exprs", toElt = "exprs_norm",
           housekeepers = HOUSEKEEPERS, ...) {
    norm_method <- match.arg(norm_method)
    switch(norm_method,
      "quant" = {
        quantileNorm(object,
          toElt = toElt, fromElt = fromElt, ...
        )
      },
      "neg" = {
        negNorm(object,
          toElt = toElt, fromElt = fromElt, ...
        )
      },
      "hk" = {
        hkNorm(object,
          toElt = toElt, fromElt = fromElt,
          housekeepers = housekeepers, ...
        )
      },
      "subtractBackground" = {
        subtractBackground(object,
          toElt = toElt,
          fromElt = fromElt, ...
        )
      }
    )
  }
)

quantileNorm <- function(object, desiredQuantile = .75, toElt, fromElt) {
  ## Get quantile of counts for each sample
  qs <- apply(exprs(object), 2, function(x) stats::quantile(x, desiredQuantile))
  ## Save the normfactors for desired quantile
  if (toElt != "exprs_norm") {
    pData(object)[[paste(toElt, "qFactors", sep = "_")]] <- qs / ngeoMean(qs)
  } else {
    pData(object)[["normFactors"]] <- qs / ngeoMean(qs)
  }
  assayDataElement(object, toElt, validate = TRUE) <- sweep(assayDataElement(object, fromElt), 2L, qs / ngeoMean(qs), FUN = "/")
  return(object)
}

negNorm <- function(object, toElt, fromElt) {
  if (!featureType(object) == "Target") {
    stop("Error: Negative Background normalization is for collapsed data set.
        Run function aggregateCounts() to collapse the probes to targets.\n")
  }
  if (is.null(fData(object)[["Module"]])) {
    stop("Error: Module is not specified in the object. Check your GeoMxSet object. \n")
  }

  if(analyte(object) == "Protein"){
    neg.names <- iggNames(object)
    
    # estimate background:
    negfactor <- apply(exprs(object[neg.names,, drop = FALSE]), 2, function(x){pmax(ngeoMean(x), 1)})
    
    pool_neg_factors <-
      negfactor/exp(mean(log(negfactor)))
    pool_counts <- as.matrix(exprs(object)) %*%
      diag(1 / pool_neg_factors)
    
    # match format of RNA pool_neg_factors
    pool_neg_factors <- t(as.matrix(pool_neg_factors))
    rownames(pool_neg_factors) <- "IgGs"
    
    pool_neg_norm <- list()
    pool_neg_norm[[1L]] <- list(normFactors = pool_neg_factors, norm_exprs = pool_counts)
  }else{
    # check if single panel
    pools <- as.list(unique(fData(object)[["Module"]]))
    pool_neg_norm <- lapply(
      pools,
      function(pool) {
        # Get pool and corresponding target counts
        pool_neg <- fData(object)[which(fData(object)$CodeClass == "Negative" &
                                          fData(object)$Module == pool), "TargetName"]
        if (length(pool_neg) < 1) {
          stop(paste0(
            "Error: No negative could be located for probe pool ",
            pool, "."
          ))
        }
        if (length(pool_neg) > 1) {
          stop(paste0(
            "Error: More than one negative was located for probe pool ",
            pool, "."
          ))
        }
        pool_targets <- fData(object)[which(fData(object)$Module == pool), "TargetName"]
        # Calculate normalization factor and normalized counts
        pool_neg_factors <-
          exprs(object[pool_neg,])/exp(mean(log(exprs(object[pool_neg,]))))
        pool_counts <- as.matrix(exprs(object[pool_targets,])) %*%
          diag(1 / pool_neg_factors[1:ncol(pool_neg_factors)])
        
        norm_list <- list(normFactors = pool_neg_factors, norm_exprs = pool_counts)
        return(norm_list)
      }
    )
  }
  
  

  ## Save the normfactors in desired pData element
  if (toElt != "exprs_norm") {
      pData(object)[[paste(toElt, "negFactors", sep = "_")]] <- t(do.call(rbind, lapply(pool_neg_norm, "[[", 1)))

  } else {
      pData(object)[["normFactors"]] <- t(do.call(rbind, lapply(pool_neg_norm, "[[", 1)))
  }

  # Collapse data back into one data frame
  neg_norm_df <- data.frame(do.call(rbind, lapply(pool_neg_norm, "[[", 2)))
  colnames(neg_norm_df) <- colnames(exprs(object))
  neg_norm_df <- neg_norm_df[rownames(exprs(object)), ]

  ## Save the exprs_norm in desired pData element
  assayDataElement(object, toElt) <- as.matrix(neg_norm_df)
  return(object)
}

hkNorm <- function(object, toElt, fromElt, housekeepers) {
  if (!featureType(object) == "Target") {
    stop("Housekeeping normalization is for collapsed data set.
            Run function aggregateCounts() to collapse the probes to targets.\n")
  } else {
    if(analyte(object) == "Protein" & all(housekeepers == HOUSEKEEPERS)){
      housekeepers <- hkNames(object)
    }
    hksubset <- subset(object, subset = TargetName %in% housekeepers)
    hks <- apply(exprs(hksubset), 2, function(x) ngeoMean(x))
    ## Save the normfactors in desired pData element
    if (toElt != "exprs_norm") {
      pData(object)[[paste(toElt, "hkFactors", sep = "_")]] <- hks / ngeoMean(hks)
    } else {
      pData(object)[["hknormFactors"]] <- hks / ngeoMean(hks)
    }
    assayDataElement(object, toElt) <- sweep(assayDataElement(object, fromElt), 2L, hks / ngeoMean(hks), FUN = "/")
    return(object)
  }
}


# subtract background
subtractBackground <- function(object, toElt, fromElt, byPanel=TRUE) {
    if (featureType(object) == "Target") {
        if(!any(fData(object)$CodeClass == "Negative")){
          stop("Error: No negative could be located for probe pool(s)")
        }
        if(analyte(object) == "Protein"){
          byPanel <- FALSE
        }
        negSet <- negativeControlSubset(object)
        if (byPanel) {
            correctedByPanel <- 
                lapply(unique(fData(object)[["Module"]]), function(currPanel) {
                    panelSet <- subset(object, subset=Module == currPanel)
                    panelNegSet <- subset(negSet, subset=Module == currPanel)
                    panelNegGeo <- 
                        apply(assayDataElement(panelNegSet, elt=fromElt), 
                              2L, ngeoMean)
                    panelCorrectExprs <- 
                        t(assayDataApply(panelSet, MARGIN=1L, 
                          FUN=`-`, t(panelNegGeo), elt = fromElt))
                    colnames(panelCorrectExprs) <- sampleNames(object)
                    return(panelCorrectExprs)
                })
            correctedMatrix <- do.call(rbind, correctedByPanel)
            assayDataElement(object, elt=toElt) <- 
                correctedMatrix[featureNames(object), sampleNames(object)]
        } else {
            allNegGeo <- 
                apply(assayDataElement(negSet, elt=fromElt), 2L, ngeoMean)
            assayDataElement(object, toElt) <-
                t(assayDataApply(object, MARGIN = 1L, 
                  FUN = `-`, t(allNegGeo), elt = fromElt))
        }
        assayDataElement(object, elt=toElt)[
            assayDataElement(object, elt=toElt) < 0] <- 0
    } else {
        warning("Data on probe-level; no background subtraction performed. ",
                "Aggregate counts prior to background subtraction.\n")
    }
    return(object)
}


#' Generate normalization factors 
#' 
#' @description 
#' For use with protein data ONLY.
#' 
#' Generate normalization factors for protein data to determine the best normalization method
#' 
#' @param object name of the object class to subset
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param igg.names names of IgGs, if NULL IgGs will be detected automatically
#' @param hk.names names of HK, if NULL HK will be detected automatically
#' @param area name of area column in annotation sheet, optional
#' @param nuclei name of nuclei column in annotation sheet, optional
#' 
#' @examples
#' proteinData <- readRDS(file= system.file("extdata","DSP_Proteogenomics_Example_Data", 
#' "proteinData.rds", package = "GeomxTools"))
#' 
#' normfactors <- computeNormalizationFactors(object = proteinData)
#' 
#' normfactors_withAreaNuclei <- computeNormalizationFactors(object = proteinData,
#' area = "AOI.Size.um2", nuclei = "Nuclei.Counts")
#' 
#' @export

computeNormalizationFactors <- function(object, igg.names = NULL, hk.names = NULL,
                                          area = NULL, nuclei = NULL) {
  
  if(analyte(object) != "Protein"){
    stop("This function is only for protein data.")
  }
  
  if(is.null(igg.names)){
    igg.names <- iggNames(object)
  }
  
  if(is.null(hk.names)){
    hk.names <- hkNames(object)
  }
  
  segmentAnnotations = sData(object)
  targetAnnotations = fData(object)
  dataset = exprs(object)
  
  # igg and hk factors:
  if (length(igg.names) > 1) {
    igg.factor <- apply(dataset[igg.names, , drop = FALSE], 2, function(x){pmax(ngeoMean(x), 1)})
  }
  if (length(hk.names) > 1) {
    hk.factor <- apply(dataset[hk.names, , drop = FALSE], 2, function(x){pmax(ngeoMean(x), 1)})
  }
  
  # area and nuclei factors:
  if (any(colnames(segmentAnnotations) == area)) {
    area.factor <- as.numeric(segmentAnnotations[[area]])
  }
  if (any(colnames(segmentAnnotations) == nuclei)) {
    nuclei.factor <- as.numeric(segmentAnnotations[[nuclei]])
  }
  
  # matrix of all available factors:
  factornames <- c("igg.factor", "hk.factor", "area.factor", "nuclei.factor")
  names(factornames) <- c("Neg geomean", "HK geomean", "Area", "Nuclei")
  factornames <- factornames[is.element(factornames, ls())]
  
  factors <- c()
  for (fname in factornames) {
    factors <- cbind(factors, get(fname))
  }
  colnames(factors) <- names(factornames)
  
  return(factors)
}

############  NOT USED OR TESTED IN DEV  ############

#' Check QC Flags in the GeoMxSet and removes the probe or sample from the object
#' @rdname checkQCFlags
#' @param object name of the NanoStringGeoMxSet object to check the QC Flags
#' @param ...  for other arguments
#' @return a NanoStringGeoMxSet object probes and samples failing QC removed
#' @export
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'   package = "GeomxTools"
#' )
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' QCobject <- checkQCFlags(demoData)
setGeneric("checkQCFlags",
  signature = c("object"),
  function(object, ...) {
    standardGeneric("checkQCFlags")
  }
)

############  NOT USED OR TESTED IN DEV  ############

#' checkQCFlags
#' @param object name of the NanoStringGeoMxSet object to check the QC Flags
#' @param removeLowLocalOutliers logical, if TRUE it sets outlier counts to zero,  default is FALSE,
#' @param ... optional arguments
#' @return NanoStringGeoMxSet
#' @export
#'
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'   package = "GeomxTools"
#' )
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' QCobject <- checkQCFlags(demoData)
setMethod(
  "checkQCFlags", "NanoStringGeoMxSet",
  function(object, removeLowLocalOutliers = FALSE, ...) {
    ## Remove all samples that failed AOI QC should it have been run
    AOIQCFlags <- protocolData(object)[["QCFlags"]]
    if (is.null(AOIQCFlags)) {
      warning("AOI QC has not been run on this data set.  Proceed with caution.\n")
    } else {
      QCResultsIndex <- which(apply(AOIQCFlags, 1L, function(x) sum(x) == 0L))
      object <- object[, QCResultsIndex]
    }
    ## Remove all low probe count and ratio probes that failed QC
    ProbeQCFlags <- fData(object)[["QCFlags"]]
    if (is.null(ProbeQCFlags)) {
      warning("Probe QC has not been run on this data set.  Proceed with caution.\n")
    } else {
      ProbeQCFlags <- ProbeQCFlags[, c("LowProbeCount", "LowProbeRatio", "GlobalOutlier")]
      probeQCResultsIndex <- which(apply(ProbeQCFlags, 1L, function(x) sum(x) == 0L))
      object <- object[probeQCResultsIndex, ]
    }
    ## Check if option to remove local outliers is set to TRUE
    if (removeLowLocalOutliers == TRUE) {
      ProbeQCFlags <- fData(object)[["QCFlags"]]
      ProbeQCFlags <- ProbeQCFlags[, grepl("HighLocalOutlier|LowLocalOutlier", names(ProbeQCFlags))]
      ## RV: This will remove all probes that has flags. Need to modify this to replace only the sample.
      probeQCResultsIndex <- which(apply(ProbeQCFlags, 1L, function(x) x == TRUE))
      object <- object[probeQCResultsIndex, ]
    }
    return(object)
  }
)
