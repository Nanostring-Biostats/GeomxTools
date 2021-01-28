DEFAULTS <- list(minSaturation=0.7, minReads=10000, minProbeRatio=0.1, 
    minimumCount=10, localOutlierAlpha=0.01, globalOutlierRatio=0.2, 
    loqCutoff=1.0, highCountCutoff=10000)

#NEO to fix the rox comments with new stuff
#' Add QC flags to feature or protocol data
#' 
#' @param object name of the object class to perform QC on
#' \enumerate{
#'     \item{NanoStringGeomxSet, use the NanoStringGeomxSet class}
#' }
#' @param dataDim the dimension of the object to QC on
#' \enumerate{
#'     \item{sample, QC data on the AOI level}
#'     \item{feature, QC data on probe or target level}
#' }
#' @param qcCutoffs list of cutoffs and thresholds to use for QC
#' 
#' @return the object that QC was performed on
#' 
#' @examples
#' 
#NEO set default cutoffs as in DA and make sure set to NULL and check for NULL if not a normal DA cutoff
setMethod("setQCFlags",
    signature(object="NanoStringGeomxSet"),
    function(object, qcCutoffs=DEFAULTS, featType="Probe", ...) {
        qcCutoffs <- checkCutoffs(qcCutoffs)
        #NEO to add featureType accessor to the class plus validation
        object <- setAOIFlags(object=object, qcCutoffs=qcCutoffs)
        if (featType == "Probe") {
            object <- setProbeFlags(object=object, qcCutoffs=qcCutoffs)
        } else if (featType == "Target") {
        #    object <- setTargetFlags(object=object, qcCutoffs=qcCutoffs)
        }
        return(object)
})

#NEO these should match the default calls in DA
setAOIFlags <- function(object, qcCutoffs=DEFAULTS) {
    object <- setSaturationFlags(object=object, 
        cutoff=qcCutoffs[["minSaturation"]])
    object <- setLowReadFlags(object=object, 
        cutoff=qcCutoffs[["minReads"]])
    return(object)
}

setProbeFlags <- function(object, qcCutoffs=DEFAULTS) {
    object <- setProbeRatioFlags(object=object, 
        cutoff=qcCutoffs[["minProbeRatio"]])
    #object <- setProbeCountFlags(object=object, 
    #    cutoff=qcCutoffs[["minimumCount"]])
    #object <- setLocalFlags(object=object, 
    #    cutoff=qcCutoffs[["localOutlierAlpha"]])
    #object <- setGlobalFlags(object=object, 
    #    cutoff=qcCutoffs[["globalOutlierRatio"]])
    return(object)
}

setTargetFlags <- function(object, qcCutoffs=DEFAULTS) {
    object <- 
        setLOQFlags(object=object, cutoff=qcCutoffs[["loqCutoff"]])
    object <- 
        setHighCountFlags(object=object, cutoff=qcCutoffs[["loqCutoff"]])
    return(object)
}

#NEO these are exported so advanced users can use filters not part of default DA pipeline
setSaturationFlags <- function(object, cutoff=DEFAULTS[["minSaturation"]]) {
    percentUnique <- 
        sData(object)["DeduplicatedReads"] / sData(object)["Aligned"]
    percentUnique <- percentUnique > cutoff
    colnames(percentUnique) <- "Saturation"
    object<- appendSampleFlags(object, percentUnique)
    return(object)
}

setLowReadFlags <- function(object, cutoff=DEFAULTS[["minReads"]]) {
    lowReads <- sData(object)["Raw"] < cutoff
    colnames(lowReads) <- "LowReads"
    object <- appendSampleFlags(object, lowReads)
    return(object)
}

# NEO gene needs to be replaced with Target throughout
setProbeRatioFlags <- 
    function(object=object, cutoff=qcCutoffs[["minProbeRatio"]]) {
        rawTargetCounts <- collapseCounts(object)
        rawTargetCounts[["Mean"]] <- 
            apply(rawTargetCounts[, sampleNames(object)], 
                MARGIN=1, FUN=ngeoMean)
        rownames(rawTargetCounts) <- rawTargetCounts[["Target"]]
        targetMeans <- rawTargetCounts[fData(object)[["Target"]], "Mean"]
        probeMeans <- apply(assayDataElement(object, elt="exprs"), 
            MARGIN=1, FUN=ngeoMean)
        probeRatioFlags <- (probeMeans / targetMeans) < cutoff
        probeRatioFlags <- data.frame("RatioFlags"=probeRatioFlags)
        object <- appendFeatureFlags(object, probeRatioFlags)
        return(object)
    }

checkCutoffs <- function(qcCutoffs) {
    if (!all(names(DEFAULTS) %in% names(qcCutoffs))) {
        qcCutoffs <- append(qcCutoffs, 
            DEFAULTS[!(names(DEFAULTS) %in% names(qcCutoffs))])
    }
    return(qcCutoffs)
}

appendSampleFlags <- function(object, currFlags) {
    if("QCFlags" %in% varLabels(protocolData(object))) {
        protocolData(object)[["QCFlags"]] <- 
            cbind(protocolData(object)[["QCFlags"]], currFlags) 
    } else {
        protocolData(object)[["QCFlags"]] <- currFlags
    }
    return(object)
}

appendFeatureFlags <- function(object, currFlags) {
    if("QCFlags" %in% varLabels(featureData(object))) {
        featureData(object)[["QCFlags"]] <- 
            cbind(featureData(object)[["QCFlags"]], currFlags) 
    } else {
        featureData(object)[["QCFlags"]] <- currFlags
    }
    return(object)
}