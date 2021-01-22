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
    function(object, dataDim="sample", qcCutoffs=DEFAULTS, ...) {
        qcCutoffs <- append(qcCutoffs, 
            DEFAULTS[!(names(defaultCutoffs) %in% names(qcCutoffs))])
        #NEO to add featureType accessor to the class plus validation
        .setGeomxFlags(object=object, dataDim=dataDim, 
            featType=featureType(object), qcCutoffs=qcCutoffs)
})

#NEO these (setGeomxFlags) should match the default calls in DA
setGeneric(".setGeomxFlags", 
    function(object, qcCutoffs, ...) standardGeneric("setGeomxFlags"))

setMethod(".setGeomxFlags", 
    signature(dataDim="sample"),
    function(object, qcCutoffs=DEFAULTS) {
        object <- setSaturationFlags(object=object, 
            cutoff=qcCutoffs[["minSaturation"]])
        object <- setLowReadFlags(object=object, 
            cutoff=qcCutoffs[["minReads"]])
})

setMethod(".setGeomxFlags",
    signature(dataDim="feature", featType="Probe"),
    function(object, qcCutoffs=DEFAULTS) {
        object <- setProbeRatioFlags(object=object, 
            cutoff=qcCutoffs[["minProbeRatio"]])
        object <- setProbeCountFlags(object=object, 
            cutoff=qcCutoffs[["minimumCount"]])
        object <- setLocalFlags(object=object, 
            cutoff=qcCutoffs[["localOutlierAlpha"]])
        object <- setGlobalFlags(object=object, 
            cutoff=qcCutoffs[["globalOutlierRatio"]])
})

setMethod(".setGeomxFlags",
    signature(dataDim="feature", featType="Target"),
    function(object, qcCutoffs=DEFAULTS) {
        object <- 
            setLOQFlags(object=object, cutoff=qcCutoffs[["loqCutoff"]])
        object <- 
            setHighCountFlags(object=object, cutoff=qcCutoffs[["loqCutoff"]])
        return(object)
})

#NEO these are exported so advanced users can use filters not part of default DA pipeline
setSaturationFlags <- function(object, cutoff=DEFAULTS[["minSaturation"]]) {
    percentUnique <- 
        sData(object)[["DeduplicatedReads"]] / sData(object)[["Aligned"]]
    protocolData(object)[["QCFlags"]][, "Saturation"] <- 
        percentUnique > cutoff
    return(object)
}

setLowReadFlags <- function(object, cutoff=DEFAULTS[["minReads"]]) {
    protocolData(object)[["QCFlags"]][, "LowReads"] <- 
        sData(object)[["Raw"]] < cutoff
    return(object)
}

# NEO gene needs to be replaced with Target throughout
setProbeRatioFlags <- function(object=object, 
    cutoff=qcCutoffs[["minProbeRatio"]]) {
        #NEO make generic aggregate counts function so this can be called many times in QC or to generate new class
        gene_assay <- reshape2::dcast(probe_assay, Gene + Pool ~ Sample_ID, 
            value.var = 'Count', fun.aggregate = ngeoMean, fill = 1)
        targetMeans <- lapply(featureNames(object), assayDataElement(object, elt="exprs"))
}

