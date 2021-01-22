DEFAULTS <- list(minPercentUnique=0.7, minReads=10000, minSaturation=50, 
    minProbeRatio=0.1, minimumCount=10, localOutlierAlpha=0.01, 
    globalOutlierRatio=0.2, loqCutoff=1.0, highCountCutoff=10000)

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

setGeneric(".setGeomxFlags", signature=c("dataDim", "featType"), 
    function(object, qcCutoffs, ...) standardGeneric("setGeomxFlags"))

setMethod(".setGeomxFlags", signature(dataDim="sample", featType="Probe"),
    function(object, qcCutoffs=DEFAULTS) {
        .setUniqueFlags(object=object, cutoff=qcCutoffs[["minUnique"]])
        .setLowReadFlags(object=object, cutoff=qcCutoffs[["minReads"]])
})

setMethod(".setGeomxFlags",
    signature(dataDim="sample", featType="Target"),
    function(object, qcCutoffs=DEFAULTS) {
        .setSaturationFlags(object=object, cutoff=qcCutoffs[["minSaturation"]])
})

setMethod(".setGeomxFlags",
    signature(dataDim="feature", featType="Probe"),
    function(object, qcCutoffs=DEFAULTS) {
        .setProbeRatioFlags(object=object, cutoff=qcCutoffs[["minProbeRatio"]])
        .setProbeCountFlags(object=object, cutoff=qcCutoffs[["minimumCount"]])
        .setLocalFlags(object=object, cutoff=qcCutoffs[["localOutlierAlpha"]])
        .setGlobalFlags(object=object,
            cutoff=qcCutoffs[["globalOutlierRatio"]])
})

setMethod(".setGeomxFlags",
    signature(dataDim="feature", featType="Target"),
    function(object, qcCutoffs=DEFAULTS) {
        .setLOQFlags(object=object, cutoff=qcCutoffs[["loqCutoff"]])
        .setHighCountFlags(object=object, cutoff=qcCutoffs[["loqCutoff"]])
})

.setLowReadFlags <- function(object, cutoff=DEFAULTS[["minPercentUnique"]]) {
    dim(demoData)[["Features"]] > 
}

.setLowReadFlags <- function(object, cutoff=DEFAULTS[["minReads"]]) {
    sData(object)[[""]]
}

.setSaturationFlags <- function(object, minSaturation=50) {

}



