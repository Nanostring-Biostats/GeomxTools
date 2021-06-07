DEFAULTS <- list(minSegmentReads=1000, percentTrimmed=80, percentStitched=80, 
    percentAligned=80, percentSaturation=50, minNegativeCount=10, 
    maxNTCCount=60, minNuclei=16000, minArea=20, minProbeCount=10, 
    minProbeRatio=0.1, outlierTestAlpha=0.01, percentFailGrubbs=20, 
    loqCutoff=1.0, highCountCutoff=10000)

#NEO to fix the rox comments with new stuff
#' Add QC flags to feature or protocol data
#' 
#' @param object name of the object class to perform QC on
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
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
    signature(object="NanoStringGeoMxSet"),
    function(object, qcCutoffs=DEFAULTS, ...) {
        qcCutoffs <- checkCutoffs(qcCutoffs)
        object <- setSeqQCFlags(object=object, qcCutoffs=qcCutoffs)
        object <- setBackgroundQCFlags(object=object, qcCutoffs=qcCutoffs)
        object <- setGeoMxQCFlags(object=object, qcCutoffs=qcCutoffs)
        if (featureType(object) == "Probe") {
            object <- setBioProbeQCFlags(object=object, qcCutoffs=qcCutoffs)
        } else if (featureType(object) == "Target") {
        #   object <- setTargetFlags(object=object, qcCutoffs=qcCutoffs)
        } else {
            valid(Object(x))
        }
        return(object)
})

#' Add sequencing QC flags to NanoStringGeoMxSet object protocol data
#' 
#' @param object name of the NanoStringGeoMxSet object to perform QC on
#' @param qcCutoffs a list of qc cutoffs to use
#' \enumerate{
#'     \item{minSegmentReads, 
#'           numeric to flag segments with less than this number of reads}
#'     \item{percentAligned, 
#'           numeric to flag segments with less than this percent of aligned reads}
#'     \item{percentSaturation, 
#'           numeric to flag segments with less than this percent of 
#'           sequencing saturation}
#' }
#' @return \code{NanoStringGeoMxSet} object with \code{QCFlags} data frame 
#'             appended to \code{protocolData}
#' 
#' @examples
#' setSeqQCFlags(demoData, 
#'                  qcCutoffs=list(minSegmentReads=1000, 
#'                                 percentAligned=80, 
#'                                 percentSaturation=50))
#' 
#' @export
#' 
setSeqQCFlags <- function(object, qcCutoffs=DEFAULTS) {
    qcCutoffs <- checkCutoffs(qcCutoffs)
    object <- setLowReadFlags(object=object, 
        cutoff=qcCutoffs[["minSegmentReads"]])
    object <- setTrimmedFlags(object=object, 
        cutoff=qcCutoffs[["percentTrimmed"]])
    object <- setStitchedFlags(object=object, 
        cutoff=qcCutoffs[["percentStitched"]])
    object <- setAlignedFlags(object=object, 
        cutoff=qcCutoffs[["percentAligned"]])
    object <- setSaturationFlags(object=object, 
        cutoff=qcCutoffs[["percentSaturation"]])
    return(object)
}

setLowReadFlags <- function(object, cutoff=DEFAULTS[["minSegmentReads"]]) {
    lowReads <- sData(object)["Raw"] < cutoff
    colnames(lowReads) <- "LowReads"
    object <- appendSampleFlags(object, lowReads)
    return(object)
}

setTrimmedFlags <- function(object, cutoff=DEFAULTS[["percentTrimmed"]]) {
    percentTrim <- 100 * (sData(object)["Trimmed"] / sData(object)["Raw"])
    percentTrim <- percentTrim < cutoff
    colnames(percentTrim) <- "LowTrimmed"
    object<- appendSampleFlags(object, percentTrim)
    return(object)
}

setStitchedFlags <- function(object, cutoff=DEFAULTS[["percentStitched"]]) {
    percentStitch <- 100 * (sData(object)["Stitched"] / sData(object)["Raw"])
    percentStitch <- percentStitch < cutoff
    colnames(percentStitch) <- "LowStitched"
    object<- appendSampleFlags(object, percentStitch)
    return(object)
}

setAlignedFlags <- function(object, cutoff=DEFAULTS[["percentAligned"]]) {
    percentAlign <- 100 * (sData(object)["Aligned"] / sData(object)["Raw"])
    percentAlign <- percentAlign < cutoff
    colnames(percentAlign) <- "LowAligned"
    object<- appendSampleFlags(object, percentAlign)
    return(object)
}

setSaturationFlags <- function(object, cutoff=DEFAULTS[["percentSaturation"]]) {
    percentSaturated <- 
        1 - sData(object)["DeduplicatedReads"] / sData(object)["Aligned"]
    percentSaturated <- percentSaturated < cutoff
    colnames(percentSaturated) <- "LowSaturation"
    object<- appendSampleFlags(object, percentSaturated)
    return(object)
}


#' Add background QC flags to NanoStringGeoMxSet object protocol data
#' 
#' @param object name of the NanoStringGeoMxSet object to perform QC on
#' @param qcCutoffs a list of qc cutoffs to use
#' \enumerate{
#'     \item{minNegativeCount, 
#'           numeric to flag segments with less than this number of counts}
#'     \item{maxNTCCount, 
#'           numeric to flag segments with more than this number of NTC counts}
#' }
#' @return \code{NanoStringGeoMxSet} object with \code{QCFlags} data frame 
#'             appended to \code{protocolData}
#' 
#' @examples
#' setBackgroundQCFlags(demoData, 
#'                  qcCutoffs=list(minNegativeCount=10, 
#'                                 maxNTCCount=60))
#' 
#' @export
#' 
setBackgroundQCFlags <- function(object, qcCutoffs=DEFAULTS) {
    qcCutoffs <- checkCutoffs(qcCutoffs)
    object <- setLowNegFlags(object=object, 
        cutoff=qcCutoffs[["minNegativeCount"]])
    object <- setHighNTCFlags(object=object, 
        cutoff=qcCutoffs[["maxNTCCount"]])
    return(object)
}

setLowNegFlags <- function(object, cutoff=DEFAULTS[["minNegativeCount"]]) {
    negativeGeoMeans <- 
        esBy(negativeControlSubset(object), 
             GROUP="Module", 
             FUN=function( x ) { 
                 assayDataApply( x, MARGIN = 2, FUN=ngeoMean, elt="exprs" ) 
             }) 
    lowNegs <- 
        data.frame("lowNegatives"=apply(negativeGeoMeans < cutoff, 1, sum) > 0)
    object <- appendSampleFlags(object, lowNegs)
    return(object)
}

setHighNTCFlags <- function(object, cutoff=DEFAULTS[["maxNTCCount"]]) {
    if (!is.null(sData(object)[["NTC"]])) {
        highNTC <- sData(object)["NTC"] > cutoff
        colnames(highNTC) <- "highNTC"
        object <- appendSampleFlags(object, highNTC)
    }
    return(object)
}

#' Add GeoMx segment QC flags to NanoStringGeoMxSet object protocol data
#' 
#' @param object name of the NanoStringGeoMxSet object to perform QC on
#' @param qcCutoffs a list of qc cutoffs to use
#' \enumerate{
#'     \item{minNuclei, 
#'           numeric to flag segments with less than this number of nuclei}
#'     \item{minArea, 
#'           numeric to flag segments with less than this um^2 area}
#' }
#' @return \code{NanoStringGeoMxSet} object with \code{QCFlags} data frame 
#'         appended to \code{protocolData}
#' 
#' @examples
#' setGeoMxQCFlags(demoData, 
#'                  qcCutoffs=list(minNuclei=16000, 
#'                                 minArea=20))
#' 
#' @export
#' 
setGeoMxQCFlags <- function(object, qcCutoffs=DEFAULTS) {
    qcCutoffs <- checkCutoffs(qcCutoffs)
    object <- setNucleiFlags(object=object, 
        cutoff=qcCutoffs[["minNuclei"]])
    object <- setAreaFlags(object=object, 
        cutoff=qcCutoffs[["minArea"]])
    return(object)
}

setNucleiFlags <- function(object, cutoff=DEFAULTS[["minNuclei"]]) {
    if (!is.null(sData(object)[["nuclei"]])) {
        lowNuclei <- sData(object)["nuclei"] < cutoff
        colnames(lowNuclei) <- "lowNuclei"
        object <- appendSampleFlags(object, lowNuclei)
    }
    return(object)
}

setAreaFlags <- function(object, cutoff=DEFAULTS[["minArea"]]) {
    if (!is.null(sData(object)[["area"]])) {
        lowArea <- sData(object)["area"] < cutoff
        colnames(lowArea) <- "lowArea"
        object <- appendSampleFlags(object, lowArea)
    }
    return(object)
}

#' Add probe QC flags to NanoStringGeoMxSet object feature data
#' 
#' @param object name of the NanoStringGeoMxSet object to perform QC on
#' @param qcCutoffs a list of qc cutoffs to use
#' \enumerate{
#'     \item{minProbeRatio, 
#'           numeric between 0 and 1 to flag probes that have 
#'           (geomean probe in all segments) / (geomean probes within target)
#'           less than or equal to this ratio}
#'     \item{percentFailGrubbs, 
#'           numeric to flag probes that fail Grubb's test in 
#'           greater than or equal this percent of segments}
#' }
#' @return \code{NanoStringGeoMxSet} object with \code{QCFlags} data frame 
#'         appended to \code{protocolData}
#' 
#' @examples
#' setBioProbeQCFlags(demoData, 
#'                    qcCutoffs=list(minProbeRatio=0.1,
#'                                   percentFailGrubbs=20))
#' 
#' @export
#' 
setBioProbeQCFlags <- function(object, qcCutoffs=DEFAULTS) {
    qcCutoffs <- checkCutoffs(qcCutoffs)
    object <- setProbeRatioFlags(object=object, 
        cutoff=qcCutoffs[["minProbeRatio"]])
    object <- setGrubbsFlags(object=object,
                             minCount=qcCutoffs[["minProbeCount"]],
                             cutoff=qcCutoffs[["outlierTestAlpha"]],
                             percFail=qcCutoffs[["percentFailGrubbs"]])
    return(object)
}

setProbeRatioFlags <- function(object, 
                               cutoff=DEFAULTS[["minProbeRatio"]]) {
    rawTargetCounts <- collapseCounts(object)
    rawTargetCounts[["Mean"]] <- 
        apply(rawTargetCounts[, sampleNames(object)], 
            MARGIN=1, FUN=ngeoMean)
    rownames(rawTargetCounts) <- rawTargetCounts[["TargetName"]]
    targetMeans <- rawTargetCounts[fData(object)[["TargetName"]], "Mean"]
    probeMeans <- apply(assayDataElement(object, elt="exprs"), 
        MARGIN=1, FUN=ngeoMean)
    probeRatioFlags <- (probeMeans / targetMeans) < cutoff
    probeRatioFlags <- data.frame("LowProbeRatio"=probeRatioFlags)
    object <- appendFeatureFlags(object, probeRatioFlags)
    return(object)
}

setGrubbsFlags <- function(object, 
                           cutoff=DEFAULTS[["outlierTestAlpha"]], 
                           minCount=DEFAULTS[["minProbeCount"]],
                           percFail=DEFAULTS[["percentFailGrubbs"]]) {
    # Skip targets with less than 3 probes
    multiProbeTable <- with(object, table(TargetName, Module)) >= 3L
    indices <- which(multiProbeTable, arr.ind=TRUE)
    if (length(indices) != dim(multiProbeTable)[1L]) {
        targs <- rownames(multiProbeTable)[as.vector(indices[, "TargetName"])]
        mods <- colnames(multiProbeTable)[as.vector(indices[, "Module"])]
        multiProbeList <- 
            unlist(lapply(seq_along(targs), function(x){
                featureNames(subset(object, 
                    subset=TargetName == targs[x] & Module == mods[x]))
            }))
        multiObject <- object[multiProbeList, ]
    } else {
        multiObject <- object
    }

    # Skip targets with all counts less than minimum probe count cutoff
    multiTab <- unique(fData(multiObject)[c("TargetName", "Module")])
    rownames(multiTab) <- NULL
    multiTab[["countPass"]] <- 
        apply(multiTab, 
              MARGIN=1L, 
              FUN=function(x){
                  max(exprs(subset(multiObject, 
                      subset=TargetName == x[["TargetName"]] & 
                             Module == x[["Module"]]))) > minCount
              })
    if (sum(multiTab[["countPass"]]) > 0) {
        multiTabPass <- multiTab[multiTab[["countPass"]], ]
    } else {
        # Skip grubb's test altogether if counts are low for remaining targets
        grubbsFlags <- data.frame(GrubbsOutlier=rep(FALSE, dim(object)[1L]),
                                  row.names=featureNames(object))
        object <- appendFeatureFlags(object, grubbsFlags)
        return(object)
    }
    if (dim(multiTabPass)[1L] != dim(multiTab)[1L]) {
        passProbeList <- 
            unlist(lapply(seq_along(rownames(multiTabPass)), function(x){
                featureNames(subset(object, 
                    subset=TargetName == multiTabPass[["TargetName"]][x] & 
                               Module == multiTabPass[["Module"]][x]))
            }))
        minPassObject <- multiObject[passProbeList, ]
    } else {
        minPassObject <- multiObject
    }

    #Perform Grubb's outlier test
    probeCounts <- 
        setDT(cbind(fData(minPassObject)[, c("RTS_ID", "TargetName", "Module")], 
            assayDataElement(minPassObject, elt="exprs")))
    probeCounts <- melt(probeCounts, 
        id.vars=c("RTS_ID", "TargetName", "Module"), 
        variable.name="Sample_ID", 
        value.name="Count", variable.factor=FALSE)
    probeCounts[, Count:=logt(Count)]
    probeCounts[, "LowLocalOutlier"] <- FALSE
    probeCounts[, "HighLocalOutlier"] <- FALSE
    probeCounts <- probeCounts[, suppressWarnings(grubbsFlag(.SD, alpha=cutoff)), 
        by=.(TargetName, Module, Sample_ID)]
    lowFlags <- 
        as.data.frame(dcast(probeCounts, 
                            RTS_ID ~ Sample_ID, value.var="LowLocalOutlier"), stringsAsFactor=FALSE)
    highFlags <- 
        as.data.frame(dcast(probeCounts, 
                            RTS_ID ~ Sample_ID, value.var="HighLocalOutlier"), stringsAsFactor=FALSE)
    rownames(lowFlags) <- lowFlags[["RTS_ID"]]
    rownames(highFlags) <- highFlags[["RTS_ID"]]
    lowFlags <- lowFlags[, colnames(lowFlags) != "RTS_ID"]
    highFlags <- highFlags[, colnames(highFlags) != "RTS_ID"]
    lowFlags[["lowOutliers"]] <- 
        (apply(lowFlags, MARGIN=1L, FUN=mean) * 100) >= percFail
    highFlags[["highOutliers"]] <- 
        (apply(highFlags, MARGIN=1L, FUN=mean) * 100) >= percFail
    grubbsFlags <- data.frame(GrubbsOutlier=rep(FALSE, dim(object)[1L]),
                              row.names=featureNames(object))
    if (sum(lowFlags[["lowOutliers"]]) > 0) {
        grubbsFlags[rownames(lowFlags[lowFlags[["lowOutliers"]],]), 
                    "GrubbsOutlier"] <- TRUE
    }
    if (sum(highFlags[["highOutliers"]]) > 0) {
        grubbsFlags[rownames(highFlags[highFlags[["highOutliers"]],]), 
                    "GrubbsOutlier"] <- TRUE
    }
    object <- appendFeatureFlags(object, grubbsFlags)
    return(object)
}

checkCutoffs <- function(qcCutoffs) {
    if (suppressWarnings(!all(lapply(qcCutoffs, is.numeric)))) {
        stop("qcCutoffs must be numeric values")
    }
    if (!all(names(DEFAULTS) %in% names(qcCutoffs))) {
        qcCutoffs <- append(qcCutoffs, 
            DEFAULTS[!(names(DEFAULTS) %in% names(qcCutoffs))])
    }
    return(qcCutoffs)
}

appendSampleFlags <- function(object, currFlags) {
    if ("QCFlags" %in% varLabels(protocolData(object))) {
        protocolData(object)[["QCFlags"]] <- 
            cbind(protocolData(object)[["QCFlags"]], currFlags)
    } else {
        protocolData(object)[["QCFlags"]] <- currFlags
    }
    return(object)
}

appendFeatureFlags <- function(object, currFlags) {
    if ("QCFlags" %in% varLabels(featureData(object))) {
        featureData(object)[["QCFlags"]] <- 
            cbind(featureData(object)[["QCFlags"]], currFlags) 
    } else {
        featureData(object)[["QCFlags"]] <- currFlags
    }
    return(object)
}

# grubbs.flag
# helper function to remove outliers using Grubbs' test given the controlled type I error alpha
# INPUT
#   x = named vector
#   alpha = indicator of the expected Type I erorr frequency
#   logt = boolean to log 10 transform data prior to outlier test
#   min_count = integer of minimum expected counts for testing
# OUTPUT
#   vector of TRUE / FALSE For flagged or not flagged
grubbsFlag <- function(countDT, alpha=0.01) {
    if (dim(countDT)[1] < 3 | all(countDT[, Count] == countDT[1][, Count])) {
        return(countDT)
    }
    grubbsResult <- outliers::grubbs.test(countDT[, Count], two.sided=TRUE)
    if (grubbsResult$p.value < alpha) {
        if (grepl("lowest", tolower(grubbsResult$alternative))) {
            countDT[Count == min(Count), "LowLocalOutlier"] <- TRUE
        } else {
            countDT[Count == max(Count), "HighLocalOutlier"] <- TRUE
        }
    } else {
        return(countDT)
    }
    return(countDT)
}


removeFlagProbes <- function(object, removeFlagCols=NULL) {
    
}

#NEO make sure that when collapsed count occurs feature data QCFlags is removed
setTargetFlags <- function(object, qcCutoffs=DEFAULTS) {
    object <- 
        setLOQFlags(object, cutoff=qcCutoffs[["loqCutoff"]])
    object <- 
        setHighCountFlags(object=object, cutoff=qcCutoffs[["highCountCutoff"]])
    return(object)
}

setLOQFlags <- function(object, cutoff=DEFAULTS[["loqCutoff"]]) {
        #NEO need negative values for LOQ
        if (featureType(object) == "Target") {
            negativeObject <- negativeControlSubset(object)
            LOQs <- esApply(negativeObject, MARGIN=2, FUN=function(x) {
                ngeoMean(x) * ngeoSD(x) ^ cutoff})
            LOQFlags <- esApply(object, MARGIN=1, FUN=function(x) {
                x > LOQs})
            object <- appendFeatureFlags(object, LOQFlags)
        } else {
            warning(paste("Incorrect feature type.",
                "Feature type should be Target."))
        }
        return(object)
}

setProbeCountFlags <- function(object, cutoff=DEFAULTS[["minProbeCount"]]) {
        probeCountFlags <- apply(assayDataElement(object, elt="exprs"), 
            MARGIN=1, FUN=function(x, minCount){
                all(x < minCount)
            }, minCount=cutoff)
        probeCountFlags <- data.frame("LowProbeCount"=probeCountFlags)
        object <- appendFeatureFlags(object, probeCountFlags)
        return(object)
}
