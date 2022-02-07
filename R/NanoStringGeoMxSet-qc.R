DEFAULTS <- list(minSegmentReads=1000, percentTrimmed=80, percentStitched=80, 
    percentAligned=80, percentSaturation=50, minNegativeCount=10, 
    maxNTCCount=60, minNuclei=200, minArea=16000, minProbeCount=10, 
    minProbeRatio=0.1, outlierTestAlpha=0.01, percentFailGrubbs=20, 
    loqCutoff=1.0, highCountCutoff=10000)

#' Add segment QC flags to protocol data
#' 
#' @param object name of the object class to perform QC on
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param qcCutoffs list of cutoffs and thresholds to use for QC
#' 
#' @return \code{NanoStringGeoMxSet} object with \code{QCFlags} data frame 
#'             appended to \code{protocolData}
#' 
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'                        package="GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' setSegmentQCFlags(demoData, 
#'                   qcCutoffs=list(minSegmentReads=1000, 
#'                                  percentAligned=80, 
#'                                  percentSaturation=50,
#'                                  minNegativeCount=10, 
#'                                  maxNTCCount=60,
#'                                  minNuclei=16000, 
#'                                  minArea=20))
#' 
#' @export
#' 
setSegmentQCFlags <- function(object, qcCutoffs=DEFAULTS) {
    qcCutoffs <- checkCutoffs(qcCutoffs)
    object <- setSeqQCFlags(object=object, qcCutoffs=qcCutoffs)
    object <- setBackgroundQCFlags(object=object, qcCutoffs=qcCutoffs)
    object <- setGeoMxQCFlags(object=object, qcCutoffs=qcCutoffs)
    return(object)
}

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
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'                        package="GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
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
    colnames(percentTrim) <- "Trimmed (%)"
    protocolData(object)[["Trimmed (%)"]] <- percentTrim
    percentTrim <- percentTrim < cutoff
    colnames(percentTrim) <- "LowTrimmed"
    object<- appendSampleFlags(object, percentTrim)
    return(object)
}

setStitchedFlags <- function(object, cutoff=DEFAULTS[["percentStitched"]]) {
    percentStitch <- 100 * (sData(object)["Stitched"] / sData(object)["Raw"])
    colnames(percentStitch) <- "Stitched (%)"
    protocolData(object)[["Stitched (%)"]] <- percentStitch
    percentStitch <- percentStitch < cutoff
    colnames(percentStitch) <- "LowStitched"
    object<- appendSampleFlags(object, percentStitch)
    return(object)
}

setAlignedFlags <- function(object, cutoff=DEFAULTS[["percentAligned"]]) {
    percentAlign <- 100 * (sData(object)["Aligned"] / sData(object)["Raw"])
    colnames(percentAlign) <- "Aligned (%)"
    protocolData(object)[["Aligned (%)"]] <- percentAlign
    percentAlign <- percentAlign < cutoff
    colnames(percentAlign) <- "LowAligned"
    object<- appendSampleFlags(object, percentAlign)
    return(object)
}

setSaturationFlags <- function(object, cutoff=DEFAULTS[["percentSaturation"]]) {
    percentSaturated <- 
        100 * (1 - sData(object)["DeduplicatedReads"] / sData(object)["Aligned"])
    colnames(percentSaturated) <- "Saturated (%)"
    protocolData(object)[["Saturated (%)"]] <- percentSaturated
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
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'                        package="GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' setBackgroundQCFlags(demoData, 
#'                  qcCutoffs=list(minNegativeCount=10, 
#'                                 maxNTCCount=60))
#' 
#' @export
#' 
setBackgroundQCFlags <- function(object, qcCutoffs=DEFAULTS) {
    qcCutoffs <- checkCutoffs(qcCutoffs)
    if(any(fData(object)$AnalyteType == "RNA")){
      object <- setLowNegFlags(object=object, 
                               cutoff=qcCutoffs[["minNegativeCount"]])
    }
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
        data.frame("LowNegatives"=apply(negativeGeoMeans < cutoff, 1, sum) > 0)
    object <- appendSampleFlags(object, lowNegs)
    return(object)
}

setHighNTCFlags <- function(object, cutoff=DEFAULTS[["maxNTCCount"]]) {
    if (!is.null(sData(object)[["NTC"]])) {
        highNTC <- sData(object)["NTC"] > cutoff
        colnames(highNTC) <- "HighNTC"
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
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'                        package="GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
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
        colnames(lowNuclei) <- "LowNuclei"
        object <- appendSampleFlags(object, lowNuclei)
    }
    return(object)
}

setAreaFlags <- function(object, cutoff=DEFAULTS[["minArea"]]) {
    if (!is.null(sData(object)[["area"]])) {
        lowArea <- sData(object)["area"] < cutoff
        colnames(lowArea) <- "LowArea"
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
#' @param removeLocalOutliers boolean to determine if 
#'     local outliers should be excluded from \code{exprs} matrix
#'
#' @return \code{NanoStringGeoMxSet} object with \code{QCFlags} data frame 
#'         appended to \code{protocolData}
#' 
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'                        package="GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' demoData <- shiftCountsOne(demoData, elt="exprs", useDALogic=TRUE)
#' setBioProbeQCFlags(demoData, 
#'                    qcCutoffs=list(minProbeRatio=0.1,
#'                                   percentFailGrubbs=20),
#'                    removeLocalOutliers=TRUE)
#' 
#' @export
#' 
setBioProbeQCFlags <- function(object, 
                               qcCutoffs=DEFAULTS, 
                               removeLocalOutliers=TRUE) {
    qcCutoffs <- checkCutoffs(qcCutoffs)
    object <- setProbeRatioFlags(object=object, 
        cutoff=qcCutoffs[["minProbeRatio"]])
    object <- setGrubbsFlags(object=object,
                             minCount=qcCutoffs[["minProbeCount"]],
                             alphaCutoff=qcCutoffs[["outlierTestAlpha"]],
                             percFail=qcCutoffs[["percentFailGrubbs"]])
    if (removeLocalOutliers) {
        object <- changeLocalsToNA(object)
    }
    return(object)
}

setProbeRatioFlags <- function(object, 
                               cutoff=DEFAULTS[["minProbeRatio"]]) {
    # Skip targets with single probes
    multiProbeTable <- with(object, table(TargetName)) > 1L
    multiProbeTargs <- 
        names(multiProbeTable)[which(multiProbeTable, arr.ind=TRUE)]
    if (length(multiProbeTargs) > 0) {
        multiObject <- 
            object[fData(object)[["TargetName"]] %in% multiProbeTargs, ]
    } else {
        warning("Object has no multi-probe targets. ",
                "No ratio testing can be performed.")
        return(object)
    }

    probeCounts <- data.table(cbind(fData(multiObject)[, c("TargetName", 
                                                           "Module")],
                                    list(RTS_ID=featureNames(multiObject)),
                                    exprs(multiObject)))
    collapsedCounts <- as.data.frame(probeCounts[, lapply(.SD, ngeoMean),
                                              .SDcols=!c("RTS_ID"),
                                              by=c("TargetName", "Module")])
    targetMeans <- 
        apply(collapsedCounts[, 
                  colnames(collapsedCounts) %in% sampleNames(multiObject)],
                  1, ngeoMean)
    names(targetMeans) <- collapsedCounts[["TargetName"]]

    probeMeans <- apply(assayDataElement(multiObject, elt="exprs"), 
        MARGIN=1, FUN=ngeoMean)
    multiProbeRatios <- probeMeans / targetMeans[probeCounts$TargetName]

    probeRatios <- data.frame("ProbeRatio"=rep(1, dim(object)[1L]), 
                              row.names=featureNames(object))
    probeRatios[names(multiProbeRatios), "ProbeRatio"] <- 
        multiProbeRatios

    featureData(object)[["ProbeRatio"]] <- probeRatios
    probeRatioFlags <- probeRatios <= cutoff
    colnames(probeRatioFlags) <- "LowProbeRatio"
    object <- appendFeatureFlags(object, probeRatioFlags)
    return(object)
}

setGrubbsFlags <- function(object, 
                           alphaCutoff=DEFAULTS[["outlierTestAlpha"]], 
                           minCount=DEFAULTS[["minProbeCount"]],
                           percFail=DEFAULTS[["percentFailGrubbs"]]) {
    # Skip targets with less than 3 probes
    multiProbeTable <- with(object, table(TargetName)) >= 3L
    multiProbeTargs <- 
        names(multiProbeTable)[which(multiProbeTable, arr.ind=TRUE)]
    if (length(multiProbeTargs) > 0) {
        multiObject <- 
            object[fData(object)[["TargetName"]] %in% multiProbeTargs, ]
    } else {
        warning("Object has no targets with at least 3 probes. ",
                "No outlier testing can be performed.")
        return(object)
    }

    if (!"ProbeRatio" %in% fvarLabels(object)) {
        warning("Probe ratio QC has not yet been performed. ",
                "Suggests running probe ratio QC then re-running Grubbs QC.")
    } else {
        multiObject <- 
            multiObject[!fData(multiObject)[["QCFlags"]][, "LowProbeRatio"], ]
    }

    grubbsResults <- 
        lapply(unique(fData(multiObject)[["Module"]]), FUN=function(x) {
            modObject <- subset(multiObject, subset=Module == x)
            assayDataElement(modObject, elt="log10") <- 
                logtBase(exprs(modObject), base=10)
            targetOutliers <- 
                esBy(modObject, GROUP="TargetName", simplify=FALSE, FUN=function(y) {
                    currExprs <- assayDataElement(y, elt="log10")
                    currResult <- apply(currExprs, 2, function(z) { 
                        if (max(z) < logtBase(minCount, base=10)) {
                            return(NA)
                        }
                        if (length(unique(z)) == 1) {
                            return(NA)
                        }
                        grubbsTest <- 
                            suppressWarnings(
                                outliers::grubbs.test(z,  
                                                      two.sided=TRUE, 
                                                      type=10))
                        if(grubbsTest$p.value < alphaCutoff & 
                               !is.null(names(grubbsTest$p.value))) {
                            if (grepl("lowest", tolower(grubbsTest$alternative))) {
                                outlierDirection <- 
                                    data.frame(LowLocalOutlier=TRUE, 
                                               HighLocalOutlier=FALSE, 
                                               RTS_ID=names(grubbsTest$p.value))
   
                            } else {
                                outlierDirection <- 
                                    data.frame(LowLocalOutlier=FALSE, 
                                               HighLocalOutlier=TRUE, 
                                               RTS_ID=names(grubbsTest$p.value))
                            }
                            return(outlierDirection)
                        } else {
                            return(NA)
                        }
                    })
                    names(currResult) <- colnames(currExprs)
                    currResult <- currResult[!is.na(currResult)]
                    appendedResult <-
                        lapply(names(currResult), function(q) {
                            toAppendResult <- currResult[[q]]
                            toAppendResult[["Sample_ID"]] <- q
                            return(toAppendResult)
                        })
                    ifelse(length(appendedResult) < 1, 
                           return(), 
                           return(do.call(rbind, appendedResult)))
                })
            targetOutliers <- 
                targetOutliers[!targetOutliers %in% list(NULL)]
            if(length(targetOutliers) > 0) {
                targetOutliers <- do.call(rbind, targetOutliers)
                return(targetOutliers)
            } else {
                return()
            }
        })
    
    featureData(object)[["OutlierFrequency"]] <- 0
    globalFlags <- data.frame(GlobalGrubbsOutlier=rep(FALSE, dim(object)[1L]),
                              row.names=featureNames(object))
    localGrubbsFlags <- 
        data.frame(matrix(nrow=dim(object)[1L], 
                          ncol=dim(object)[2L], 
                          dimnames=list(featureNames(object), sampleNames(object))),
                          check.names=FALSE)
    
    grubbsResults <- grubbsResults[!grubbsResults %in% list(NULL)]
    if(length(grubbsResults) > 0) {
        grubbsResults <- do.call(rbind, grubbsResults)
        
        outlierFreqs <- table(grubbsResults[["RTS_ID"]]) / dim(object)[2L]
        fData(object)[names(outlierFreqs), "OutlierFrequency"] <- outlierFreqs
        globalFlags[["GlobalGrubbsOutlier"]] <- 
            100 * fData(object)[["OutlierFrequency"]] >= percFail
       
        subsetOfLocals <- 
            lapply(unique(grubbsResults[["RTS_ID"]]), function(currOut) {
                currFlags <- localGrubbsFlags[currOut, , drop=FALSE]
                lowAOIs <- 
                    grubbsResults[(grubbsResults[["RTS_ID"]] == currOut & 
                                       grubbsResults[["LowLocalOutlier"]]), 
                                  "Sample_ID"]
                highAOIs <- 
                     grubbsResults[(grubbsResults[["RTS_ID"]] == currOut & 
                                       grubbsResults[["HighLocalOutlier"]]), 
                                   "Sample_ID"]
                if (length(lowAOIs) > 0) {currFlags[, lowAOIs] <- 
                    rep(TRUE, length(lowAOIs))}
                if (length(highAOIs) > 0) {currFlags[, highAOIs] <- 
                    rep(TRUE, length(highAOIs))}
                return(currFlags)
        })
        subsetOfLocals <- do.call(rbind, subsetOfLocals)
        localGrubbsFlags[rownames(subsetOfLocals), 
                         colnames(subsetOfLocals)] <- subsetOfLocals
        localGrubbsFlags <- 
            localGrubbsFlags[featureNames(object), sampleNames(object)]
    }
    localGrubbsFlags[is.na(localGrubbsFlags)] <- FALSE
    LocalGrubbsOutlier <- 
        data.frame(LocalGrubbsOutlier=localGrubbsFlags, check.names=FALSE)
    object <- appendFeatureFlags(object, globalFlags)
    object <- appendFeatureFlags(object, LocalGrubbsOutlier)
    return(object)
}

changeLocalsToNA <- function(object) {
    localColumns <- grepl("LocalGrubbs", colnames(fData(object)[["QCFlags"]]))
    if (sum(localColumns) > 0L) {
        assayDataElement(object, elt="preLocalRemoval") <- exprs(object)
        exprs(object)[fData(object)[["QCFlags"]][, localColumns] == TRUE] <- NA
    } else {
        warning("Local outlier test not performed yet. ",
                "No change in local outlier values performed.")
    }
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
    currQCName <- colnames(currFlags)
    if ("QCFlags" %in% varLabels(protocolData(object))) {
        if (currQCName %in% colnames(protocolData(object)[["QCFlags"]])) {
            protocolData(object)[["QCFlags"]] <- 
                protocolData(object)[["QCFlags"]][, 
                    !colnames(protocolData(object)[["QCFlags"]]) %in% currQCName]
        }
        protocolData(object)[["QCFlags"]] <- 
            cbind(protocolData(object)[["QCFlags"]], currFlags)
    } else {
        protocolData(object)[["QCFlags"]] <- currFlags
    }
    return(object)
}

appendFeatureFlags <- function(object, currFlags) {
    currQCName <- colnames(currFlags)
    if ("QCFlags" %in% varLabels(featureData(object))) {
        if (any(currQCName %in% colnames(featureData(object)[["QCFlags"]]))) {
            featureData(object)[["QCFlags"]] <- 
                featureData(object)[["QCFlags"]][, 
                    !colnames(featureData(object)[["QCFlags"]]) %in% currQCName]
        }
        featureData(object)[["QCFlags"]] <- 
            cbind(featureData(object)[["QCFlags"]], currFlags) 
    } else {
        featureData(object)[["QCFlags"]] <- currFlags
    }
    return(object)
}

############Not used##########################################################
#' Add QC flags to feature and protocol data simultaneously
#' 
#' @param object name of the object class to perform QC on
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param qcCutoffs list of cutoffs and thresholds to use for QC
#' @param ... optional parameters to pass
#' 
#' @return the object that QC was performed on
#' 
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'                        package="GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' setQCFlags(demoData)
#' 
#' @export
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
        }
        return(object)
})

removeFlagProbes <- function(object, removeFlagCols=NULL) {
    
}

#NEO make sure that when collapsed count occurs feature data QCFlags is removed
setTargetFlags <- function(object, qcCutoffs=DEFAULTS) {
    if(any(fData(object)$AnalyteType == "RNA")){
      object <- 
        setLOQFlags(object, cutoff=qcCutoffs[["loqCutoff"]])
    }  
    object <- 
        setHighCountFlags(object=object, cutoff=qcCutoffs[["highCountCutoff"]])
    return(object)
}

setHighCountFlags <- function(object, cutoff=DEFAULTS[["highCountCutoff"]]) {
    cutoff <- checkCutoffs(cutoff)
    return(object)
}

setLOQFlags <- function(object, cutoff=DEFAULTS[["loqCutoff"]]) {
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


## Protein QC

concordance_plot <- function(mat, col = rgb(0, 0, 0, 0.5),
                             collegend = NULL, legend.main = NULL,
                             main = "", pch = 16, cex = 1.5, ...) {

    # subsidiary function for printing the SD of a ratio, for use by pairs():
    print.sd.log.ratio <- function(x, y, digits = 2, prefix = "", cex.cor = 1, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        pairwise.stat <- sd(log2(pmax(x, 1)) - log2(pmax(y, 1)), na.rm = TRUE)
        txt <- round(pairwise.stat, 2)
        txt <- paste0(prefix, txt)
        legend("center", legend = txt, box.col = "white", cex = 1.5)
    }

    # draw pairs plot:
    par(mar = c(4, 4, 2, 1))

    pairs(mat,
          log = "xy",
          upper.panel = points,
          lower.panel = print.sd.log.ratio,
          oma = c(3, 3, 3, 35),
          col = col,
          pch = pch,
          cex = cex,
          labels = colnames(mat),
          main = main
    )
    # draw a legend:
    if (length(collegend) > 0) {
        legend("right",
               col = c(NA, collegend, NA, NA),
               legend = c(legend.main, names(collegend), "", "Numbers show SD(log2(ratios))"),
               pch = 16
        )
    }
}

#' Generate Protein QC signal boxplot figure
#' 
#' @description 
#' For use with protein data ONLY.
#' 
#' @param object name of the object class to subset
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param neg.names names of IgGs, if NULL IgGs will be detected automatically
#' 
#' @examples
#' testData <- readRDS(file= system.file("extdata","DSP_Proteogenomics_Example_Data", 
#' "proteinData.rds", package = "GeomxTools"))
#' 
#' proteinData <- analyteSubset(object = aggTestData, analyte = "protein")
#' 
#' igg.names <- igg_names(proteinData)
#' 
#' qc_protein_signal(object = proteinData, neg.names = igg.names)
#' 
#' @export

qc_protein_signal <- function(object, neg.names=NULL) {
    
    if(any(fData(object)$AnalyteType == "Protein") == FALSE){
      stop("This figure is only meant for protein data")
    }  
  
    if(is.null(neg.names)){
      neg.names <- igg_names(object)
    }
  
    object <- object[fData(object)$AnalyteType == "Protein",]
    raw <- exprs(object)
  
    # estimate background:
    negfactor <- apply(raw[neg.names, , drop = FALSE], 2, function(x){pmax(ngeoMean(x), 1)})

    # calc snr
    snr <- sweep(raw, 2, negfactor, "/")

    igginds <- which(is.element(rownames(snr), neg.names))
    o <- c(igginds, setdiff(order(apply(snr, 1, median)), igginds))

    protnames <- rownames(snr)
    
    par(mar = c(11, 4, 2, 1))
    boxplot(t(log2(snr[o, ])),
            las = 2,
            outline = FALSE,
            ylim = range(log2(snr+1)),
            names = protnames[o],
            ylab = "Log2 signal-to-background ratio",
            cex.axis = .85 - 0.3 * (nrow(snr) > 60)
    )
    axis(2, at = 1, labels = 1, las = 2, cex = 0.5)
    points(jitter(rep(1:nrow(snr), ncol(snr))),
           log2(snr[o, ]),
           col = "#00008B80",
           pch = 16, cex = 0.5
    )
    abline(h = 0)
    abline(v = length(igginds) + 0.5, lty = 2)
    abline(h = 1, lty = 2)
}

#' Generate concordance figure of targets based on user provided factors
#' 
#' @description  
#' Upper panels are the concordance plot.
#' Lower panels are the standard deviation of the log2-ratios between the targets 
#' 
#' @param targetList names of targets to plot concordance, normally IgGs. 
#' @param object name of the object class to subset
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param plot_factors segment factor to color the plot by
#' 
#' @examples
#' testData <- readRDS(file= system.file("extdata","DSP_Proteogenomics_Example_Data", 
#' "proteinData.rds", package = "GeomxTools"))
#' 
#' proteinData <- analyteSubset(object = aggTestData, analyte = "protein")
#' RNAData <- analyteSubset(object = aggTestData, analyte = "RNA")
#' 
#' igg.names <- igg_names(proteinData)
#' 
#' plot_concordance(targetList = igg.names, object = proteinData, plot_factors = "Segment_Type)
#' plot_concordance(targetList = c("C1orf43", "GPI", "OAZ1"), object = RNAData, plot_factors = "Segment_Type)
#' 
#' @export

plot_concordance <- function(targetList, object, plot_factors){
  
  if(!plot_factors %in% colnames(sData(object))){
    stop("Given plot_factors are not in dataset, spelling and capitalization matter")
  }
  
  cols <- assign_colors(annot = sData(object)[, plot_factors, drop = FALSE])
  
  if (length(targetList) > 1) {
    par(mar = c(4, 4, 4, 1))
    for (varname in plot_factors) {
      tempmat <- t(pmax(exprs(object)[targetList, ], 1))
      colnames(tempmat) <- paste0(colnames(tempmat), " counts")
      concordance_plot(
        mat = tempmat,
        col = cols[[varname]][as.character(sData(object)[, varname])],
        collegend = cols[[varname]],
        legend.main = varname
      )
    }
  }
}


#' Generate concordance figure of normalization factors based on user provided factors
#' 
#' @description 
#' For use with protein data ONLY.
#' 
#' Upper panels are the concordance plot.
#' Lower panels are the standard deviation of the log2-ratios between the normalization factors 
#' 
#' @param object name of the object class to subset
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param plot_factors segment factor to color the plot by
#' @param normfactors normalization factors from compute_normalization_factors(). If NULL these are calculated automatically. 
#' 
#' @examples
#' testData <- readRDS(file= system.file("extdata","DSP_Proteogenomics_Example_Data", 
#' "proteinData.rds", package = "GeomxTools"))
#' 
#' proteinData <- analyteSubset(object = aggTestData, analyte = "protein")
#' 
#' plot_normFactor_concordance(object = proteinData, plot_factors = "Segment_Type")
#' 
#' @export

plot_normFactor_concordance <- function(object, plot_factors, normfactors = NULL){
  if(!plot_factors %in% colnames(sData(object))){
    stop("Given plot_factors are not in dataset, spelling and capitalization matter")
  }
  
  cols <- assign_colors(annot = sData(object)[, plot_factors, drop = FALSE])
  
  if(is.null(normfactors)){
    normfactors <- compute_normalization_factors(object)
  }
  
  if (ncol(normfactors) > 1) {
      par(mar = c(4, 4, 2, 1))
      # pairs plots:
      for (varname in names(cols)) {
          tempmat <- pmax(normfactors, 1)
          colnames(tempmat)[colnames(tempmat) == "HK geomean"] <- "HK geomean\n(counts)"
          colnames(tempmat)[colnames(tempmat) == "Neg geomean"] <- "Neg geomean\n(counts)"
          colnames(tempmat)[colnames(tempmat) == "Nuclei"] <- "Nuclei\n(counts)"
          colnames(tempmat)[colnames(tempmat) == "Area"] <- "Area (microns squared)"

          concordance_plot(
              mat = tempmat,
              col = cols[[varname]][as.character(sData(object)[, varname])],
              collegend = cols[[varname]],
              legend.main = varname,
              main = "Normalization factors"
          )
      }
  }
}

