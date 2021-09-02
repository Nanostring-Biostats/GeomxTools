#' Aggregate probe counts to target level for feature data
#' 
#' @param object name of the NanoStringGeoMxSet object to aggregate
#' @param FUN function to use for count aggregation
#' 
#' @return a NanoStringGeoMxSet object with targets as features
#' 
#' @importFrom NanoStringRccSet esBy
#' 
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'                        package="GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' targetGeoMxSet <- aggregateCounts(demoData)
#' 
#' @export
#' 
aggregateCounts <- function(object, FUN=ngeoMean) {
    if(featureType(object) == "Target") {
        stop("GeoMxSet object feature type is already target-level. ",
             "No further aggregation can be performed.")
    }

    # Skip targets with single probe
    multiProbeTable <- with(object, table(TargetName)) > 1L
    multiProbeTargs <- 
        names(multiProbeTable)[which(multiProbeTable, arr.ind=TRUE)]
    if (length(multiProbeTargs) > 0) {
        multiObject <- 
            object[fData(object)[["TargetName"]] %in% multiProbeTargs, ]
        if ("Negative" %in% unique(fData(object)[["CodeClass"]])) {
            object <- summarizeNegatives(object)
        } else {
            warning("Object has no negatives. ",
                "No summary statistics for negatives will be calculated.")
        }
    } else {
        warning("Object has no multiprobe targets. ",
                "No aggregation was performed.")
        featureNames(object) <- fData(object)[["TargetName"]]
        featureType(object) <- "Target"
        return(object)
    }
    
    targetCounts <- do.call(rbind, esBy(multiObject, GROUP = "TargetName", 
        FUN=function(x) {esApply(x, 2, FUN)}, simplify=FALSE))
    singleProbeObject <- subset(object, 
                                subset=!TargetName %in% multiProbeTargs)
    singleProbeCounts <- exprs(singleProbeObject)
    rownames(singleProbeCounts) <- fData(singleProbeObject)[["TargetName"]]
    targetCounts <- rbind(targetCounts, singleProbeCounts)
    targetCounts <- targetCounts[unique(fData(object)[["TargetName"]]), ]
        
    targetFeats <- fData(object)
    targetFeats <- 
        targetFeats[!duplicated(targetFeats[["TargetName"]]), ]
    rownames(targetFeats) <- targetFeats[, "TargetName"]
    probeColumns <- c("RTS_ID", "QCFlags", "ProbeID", 
                      "ProbeRatio", "OutlierFrequency")
    targetFeats <- 
        targetFeats[, !colnames(targetFeats) %in% probeColumns]
    targetFeats <- 
         AnnotatedDataFrame(targetFeats[rownames(targetCounts), ], 
                            dimLabels = c("featureNames", "featureColumns"))

    targetObject <- NanoStringGeoMxSet(
        assayData = targetCounts,
        phenoData = phenoData(object),
        featureData = targetFeats,
        experimentData = experimentData(object),
        annotation = annotation(object),
        protocolData = protocolData(object),
        featureType = "Target",
        check = FALSE)
    return(targetObject)
}

#' Calculate negative probe summary stats
#' 
#' @param object name of the NanoStringGeoMxSet object to summarize
#' @param functionList optional list of additional functions to calculate negative 
#' probe stats, list element names should correspond to expected stat column header
#' 
#' @return a NanoStringGeoMxSet object with negative probe summary stats 
#' appended to sample data
#' 
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'                        package="GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' demoData <- 
#'     summarizeNegatives(demoData, 
#'                        functionList=c(mean=mean, min=min, max=max))
#' 
#' @export
#' 
summarizeNegatives <- 
    function(object, functionList=c()) {
        functionList <- 
            append(c(NegGeoMean=ngeoMean, NegGeoSD=ngeoSD), functionList)
        negObject <- negativeControlSubset(object)
        summaryList <- 
            lapply(functionList, function(currFunc) {
                esBy(negativeControlSubset(object), 
                     GROUP="Module", 
                     FUN=function(x) { 
                         assayDataApply(x, 
                                        MARGIN = 2,
                                        FUN=currFunc,
                                        elt="exprs") 
                     })
            })
        summaryListNames <- 
            unlist(lapply(names(summaryList), 
                       function (x) {
                           paste0(x, "_", colnames(summaryList[[x]]))
                        }))
        summaryDF <- do.call(cbind, summaryList)
        colnames(summaryDF) <- summaryListNames
        summaryDF <- summaryDF[sampleNames(object), ]
        pData(object) <- cbind(pData(object), summaryDF)
        return(object)
}
