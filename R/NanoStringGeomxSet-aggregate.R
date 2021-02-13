

aggregateCounts <- function(object, FUN=ngeoMean) {
    targetCounts <- do.call(rbind, esBy(object, GROUP = "TargetName", 
        FUN=function(x) {esApply(x, 2, FUN)}, simplify=FALSE))
    targetFeats <- featureData(object)@data
    targetFeats <- 
        targetFeats[!duplicated(targetFeats[, c("TargetName", "Module")]), ]
    rownames(targetFeats) <- targetFeats[, "TargetName"]
    targetFeats <- 
        targetFeats[, !colnames(targetFeats) %in% c("RTS_ID", "QCFlags")]
    targetFeats <- 
         AnnotatedDataFrame(targetFeats[rownames(targetCounts), ], 
                            dimLabels = c("featureNames", "featureColumns"))
    targetObject <- NanoStringGeomxSet(
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


