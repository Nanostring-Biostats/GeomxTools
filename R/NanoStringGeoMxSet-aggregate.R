#' Aggregate probe counts to target level for feature data
#' 
#' @param object name of the NanoStringGeoMxSet object to aggregate
#' @param FUN function to use for count aggregation
#' 
#' @return a NanoStringGeoMxSet object with targets as features
#' 
#' @examples
#' targetGeoMxSet <- aggregateCounts(demoData)
#' 
#' @export
#' 
aggregateCounts <- function(object, FUN=ngeoMean) {
    object <- summarizeNegatives(object)
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
#' demoData <- summarizeNegatives(demoData, functionList = c(mean, min, max))
#' 
#' @export
#' 
summarizeNegatives <- 
    function(object, functionList=c()) {
        functionList <- 
            append(c(NegGeoMean=ngeoMean, NegGeoSD=ngeoSD), functionList)
        negObject <- object[fData(object)[, "Negative"], ]
        summaryList <- 
            lapply(functionList, function(x) {
                esApply(negObject, MARGIN=2, FUN=x)})
        summaryDF <- do.call(cbind, summaryList)
        colnames(summaryDF) <- names(functionList)
        summaryDF <- summaryDF[sampleNames(object), ]
        pData(object) <- cbind(pData(object), summaryDF)
        return(object)
}
