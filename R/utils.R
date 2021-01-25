ngeoMean <- function(v) {
    v[v == 0] <- 1
    return(EnvStats::geoMean(v, na.rm = TRUE))
}

aggregateCounts <- function(probeCounts) {

  for (index in seq_len(length(data))) {
    countMat <- data[[index]]$Code_Summary
    #NEO this should be a check and we should stop if not all targets found in pkc?
    countMat <- countMat[which(rownames(countMat) %in% rownames(pkcData)), , drop = FALSE]
    countMat$Target <- pkcData$Target[match(countMat$RTS_ID, 
                                        pkcData$RTS_ID)]
    
    countMat$Pool <- pkcData$Module[match(countMat$RTS_ID, 
                                          pkcData$RTS_ID)]
    data[[index]]$Code_Summary <- countMat
  }
  
    targCounts <- reshape2::dcast(probeCounts, Target + Pool ~ Sample_ID, 
        value.var = 'Count', fun.aggregate = ngeoMean, fill = 1)

    if (length(unique(targCounts$Target)) != nrow(targCounts)) {
        duplicatedTargets <- 
            targCounts$Target[which(duplicated(targCounts$Target))]
        warning(
            sprintf('Some targets are listed in multiple pools including %s', 
            paste0(duplicatedTargets, collapse = ",")))
        for (duplicatedTarget in duplicatedTargets) {
            targCounts$Target[which(targCounts$Target == duplicatedTarget)] <- 
                sapply(which(targCounts$Target == duplicatedTarget), 
                    function(index) {
                        paste0(targCounts[index, c("Target", "Pool")], 
                            collapse = "_")
                    })
        }
    }

    #NEO modify this portion still
    rownames(targetAssay) <- targetAssay[, "Target"]
    assay <- as.matrix(targetAssay[, -seq_len(2)])
  
    # Create featureData
    feature <- targetAssay[, "Target", drop = FALSE]
    rownames(feature) <- feature[["Target"]]
  
    feature <- AnnotatedDataFrame(feature,
                                  dimLabels = c("featureNames", "featureColumns"))
}