library(GeomxTools)
library(testthat)
library(EnvStats)

testData <- 
    readRDS(file= system.file("extdata", "DSP_NGS_Example_Data", 
                              "demoData.rds", package = "GeomxTools"))
testData <- (shiftCountsOne(testData, elt="exprs", useDALogic=TRUE))

subData <- 
    subset(testData, subset=Module == "VnV_GeoMx_Hs_CTA_v1.2", 
           select=sampleNames(testData) %in% 
                      sample(sampleNames(testData), 
                  10, replace=FALSE))
subData <- 
    subset(subData, 
           subset=TargetName %in% 
                      sample(unique(fData(subData)[["TargetName"]]), 
                  10, replace=FALSE))

subAggd <- suppressWarnings(aggregateCounts(subData))

testthat::test_that("Feature type changed after aggregation", {
    expect_true(featureType(subData) == "Probe")
    expect_true(featureType(subAggd) == "Target")
})

testthat::test_that("Aggregated object has target dimension and annotations", {
    expect_true(all(fData(subData)[["TargetName"]] %in% featureNames(subAggd)))
    expect_true(length(unique(fData(subData)[["TargetName"]])) == 
                    dim(subAggd)[[1L]])
    expect_true(dim(subData)[[2L]] == dim(subAggd)[[2L]])

    targetLabels <- 
        fvarLabels(subData)[!fvarLabels(subData) %in% 
                                c("RTS_ID", "QCFlags", "ProbeID")]
    expect_true(all(targetLabels %in% fvarLabels(subAggd)))
    expect_true(all(svarLabels(subData) %in% svarLabels(subAggd)))
})

testthat::test_that("Target expression matrix contains aggregated counts", {
    expect_true(all(colnames(exprs(subData)) == colnames(exprs(subAggd))))
    sameTargs <- intersect(fData(subData)[["TargetName"]],
                           rownames(exprs(subAggd)))
    expect_true(length(sameTargs) == 
                    length(unique(fData(subData)[["TargetName"]])))
    subList <- lapply(featureNames(subAggd), function(currTarg) {
        aggdCounts <- exprs(subAggd)[currTarg, sampleNames(subData)]
        targData <- exprs(subset(subData, subset=TargetName == currTarg))
        targMeans <- apply(targData, 2L, function(sampCount) {
            ifelse(length(sampCount) > 1, ngeoMean(sampCount), sampCount)})
        return(all(aggdCounts == targMeans))
    })
    expect_true(all(unlist(subList)))
})

testthat::test_that("Other aggregation functions work", {
    subSum <- suppressWarnings(aggregateCounts(subData, FUN=sum))
    expect_true(all(colnames(exprs(subData)) == colnames(exprs(subSum))))
    sameTargs <- intersect(fData(subData)[["TargetName"]],
                           rownames(exprs(subSum)))
    expect_true(length(sameTargs) == 
                    length(unique(fData(subData)[["TargetName"]])))
    subList <- lapply(featureNames(subSum), function(currTarg) {
        aggdCounts <- exprs(subSum)[currTarg, sampleNames(subData)]
        targData <- exprs(subset(subData, subset=TargetName == currTarg))
        targSums <- apply(targData, 2L, function(sampCount) {
            ifelse(length(sampCount) > 1, sum(sampCount), sampCount)})
        return(all(aggdCounts == targSums))
    })
    expect_true(all(unlist(subList)))
})