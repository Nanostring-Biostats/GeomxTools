library(GeomxTools)
library(testthat)
library(EnvStats)

datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))

testData <- readRDS(file= system.file("extdata", "DSP_NGS_Example_Data", 
                            "demoData.rds", package = "GeomxTools"))

aggTestData <- aggregateCounts(testData)

DCCFiles <- DCCFiles[!basename(DCCFiles) %in% unique(sData(testData)$NTC_ID)]

PKC <- readPKCFile(PKCFiles)
PKC$RTS_ID <- gsub("RNA", "RTS00", PKC$RTS_ID)

numDCC <- 10

#random subset of 10 DCC files
DCCFiles <- DCCFiles[sample(1:length(DCCFiles), numDCC)]

matches <- NULL
dcc <- 1
while(dcc <= length(DCCFiles) & all(matches == TRUE)){
  DCC_file <- DCCFiles[dcc]
  DCC <- suppressWarnings(readDccFile(DCC_file))
  DCC_file <- basename(DCC_file)
  rownames(DCC$Code_Summary) <- gsub("RNA", "RTS00", rownames(DCC$Code_Summary))
  
  for(i in unique(fData(aggTestData)[["TargetName"]])){
    probes <- PKC$RTS_ID[which(PKC$Target == i)]
    if(sum(is.na(DCC$Code_Summary[probes,"Count"])) > 0){
      #NAs are changed to 0 counts
      DCC$Code_Summary[probes, "Count"][is.na(DCC$Code_Summary[probes, "Count"])] <- 0

      #0 count doesn't meet minimum threshold so all counts are increased by threshold
      #0.5 is default threshold from thresholdValues()
      if(length(probes) != 1){
        DCC$Code_Summary[probes, "Count"] <- DCC$Code_Summary[probes, "Count"]+0.5
      }
    }
    
    if(length(probes) == 1){
      matches <- c(matches, DCC$Code_Summary[probes,"Count"] == 
                     aggTestData@assayData$exprs[i,DCC_file])
    }else{
      matches <- c(matches, EnvStats::geoMean(DCC$Code_Summary[probes,"Count"]) == 
                     aggTestData@assayData$exprs[i,DCC_file])
    }
  }
  
  dcc <- dcc + 1
}


# req 1: test that collapsed probe count is equal to the geomean of probes for each target:------
testthat::test_that("test that collapsed probe count is equal to the geomean of probes for each target", {
  expect_true(sum(matches) == nrow(aggTestData@assayData$exprs)*numDCC)
})

matches <- NULL
dcc <- 1
while(dcc <= length(DCCFiles) & all(matches == TRUE)){
  DCC_file <- DCCFiles[dcc]
  DCC <- suppressWarnings(readDccFile(DCC_file))
  DCC_file <- basename(DCC_file)
  rownames(DCC$Code_Summary) <- gsub("RNA", "RTS00", rownames(DCC$Code_Summary))
  
  for(i in unique(PKC$Module)){
    negs <- PKC$RTS_ID[which(PKC$CodeClass == "Negative" & PKC$Module == i)]
    
    if(sum(is.na(DCC$Code_Summary[negs,"Count"])) > 0){
      #NAs are changed to 0 counts
      DCC$Code_Summary[negs, "Count"][is.na(DCC$Code_Summary[negs, "Count"])] <- 0
      
      #0 count doesn't meet minimum threshold so all counts are increased by threshold
      #0.5 is default threshold from thresholdValues()
      if(length(negs) != 1){
        DCC$Code_Summary[negs, "Count"] <- DCC$Code_Summary[negs, "Count"]+0.5
      }
    }
    
    matches <- c(matches, EnvStats::geoMean(DCC$Code_Summary[negs,"Count"]) == 
                   pData(aggTestData)[DCC_file, paste0("NegGeoMean_",i)])
    
    matches <- c(matches, EnvStats::geoSD(DCC$Code_Summary[negs,"Count"]) == 
                   pData(aggTestData)[DCC_file, paste0("NegGeoSD_",i)])
  }
  
  dcc <- dcc + 1
}


# req 2: test that geomean and geosd of negatives is correct
testthat::test_that("test that the geomean and geosd of negatives is correct", {
  #geomean and geosd (2) for each DCC (numDCC) * number of modules
  expect_true(sum(matches) == length(unique(PKC$Module))*numDCC*2)
})


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

# req 3: featureType is changed after aggregation
testthat::test_that("Feature type changed after aggregation", {
    expect_true(featureType(subData) == "Probe")
    expect_true(featureType(subAggd) == "Target")

})

# req 4: Aggregated object has target dimension and annotations
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

# req 5: Target expression matrix contains aggregated counts
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

# req 6: Other aggregation functions work
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


# req 7: Counts can't be aggregated twice
testthat::test_that("Counts can't be aggregated twice", {
  expect_error(aggregateCounts(aggTestData))
})

oneProbe <- table(fData(testData)$TargetName)
oneProbe <- names(oneProbe)[which(oneProbe == 1)]
oneProbeIDs <- fData(testData)$RTS_ID[match(oneProbe,fData(testData)$TargetName)]

# req 8: Single Probes are not aggregated
testthat::test_that("Single Probes are not aggregated", {
  expect_true(all(aggTestData@assayData$exprs[oneProbe,] == testData@assayData$exprs[oneProbeIDs,]))
})

negs <- which(fData(testData)$CodeClass == "Negative")
# req 9: Warning given when data doesn't have negatives
testthat::test_that("Warning given when data doesn't have negatives", {
  expect_warning(aggregateCounts(testData[-negs,]))
})

# req 10: Warning given when data doesn't have multi-probe targets
testthat::test_that("Warning given when data doesn't have multi-probe targets", {
  expect_warning(aggregateCounts(testData[!duplicated(fData(testData)$TargetName),]))
})

negs <- summarizeNegatives(testData)
negs2 <- summarizeNegatives(negs)

# req 11: summarizeNegatives can be rerun
testthat::test_that("summarizeNegatives can be rerun", {
    expect_identical(negs, negs2)
})

