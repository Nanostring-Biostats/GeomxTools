# sum(matches) == nrow(testData@assayData$exprs) should result in TRUE

library(GeomxTools)
library(testthat)
library(EnvStats)

datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
SampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

testData <-
  suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles, # QuickBase: readNanoStringGeomxSet, need to change it.
                                          pkcFiles = PKCFiles,
                                          phenoDataFile = SampleAnnotationFile,
                                          phenoDataSheet = "CW005",
                                          phenoDataDccColName = "Sample_ID",
                                          protocolDataColNames = c("aoi",
                                                                   "cell_line",
                                                                   "roi_rep",
                                                                   "pool_rep",
                                                                   "slide_rep"),
                                          experimentDataColNames = c("panel")))


test_Data <- aggregateCounts(testData)

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
  
  for(i in unique(fData(test_Data)[["TargetName"]])){
    probes <- PKC$RTS_ID[which(PKC$Target == i)]
    if(sum(is.na(DCC$Code_Summary[probes,"Count"])) > 0){
      #NAs are changed to 0 counts
      DCC$Code_Summary[probes, "Count"][is.na(DCC$Code_Summary[probes, "Count"])] <- 0
      
      #0 count doesn't meet minimum threshold so all counts are increased by threshold
      #0.5 is default threshold from thresholdValues()
      DCC$Code_Summary[probes, "Count"] <- DCC$Code_Summary[probes, "Count"]+0.5
    }
    
    matches <- c(matches, EnvStats::geoMean(DCC$Code_Summary[probes,"Count"]) == 
                   test_Data@assayData$exprs[i,DCC_file])
  }
  
  dcc <- dcc + 1
}




# req 1: test that the number of collapsed probe is correct:------
testthat::test_that("test that the number of collapsed probe is correct", {
  expect_true(sum(matches) == nrow(test_Data@assayData$exprs)*numDCC)
})



# a = PKC$RTS_ID[which(PKC$Target == "ANGPT2")]
# EnvStats::geoMean(DCC$Code_Summary[a,"Count"])
# EnvStats::geoMean(DCC$Code_Summary[a,"Count"],na.rm = T)
# test_Data@assayData$exprs["ANGPT2",DCC_file]
# DCC$Code_Summary[a,"Count"]
# EnvStats::geoMean(c(9,1,2,1))
# EnvStats::geoMean(c(9,1,2,1,1))