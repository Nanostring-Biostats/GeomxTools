# sum(matches) == nrow(testData@assayData$exprs) should result in TRUE

library(NanoStringNCTools)
library(GeomxTools)
library(testthat)
library(EnvStats)

datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
SampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

testData <-
  suppressWarnings(GeomxTools::readNanoStringGeoMxSet(dccFiles = DCCFiles, # QuickBase: readNanoStringGeomxSet, need to change it.
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


PKC <- readPKCFile(PKCFiles)
PKC$RTS_ID <- gsub("RNA", "RTS00", PKC$RTS_ID)
DCC_file <- DCCFiles[2]
DCC <- suppressWarnings(readDccFile(DCC_file))
DCC_file <- basename(DCC_file)
rownames(DCC$Code_Summary) <- gsub("RNA", "RTS00", rownames(DCC$Code_Summary))

matches <- NULL 
for(i in unique(fData(test_Data)[["TargetName"]])){
  probes <- PKC$RTS_ID[which(PKC$Target == i)]
  if(sum(is.na(DCC$Code_Summary[probes,"Count"])) > 0){
    DCC$Code_Summary[probes, "Count"][is.na(DCC$Code_Summary[probes, "Count"])] <- 1
  }
  matches <- c(matches, EnvStats::geoMean(DCC$Code_Summary[probes,"Count"]) == 
                 test_Data@assayData$exprs[i,DCC_file])
}



# req 1: test that the number of collapsed probe is correct:------
testthat::test_that("test that the number of collapsed probe is correct", {
  expect_true(sum(matches) == nrow(test_Data@assayData$exprs))
})



# a = PKC$RTS_ID[which(PKC$Target == "ANGPT2")]
# EnvStats::geoMean(DCC$Code_Summary[a,"Count"])
# EnvStats::geoMean(DCC$Code_Summary[a,"Count"],na.rm = T)
# test_Data@assayData$exprs["ANGPT2",DCC_file]
# DCC$Code_Summary[a,"Count"]
# EnvStats::geoMean(c(9,1,2,1))
# EnvStats::geoMean(c(9,1,2,1,1))