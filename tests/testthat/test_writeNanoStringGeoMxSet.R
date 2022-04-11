### CONFIGURATION SECTION ###
library(GeomxTools)
library(testthat)

datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
SampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

demoData <- readRDS(file= system.file("extdata", "DSP_NGS_Example_Data", 
                                      "demoData.rds", package = "GeomxTools"))


writeDir <-  "writeTest/"
if(dir.exists(writeDir)){
  unlink(writeDir, recursive = TRUE, force = TRUE)
}

NTCs <- unique(demoData@protocolData@data$NTC_ID)
writtenDCCs <- basename(DCCFiles)
writtenDCCs <- writtenDCCs[!writtenDCCs %in% NTCs]

writeNanoStringGeoMxSet(demoData, dir = writeDir)

# req 1: test DCC files are written after read in
testthat::test_that("test DCC files are written after read in", {
  expect_true(dir.exists(writeDir))
  expect_true(all(dir(writeDir) %in% writtenDCCs))
  expect_true(all(writtenDCCs %in% dir(writeDir)))
})

testData <-
  suppressWarnings(readNanoStringGeoMxSet(dccFiles = dir(writeDir, full.names = TRUE), 
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

# req 2: test GeomxSet object from written DCC is identical to original
testthat::test_that("test GeomxSet object from written DCC is identical to original", {
  expect_identical(pData(testData), pData(demoData))
  expect_identical(fData(testData), fData(demoData))
  expect_identical(sData(testData)[,-which(colnames(sData(testData)) == "NTC")], 
                   sData(demoData)[,-which(colnames(sData(demoData)) == "NTC")])
  expect_identical(exprs(testData), exprs(demoData))
  expect_identical(testData@experimentData, demoData@experimentData)
})

rm(testData)
unlink(writeDir, recursive = TRUE, force = TRUE)

#Shift counts to one to mimic how DSPDA handles zero counts
demoData <- shiftCountsOne(demoData, elt="exprs", useDALogic=TRUE) 
writeNanoStringGeoMxSet(demoData, dir = writeDir)

# req 3: test DCC files are written after shifting counts by one
testthat::test_that("test DCC files are written after shifting counts by one", {
  expect_true(dir.exists(writeDir))
  expect_true(all(dir(writeDir) %in% writtenDCCs))
  expect_true(all(writtenDCCs %in% dir(writeDir)))
})
unlink(writeDir, recursive = TRUE, force = TRUE)

#QC flags
demoData <- setSegmentQCFlags(demoData)
writeNanoStringGeoMxSet(demoData, dir = writeDir)

# req 4: test DCC files are written after SegmentQC
testthat::test_that("test DCC files are written after SegmentQC", {
  expect_true(dir.exists(writeDir))
  expect_true(all(dir(writeDir) %in% writtenDCCs))
  expect_true(all(writtenDCCs %in% dir(writeDir)))
})
unlink(writeDir, recursive = TRUE, force = TRUE)

demoData <- setBioProbeQCFlags(demoData)
writeNanoStringGeoMxSet(demoData, dir = writeDir)

# req 5: test DCC files are written after Probe QC
testthat::test_that("test DCC files are written after Probe QC", {
  expect_true(dir.exists(writeDir))
  expect_true(all(dir(writeDir) %in% writtenDCCs))
  expect_true(all(writtenDCCs %in% dir(writeDir)))
})
unlink(writeDir, recursive = TRUE, force = TRUE)

#aggregate counts
demoData <- aggregateCounts(demoData)


# req 6: test error occurs when writing set after aggregating counts
testthat::test_that("test error occurs when writing set after aggregating counts", {
  expect_error(writeNanoStringGeoMxSet(demoData, dir = writeDir))
})
unlink(writeDir, recursive = TRUE, force = TRUE)

#normalize
demoData <- normalize(demoData , norm_method="quant",
                      desiredQuantile = .9, toElt = "q_norm")

# req 7: test error occurs when writing set after aggregating counts
testthat::test_that("test error occurs when writing set after manipulating target level counts", {
  expect_error(writeNanoStringGeoMxSet(demoData, dir = writeDir))
})
unlink(writeDir, recursive = TRUE, force = TRUE)

