# devtools::install_github("https://github.com/Nanostring-Biostats/GeomxTools/", ref="dev", force = TRUE, dependencies = FALSE)

library(GeomxTools)
library(testthat)


datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
SampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

protocolDataColNames <- c("aoi",
                          "cell_line",
                          "roi_rep",
                          "pool_rep",
                          "slide_rep")

testData <-
  suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
                                          pkcFiles = PKCFiles,
                                          phenoDataFile = SampleAnnotationFile,
                                          phenoDataSheet = "CW005",
                                          phenoDataDccColName = "Sample_ID",
                                          protocolDataColNames = protocolDataColNames,
                                          experimentDataColNames = c("panel")))



DCCFiles <- DCCFiles[!basename(DCCFiles) %in% unique(sData(testData)$NTC_ID)]
testData2 <- shiftCountsOne(testData, elt="exprs", useDALogic=TRUE) 
testData3 <- shiftCountsOne(testData, elt = "exprs", useDALogic = FALSE)



# req 1: test that the counts from shiftCountsOne(..., useDALogic = TRUE) are correct:------

fakeData = testData
testthat::test_that("test that the counts from shiftCountsOne(..., useDALogic = TRUE) are correct", {
  counts_0_add_1 = apply(exprs(fakeData), 1:2, function(x) ifelse(x==0, x+1, x))
  expect_true(identical(counts_0_add_1, exprs(testData2)))
})



# req 2: test that the counts from shiftCountsOne(..., useDALogic = FALSE) are correct:------
testthat::test_that("test that the counts from shiftCountsOne(..., useDALogic = FALSE) are correct", {
  counts_all_add_1 = apply(exprs(testData), 1:2, function(x) x+1)
  expect_true(identical(counts_all_add_1, exprs(testData3)))
})

# req 3: tests that countsShiftedByOne returns true when counts have been shifted:------
testthat::test_that("tests that countsShiftedByOne returns correct value", {
  expect_false(countsShiftedByOne(testData))
  expect_true(countsShiftedByOne(testData2))
  expect_true(countsShiftedByOne(testData3))
})

