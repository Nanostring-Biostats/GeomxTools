# dimLabels(testData) should result in "A" "B"
# dimLabels(testData) == labs should result in FALSE FALSE
# design(testData) should result in x ~ y
# design(testData) == des should result in logical(0)


library(GeomxTools)
library(testthat)


datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
SampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

testData <-
  suppressWarnings(GeomxTools::readNanoStringGeoMxSet(dccFiles = DCCFiles, 
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




# req 1: test that dimLabels(testData) matches the value assigned:------
testthat::test_that("test that dimLabels(testData) matches the value assigned", {
  labs <- dimLabels(testData)
  dimLabels(testData) <- c("A", "B")
  expect_true(all(dimLabels(testData) == c("A", "B")))
  expect_false(all(dimLabels(testData) == labs))
})




# req 2: test that design(testData) matches the value assigned:------
testthat::test_that("test that dimLabels(testData) matches the value assigned", {
  des <- design(testData)
  expect_null(des)
  design(testData) <- "x ~ y"
  expect_true(all(design(testData) == c("x ~ y")))
  expect_false(identical(design(testData), des))
})

# req 3: test that featureType(testData) matches the value assigned:------
testthat::test_that("test that dimLabels(testData) matches the value assigned", {
  featType <- featureType(testData)
  featureType(testData) <- "Target"
  expect_true(all(featureType(testData) == "Target"))
  expect_false(identical(featureType(testData), featType))
})

