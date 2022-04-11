# Return result from sData method and compare to expected
# Return result from svarLabels method and compare to expected
# Return result from dimLabels method and compare to expected


library(GeomxTools)
library(testthat)

datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)

protocolDataColNames <- c("aoi",
                          "cell_line",
                          "roi_rep",
                          "pool_rep",
                          "slide_rep")

testData <- readRDS(file= system.file("extdata", "DSP_NGS_Example_Data", 
                                      "demoData.rds", package = "GeomxTools"))


testData_agg <- aggregateCounts(testData)

DCCFiles <- DCCFiles[!basename(DCCFiles) %in% unique(sData(testData)$NTC_ID)]

# Spec 1: test that the rownames and column names in sData are correct:------
testthat::test_that("test that the rownames and column names in sData are correct", {
  expect_true(all(basename(DCCFiles) %in% rownames(sData(testData))))
  expect_true(all(c(names(testData@phenoData@data),
                    protocolDataColNames) %in% colnames(sData(testData))))
})


# Spec 2: test that the svarLabels method gives the correct results:------
testthat::test_that("test that the svarLabels method gives the correct results", {
  expect_true(all(svarLabels(testData) == colnames(sData(testData))))
})




# Spec 3: test that the dimLabels method gives the correct results:------
testthat::test_that("test that the dimLabels method gives the correct results", {
  expect_true(length(dimLabels(testData)) == 2)
  expect_true(all(paste0(sData(testData)[[dimLabels(testData)[2]]], ".dcc") == colnames(testData@assayData$exprs)))
  expect_true(all(testData@featureData@data[[dimLabels(testData)[1]]] == rownames(testData@assayData$exprs)))
})



# Spec 4: test that the design method gives the correct results:------
testthat::test_that("test that the design method gives the correct results", {
  expect_true(is.null(design(testData)))
})



# Spec 5: test that the featureType method gives the correct results:------
testthat::test_that("test that the featureType method gives the correct results", {
  expect_true(featureType(testData_agg) == "Target")
  expect_true(featureType(testData) == "Probe")
  expect_false(featureType(testData_agg) == featureType(testData))
})

# Spec 6: test that the countsShiftedByOne method gives the correct results:------
testthat::test_that("test that the countsShiftedByOne method gives the correct results", {
  expect_false(countsShiftedByOne(testData_agg))
  expect_false(countsShiftedByOne(testData))
  expect_true(countsShiftedByOne(shiftCountsOne(testData)))
})

proteinData <- readRDS(file= system.file("extdata","DSP_Proteogenomics_Example_Data", 
                                         "proteinData.rds", package = "GeomxTools"))

# Spec 7: test that the analyte method gives the correct results:------
testthat::test_that("test that the analyte method gives the correct results", {
  expect_true(analyte(testData_agg) == "RNA")
  expect_true(analyte(testData) == "RNA")
  expect_true(analyte(testData_agg) == analyte(testData))
  
  expect_true(analyte(proteinData) == "Protein")
  expect_false(analyte(proteinData) == analyte(testData))
})

