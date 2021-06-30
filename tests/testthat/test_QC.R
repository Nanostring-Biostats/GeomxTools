# devtools::install_github("https://github.com/Nanostring-Biostats/GeomxTools/", ref="dev", force = TRUE, dependencies = FALSE)

library(NanoStringNCTools)
library(GeomxTools)
library(testthat)
library(EnvStats)
library(ggiraph)


# This dataset is the same as what is used in Analyzing-GeoMx-NGS-Data-with-GeomxTools.html. 
# The original file is located at: datadir <- file.path( "/home" , "rstudio" , "NAS_data", "rvitancol", "kidney_demo")

datadir <- file.path( "~/NAS_data", "rvitancol", "kidney_demo")
# datadir <- file.path("~/YREN", "GeomxTools", "GeoMx_SampleData")
DCCFiles <- list.files(file.path( datadir , "DCC_files"), pattern=".dcc$", full.names=TRUE)
PKCFiles <- list.files(file.path(datadir), pattern=".pkc$", full.names=TRUE)
SampleAnnotationFile <- file.path(datadir, "kidney_demo_AOI_Annotations.xlsx")

testData <-
  suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                          pkcFiles = PKCFiles,
                                          phenoDataFile = SampleAnnotationFile,
                                          phenoDataSheet = "Template",
                                          phenoDataDccColName = "Sample_ID",
                                          protocolDataColNames = c("aoi", 
                                                                   "roi",
                                                                   "slide name"),
                                          experimentDataColNames = c("panel", "instrument_type")))

#Shift counts to one to mimic how DSPDA handles zero counts
testData <- shiftCountsOne(testData, elt="exprs", useDALogic=TRUE) 



##### Technical Signal QC #####

testData <- setSeqQCFlags(testData, 
                          qcCutoffs=list(minSegmentReads=1000, 
                                         percentAligned=80, 
                                         percentSaturation=50))
prData <- protocolData(testData)
TechSigQC <- as.data.frame(prData[["QCFlags"]])

# req 1: test that the number of Low Reads is correct:------
testthat::test_that("test that the number of Low Reads is correct", {
  expect_true(sum(TechSigQC$LowReads) == sum(prData@data$Raw < 1000))
})


# req 2: test that the number of Low Percent Trimmed is correct:------
testthat::test_that("test that the number of Low Percent Trimmed is correct", {
  expect_true(sum(TechSigQC$LowTrimmed) == sum(prData@data$`Trimmed (%)` < 80))
})


# req 3: test that the number of Low Percent Stiched is correct:------
testthat::test_that("test that the number of Low Percent Stiched is correct", {
  expect_true(sum(TechSigQC$LowStitched) == sum(prData@data$`Stitched (%)` < 80))
})


# req 4: test that the number of Low Percent Aligned is correct:------
testthat::test_that("test that the number of Low Percent Aligned is correct", {
  expect_true(sum(TechSigQC$LowAligned) == sum(prData@data$`Aligned (%)` < 80))
})


# req 5: test that the number of Low Percent Saturation is correct:------
testthat::test_that("test that the number of Low Percent Saturation is correct", {
  expect_true(sum(TechSigQC$LowSaturation) == sum(prData@data$`Saturated (%)` < 50))
})




##### Technical Background QC #####

neg_set <- exprs(negativeControlSubset(testData))

testData <- setBackgroundQCFlags(testData, 
                                 qcCutoffs=list(minNegativeCount=10, 
                                                maxNTCCount=60))


prData <- protocolData(testData)
TechBgQC <- as.data.frame(prData[["QCFlags"]])

# req 1: test that the number of Low Negatives is correct:------
testthat::test_that("test that the number of Low Negatives is correct", {
  expect_true(sum(TechBgQC$LowNegatives) == sum(prData@data$NegGeoMean < 10)) 
})


# req 2: test that the number of High NTC is correct:------
testthat::test_that("test that the number of High NTC is correct", {
  expect_true(sum(TechBgQC$HighNTC) == sum(prData@data$NTC > 60))
})




##### Segment QC #####

testData <- setGeoMxQCFlags(testData, 
                            qcCutoffs=list(minNuclei=16000,
                                           minArea=20))
prData <- protocolData(testData)
segQC <- as.data.frame(prData[["QCFlags"]])

# req 1: test that the number of Low Area is correct:------
testthat::test_that("test that the number of Low Area is correct", {
  expect_true(sum(segQC$LowArea) == sum(sData(testData)[c("area")] < 20)) 
})


# req 2: test that the number of Low Nuclei is correct:------
testthat::test_that("test that the number of Low Nuclei is correct", {
  expect_true(sum(segQC$LowNuclei) == sum(testData@phenoData@data$nuclei < 16000))
})



##### Biological Probe QC #####

testData <- setBioProbeQCFlags(testData, 
                               qcCutoffs=list(minProbeRatio=0.1,
                                              percentFailGrubbs=20))


probeQC <- fData(testData)[["QCFlags"]]
probeQCResults <- data.frame(fData(testData)[["QCFlags"]])

# req 1: test that the number of Low Probe Ratio is correct:------
testthat::test_that("test that the number of Low Probe Ratio is correct", {
  expect_true(sum(probeQC$LowProbeRatio) == sum(fData(testData)["ProbeRatio"] < 0.1)) 
})



# req 2: test that the Low Grubbs Outlier is correct:------ ?????????????
neg_set <- exprs(negativeControlSubset(testData))
# globalOutlier = rownames(probeQCResults)[which(probeQCResults$GrubbsOutlier==TRUE)]

outliers::grubbs.test(neg_set[,1], two.sided = TRUE)
probeQCResults["RTS0039349", "LocalGrubbsOutlier.DSP.1001250007851.H.A02.dcc"]
