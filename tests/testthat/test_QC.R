# devtools::install_github("https://github.com/Nanostring-Biostats/GeomxTools/", ref="dev", force = TRUE, dependencies = FALSE)

library(NanoStringNCTools)
library(GeomxTools)
library(testthat)
library(EnvStats)
library(ggiraph)


testDataRaw <- readRDS(file= system.file("extdata","DSP_NGS_Example_Data", "demoData.rds", package = "GeomxTools"))

#Shift counts to one to mimic how DSPDA handles zero counts
testData <- shiftCountsOne(testDataRaw, elt="exprs", useDALogic=TRUE) 

##### Technical Signal QC #####

testData <- setSeqQCFlags(testDataRaw, 
                          qcCutoffs=list(minSegmentReads=1000, 
                                         percentAligned=80, 
                                         percentSaturation=50))
prData <- protocolData(testData)
TechSigQC <- as.data.frame(prData[["QCFlags"]])

# Spec 1: test that the number of Low Reads is correct:------
testthat::test_that("test that the number of Low Reads is correct", {
  expect_true(sum(TechSigQC$LowReads) == sum(prData@data$Raw < 1000))
})


# Spec 2: test that the number of Low Percent Trimmed is correct:------
testthat::test_that("test that the number of Low Percent Trimmed is correct", {
  expect_true(sum(TechSigQC$LowTrimmed) == sum(prData@data$`Trimmed (%)` < 80))
})


# Spec 3: test that the number of Low Percent Stiched is correct:------
testthat::test_that("test that the number of Low Percent Stiched is correct", {
  expect_true(sum(TechSigQC$LowStitched) == sum(prData@data$`Stitched (%)` < 80))
})


# Spec 4: test that the number of Low Percent Aligned is correct:------
testthat::test_that("test that the number of Low Percent Aligned is correct", {
  expect_true(sum(TechSigQC$LowAligned) == sum(prData@data$`Aligned (%)` < 80))
})


# Spec 5: test that the number of Low Percent Saturation is correct:------
testthat::test_that("test that the number of Low Percent Saturation is correct", {
  expect_true(sum(TechSigQC$LowSaturation) == sum(prData@data$`Saturated (%)` < 50))
})

##### Technical Background QC #####

neg_set <- summarizeNegatives(testDataRaw)
neg_set <- pData(neg_set)[grep("NegGeoMean", colnames(pData(neg_set)))]

testData <- setBackgroundQCFlags(testDataRaw, 
                                 qcCutoffs=list(minNegativeCount=10, 
                                                maxNTCCount=60))


testData <- summarizeNegatives(testData)

prData <- protocolData(testData)
TechBgQC <- as.data.frame(prData[["QCFlags"]])


# Spec 6: test that the number of Low Negatives is correct:------
testthat::test_that("test that the number of Low Negatives is correct", {

  expect_true(sum(TechBgQC$LowNegatives) == 
              sum(apply(neg_set < 10, 1, sum) > 0)) 

    numberFail <- (pData(testData)[, "NegGeoMean_VnV_GeoMx_Hs_CTA_v1.2"] < 10) + 
                  (pData(testData)[, "NegGeoMean_Six-gene_test_v1_v1.1"] < 10)
    expect_true(sum(TechBgQC$LowNegatives) ==  sum(numberFail > 0)) 

})


# Spec 7: test that the number of High NTC is correct:------
testthat::test_that("test that the number of High NTC is correct", {
  expect_true(sum(TechBgQC$HighNTC) == sum(prData@data$NTC > 60))
})

##### Segment QC #####

testData <- setGeoMxQCFlags(testDataRaw, 
                            qcCutoffs=list(minArea=16000,
                                           minNuclei=20))
prData <- protocolData(testData)
segQC <- as.data.frame(prData[["QCFlags"]])

# Spec 8: test that the number of Low Area is correct:------
testthat::test_that("test that the number of Low Area is correct", {
  expect_true(sum(segQC$LowArea) == sum(sData(testData)[c("area")] < 16000)) 
})


# Spec 9: test that the number of Low Nuclei is correct:------
testthat::test_that("test that the number of Low Nuclei is correct", {
  expect_true(sum(segQC$LowNuclei) == sum(testData@phenoData@data$nuclei < 20))
})

noArea <- testDataRaw
pData(noArea) <- pData(noArea)[,1:4]

# Spec 10: test that no error occurs when no nuclei or area:------
testthat::test_that("test that no error or QC flags occur when no nuclei or area are in data", {
  expect_identical(setGeoMxQCFlags(noArea, 
                                   qcCutoffs=list(minArea=16000,
                                                  minNuclei=20)),
                   noArea)
})

##### Segment QC Flags #####

testData <- setSegmentQCFlags(testDataRaw, qcCutoffs=list(minNuclei=20,
                                                      minArea=16000,
                                                      minNegativeCount=10, 
                                                      maxNTCCount=60,
                                                      minSegmentReads=1000, 
                                                      percentAligned=80, 
                                                      percentSaturation=50))

prData <- protocolData(testData)
segQC_all <- as.data.frame(prData[["QCFlags"]])

# Spec 18: test that flags match after setSegmentQCFlags:------
testthat::test_that("test that the setGeoMxQCFlags match after setSegmentQCFlags", {
  expect_true(all(segQC_all[,colnames(segQC)] == segQC))
})

testthat::test_that("test that the setBackgroundQCFlags match after setSegmentQCFlags", {
  expect_true(all(segQC_all[,colnames(TechBgQC)] == TechBgQC))
})

testthat::test_that("test that the setSeqQCFlags match after setSegmentQCFlags", {
  expect_true(all(segQC_all[,colnames(TechSigQC)] == TechSigQC))
})


##### Biological Probe QC #####

testData <- setBioProbeQCFlags(testDataRaw, 
                               qcCutoffs=list(minProbeRatio=0.1,
                                              percentFailGrubbs=20))


probeQC <- fData(testData)[["QCFlags"]]
subTest <- subset(testData, subset=Module == gsub(".pkc", "", annotation(testData)[[2]]))
neg_set <- negativeControlSubset(subTest)
probeQCResults <- data.frame(fData(neg_set)[["QCFlags"]], check.names = FALSE)

# Spec 11: test that the number of Low Probe Ratio is correct:------
testthat::test_that("test that the number of Low Probe Ratio is correct", {
  expect_true(sum(probeQC$LowProbeRatio) == sum(fData(testData)["ProbeRatio"] <= 0.1)) 
})



# Spec 12: test that the Local Grubbs Outlier is correct:
neg_set <- assayDataElement(neg_set, elt="preLocalRemoval")
neg_set <- neg_set[, !apply(neg_set, 2, function(x) {max(x) < 10L})]
neg_set <- logtBase(neg_set, base=10L)

# results from outliers::grubbs.test before deprecation
test_results <- readRDS("testData/outliersGrubbsTest.RDS")

test_that("copied grubbs test function works as expected",{
    expect_identical(test_results, apply(neg_set, 2, grubbs.test, two.sided=TRUE))
})

test_list <- lapply(test_results, function(x) {x$p.value < 0.01})
test_outliers <- test_list[test_list == TRUE]
test_nonoutliers <- test_list[test_list == FALSE]
test_outliers <- lapply(test_outliers, function(x) {names(x)})
outlier_flags <- 
  lapply(names(test_outliers), 
         function(x) {probeQCResults[test_outliers[[x]], 
                                     paste0("LocalGrubbsOutlier.", x)]})
testthat::test_that("test that outliers detected are correct", {
  testthat::expect_true(all(unlist(outlier_flags))) 
})
outlier_count <- 
  apply(probeQCResults[, paste0("LocalGrubbsOutlier.", names(test_outliers))],
        2, 
        sum)
testthat::test_that("test that flagged sample/target sets have one outlier", {
  testthat::expect_true(all(outlier_count == 1)) 
})
nonoutlier_count <- 
  apply(probeQCResults[, paste0("LocalGrubbsOutlier.", names(test_nonoutliers))],
        2, 
        sum)
testthat::test_that("test that sample/target sets not flagged are correct", {
  testthat::expect_true(all(nonoutlier_count == 0)) 
})


# Spec 13: test that the Global Grubbs Outlier is correct:
globals_flagged <- rownames(probeQCResults)[which(probeQCResults$GlobalGrubbsOutlier==TRUE)]
probe_flag_percent <- 100 * table(t(as.data.frame(test_outliers))) / dim(testData)[["Samples"]]
global_outliers <- names(probe_flag_percent[probe_flag_percent >= 20])
testthat::test_that("flagged global outliers are correct", {
  testthat::expect_true(all(global_outliers %in% globals_flagged))
  testthat::expect_true(all(globals_flagged %in% global_outliers))
})

# Spec 14: test that genes with less than 3 probes get no grubbs flag:
fewProbes <- names(which(table(fData(testData)$TargetName) < 3))
fewProbes <- fData(testData)$RTS_ID[fData(testData)$TargetName %in% fewProbes] 
grubbsCols <- grep("Grubbs", colnames(probeQC))
testthat::test_that("genes with less than 3 probes get no grubbs flag", {
  testthat::expect_true(all(apply(probeQC[fewProbes,grubbsCols], 2, sum)) == 0)
})

