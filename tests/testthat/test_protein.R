library(GeomxTools)
library(testthat)

testData <- readRDS(file= system.file("extdata","DSP_Proteogenomics_Example_Data", 
                                      "proteinData.rds", package = "GeomxTools"))
aggTestData <- aggregateCounts(testData)

proteinData <- analyteSubset(object = aggTestData, analyte = "protein")
RNAData <- analyteSubset(object = aggTestData, analyte = "RNA")

testthat::test_that("Datasets can be subset by given analyte", {
  expect_true(all(fData(proteinData)$AnalyteType == "Protein"))
  expect_true(all(fData(proteinData)$AnalyteType != "RNA"))
  expect_true(all(fData(RNAData)$AnalyteType == "RNA"))
  expect_true(all(fData(RNAData)$AnalyteType != "Protein"))
  
  expect_error(analyteSubset(aggTestData, "TEST"))
  
  expect_true(all(proteinData@annotation %in% paste0(unique(fData(proteinData)$Module), ".pkc")))
  expect_true(all(RNAData@annotation %in% paste0(unique(fData(RNAData)$Module), ".pkc")))
})


proteinData_seqQC <- setSegmentQCFlags(proteinData, 
                             qcCutoffs=list(minSegmentReads=1000, 
                                            percentAligned=80, 
                                            percentSaturation=50))

proteinData_background <- setBackgroundQCFlags(proteinData_seqQC, 
                                    qcCutoffs=list(minNegativeCount=10, 
                                                   maxNTCCount=60))

proteinData <- setSegmentQCFlags(proteinData, 
                                 qcCutoffs=list(minNuclei=20,
                                                minArea=16000,
                                                minNegativeCount=10,
                                                maxNTCCount=60,
                                                minSegmentReads=1000, 
                                                percentAligned=80, 
                                                percentSaturation=50))

aggTestData_seqQC <- setSeqQCFlags(aggTestData, 
                                   qcCutoffs=list(minSegmentReads=1000, 
                                                  percentAligned=80, 
                                                  percentSaturation=50))

aggTestData_background <- setBackgroundQCFlags(aggTestData_seqQC, 
                                               qcCutoffs=list(minNegativeCount=10, 
                                                              maxNTCCount=60))

aggTestData <- setSegmentQCFlags(aggTestData, 
                                 qcCutoffs=list(minNuclei=20,
                                                minArea=16000,
                                                minNegativeCount=10,
                                                maxNTCCount=60,
                                                minSegmentReads=1000, 
                                                percentAligned=80, 
                                                percentSaturation=50))

test_that("Correct segment flags are added to protein data",{
  expect_equal(ncol(protocolData(proteinData_seqQC)$QCFlags), ncol(protocolData(proteinData_background)$QCFlags))
  expect_equal(ncol(protocolData(proteinData_seqQC)$QCFlags), ncol(protocolData(proteinData)$QCFlags))
  expect_gt(ncol(protocolData(aggTestData_background)$QCFlags), ncol(protocolData(aggTestData_seqQC)$QCFlags))
  expect_equal(ncol(protocolData(aggTestData_background)$QCFlags), ncol(protocolData(aggTestData)$QCFlags))
})


test_that("Warnings are given if protein data is run through setBioProbeQCFlags", {
  expect_warning(expect_warning(expect_warning(setBioProbeQCFlags(proteinData, 
                                    qcCutoffs=list(minProbeRatio=0.1,
                                                   percentFailGrubbs=20)))))
})


igg.names <- igg_names(proteinData)
hk.names <- hk_names(proteinData)

test_that("Expected IgGs are returned",{
  expect_true(all(grepl(pattern = "IgG", igg.names)))
  expect_equal(length(igg.names), length(grep(pattern = "IgG", fData(proteinData)$TargetName)))
  expect_true(all(igg.names == fData(proteinData)$TargetName[fData(proteinData)$CodeClass == "Negative"]))
  
  expect_warning(length(igg_names(RNAData)) == 0)
  expect_equal(igg_names(aggTestData), igg.names)
})

test_that("Expected HK are returned",{
  expect_true(all(hk.names == fData(proteinData)$TargetName[fData(proteinData)$CodeClass == "Control"]))
  
  expect_warning(length(hk_names(RNAData)) == 0)
  expect_equal(hk_names(aggTestData), hk.names)
})

test_that("Concordance plots are plotted", {
  expect_error(plot_concordance(igg.names, proteinData, "Segment_Typ"))
  expect_error(plot_concordance(igg.names, proteinData, "Segment_Type"), NA)
  expect_error(plot_concordance(igg.names, proteinData, "Tissue"), NA)
})

proteinData <- normalize(proteinData, norm_method = "hk", toElt="hk_norm")
proteinData <- normalize(proteinData, norm_method = "hk", toElt="hk_norm_givenHK", housekeepers = hk.names)
proteinData <- normalize(proteinData, norm_method = "neg", toElt="neg_norm")
proteinData <- normalize(proteinData, norm_method = "quant", toElt="q3_norm")

test_that("Protein data is normalized",{
  expect_true(all(proteinData@assayData$exprs != proteinData@assayData$hk_norm))
  expect_true(all(proteinData@assayData$exprs != proteinData@assayData$neg_norm))
  expect_true(all(proteinData@assayData$exprs != proteinData@assayData$q3_norm))
  
  expect_equal(proteinData@assayData$hk_norm, proteinData@assayData$hk_norm_givenHK)
  
  expect_true(all(c("hk_norm_hkFactors", "neg_norm_negFactors", "q3_norm_qFactors") %in% 
                    colnames(pData(proteinData))))
  
  expect_error(normalize(aggTestData, norm_method = "neg"))
  expect_error(normalize(aggTestData, norm_method = "hk"))
})


normfactors <- compute_normalization_factors(object = proteinData)

normfactors_area <- compute_normalization_factors(object = proteinData,
                                                  area = "AOI.Size.um2")

normfactors_area_nuc <- compute_normalization_factors(object = proteinData,
                                                      area = "AOI.Size.um2",
                                                      nuclei = "Nuclei.Counts")

test_that("compute_normalization_factors vs normalization",{
  expect_true(all(pData(proteinData)$hk_norm_hkFactors == 
                    as.data.frame(normfactors)$`HK geomean`/exp(mean(log(as.data.frame(normfactors)$`HK geomean`)))))
  expect_true(all(pData(proteinData)$neg_norm_negFactors[,1] == 
                    as.data.frame(normfactors)$`Neg geomean`/exp(mean(log(as.data.frame(normfactors)$`Neg geomean`)))))

  expect_equal(as.numeric(dim(normfactors)), as.numeric(c(ncol(proteinData), 2)))
  expect_equal(as.numeric(dim(normfactors_area)), as.numeric(c(ncol(proteinData), 3)))
  expect_equal(as.numeric(dim(normfactors_area_nuc)), as.numeric(c(ncol(proteinData), 4)))
})

test_that("QC plots are plotted", {
  expect_error(plot_normFactor_concordance(object = proteinData, plot_factors = "Segment_Type", 
                                           normfactors = normfactors), NA)
  expect_error(qc_protein_signal(object = proteinData,
                                 neg.names = igg.names), NA)
})


