library(GeomxTools)
library(testthat)

datadir <- system.file("extdata","DSP_Proteogenomics_Example_Data",
                       package = "GeomxTools")

DCCFiles <- unzip(zipfile = file.path(datadir,  "/DCCs.zip"))
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
SampleAnnotationFile <- file.path(datadir, "Annotation.xlsx")


RNAData <- suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                                   pkcFiles = PKCFiles,
                                                   phenoDataFile = SampleAnnotationFile,
                                                   phenoDataSheet = "Annotations",
                                                   phenoDataDccColName = "Sample_ID",
                                                   protocolDataColNames = c("Tissue", 
                                                                            "Segment_Type", 
                                                                            "ROI.Size"),
                                                   configFile = NULL,
                                                   analyte = "RNA",
                                                   phenoDataColPrefix = "",
                                                   experimentDataColNames = NULL))

proteinData <- suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                                       pkcFiles = PKCFiles,
                                                       phenoDataFile = SampleAnnotationFile,
                                                       phenoDataSheet = "Annotations",
                                                       phenoDataDccColName = "Sample_ID",
                                                       protocolDataColNames = c("Tissue", 
                                                                                "Segment_Type", 
                                                                                "ROI.Size"),
                                                       configFile = NULL,
                                                       analyte = "protein",
                                                       phenoDataColPrefix = "",
                                                       experimentDataColNames = NULL))

# Spec 11: Only a single analyte is read in to a GeomMxSet object.:------
testthat::test_that("Datasets can be subset by given analyte", {
  expect_true(analyte(proteinData) == "Protein")
  expect_true(analyte(proteinData) != "RNA")
  expect_true(analyte(RNAData) == "RNA")
  expect_true(analyte(RNAData) != "Protein")
  
  expect_true(featureType(RNAData) == "Probe")
  expect_true(featureType(proteinData) == "Target")
  
  expect_error(suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                      pkcFiles = PKCFiles,
                                      phenoDataFile = SampleAnnotationFile,
                                      phenoDataSheet = "Annotations",
                                      phenoDataDccColName = "Sample_ID",
                                      protocolDataColNames = c("Tissue", 
                                                               "Segment_Type", 
                                                               "ROI.Size"),
                                      configFile = NULL,
                                      analyte = "TEST",
                                      phenoDataColPrefix = "",
                                      experimentDataColNames = NULL)))
  
  expect_true(all(proteinData@annotation %in% paste0(unique(fData(proteinData)$Module), ".pkc")))
  expect_true(all(RNAData@annotation %in% paste0(unique(fData(RNAData)$Module), ".pkc")))
})

RNAData <- aggregateCounts(RNAData)

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

RNAData_seqQC <- setSeqQCFlags(RNAData, 
                                   qcCutoffs=list(minSegmentReads=1000, 
                                                  percentAligned=80, 
                                                  percentSaturation=50))

RNAData_background <- setBackgroundQCFlags(RNAData_seqQC, 
                                               qcCutoffs=list(minNegativeCount=10, 
                                                              maxNTCCount=60))

RNAData <- setSegmentQCFlags(RNAData, 
                                 qcCutoffs=list(minNuclei=20,
                                                minArea=16000,
                                                minNegativeCount=10,
                                                maxNTCCount=60,
                                                minSegmentReads=1000, 
                                                percentAligned=80, 
                                                percentSaturation=50))

# Spec 20: The appropriate segment QC flags shall be added to a protein dataset.:------
test_that("Correct segment flags are added to protein data",{
  expect_false("LowNegatives" %in% colnames(protocolData(proteinData)$QCFlags))
  expect_false("LowNegatives" %in% colnames(protocolData(proteinData_seqQC)$QCFlags))
  expect_false("LowNegatives" %in% colnames(protocolData(proteinData_background)$QCFlags))
  expect_equal(ncol(protocolData(proteinData_seqQC)$QCFlags), ncol(protocolData(proteinData_background)$QCFlags))
  expect_equal(ncol(protocolData(proteinData_seqQC)$QCFlags), ncol(protocolData(proteinData)$QCFlags))
  expect_gt(ncol(protocolData(RNAData_background)$QCFlags), ncol(protocolData(RNAData_seqQC)$QCFlags))
  expect_equal(ncol(protocolData(RNAData_background)$QCFlags), ncol(protocolData(RNAData)$QCFlags))
})

# Spec 21: A warning is given if protein data is run through setBioProbeQCFlags.:------
test_that("Warnings are given if protein data is run through setBioProbeQCFlags", {
  expect_warning(expect_warning(expect_warning(setBioProbeQCFlags(proteinData, 
                                    qcCutoffs=list(minProbeRatio=0.1,
                                                   percentFailGrubbs=20)))))
})


igg.names <- iggNames(proteinData)
hk.names <- hkNames(proteinData)

HOUSEKEEPERS <- c(
  "C1orf43", "GPI", "OAZ1", "POLR2A", "PSMB2", "RAB7A",
  "SDHA", "SNRPD3", "TBC1D10B", "TPM4", "TUBB", "UBB"
)

# Spec 1: igg.names shall return the expected target names.:------
test_that("Expected IgGs are returned",{
  expect_true(all(grepl(pattern = "IgG", igg.names)))
  expect_equal(length(igg.names), length(grep(pattern = "IgG", fData(proteinData)$TargetName)))
  expect_true(all(igg.names == fData(proteinData)$TargetName[fData(proteinData)$CodeClass == "Negative"]))
  
  expect_warning(iggNames(RNAData))
})

# Spec 2: hk.names shall return the expected target names.:------
test_that("Expected HK are returned",{
  expect_true(all(hk.names == fData(proteinData)$TargetName[fData(proteinData)$CodeClass == "Control"]))
  
  expect_warning(hkNames(RNAData))
})


test_that("Concordance plots are plotted", {
  expect_error(plotConcordance(igg.names, proteinData, "Segment_Typ"))
  expect_error(fig <- plotConcordance(igg.names, proteinData, "Segment_Type"), NA)
  expect_true(class(fig)[1] == "gg")
  expect_error(fig, NA)
  
  expect_error(fig <- plotConcordance(igg.names, proteinData, "Tissue"), NA)
  expect_true(class(fig)[1] == "gg")
  expect_error(fig, NA)
  
  expect_error(fig <- plotConcordance(HOUSEKEEPERS[1:4], RNAData, "Segment_Type"), NA)
  expect_true(class(fig)[1] == "gg")
  expect_error(fig, NA)
})

proteinData <- normalize(proteinData, norm_method = "hk", toElt="hk_norm")
proteinData <- normalize(proteinData, norm_method = "hk", toElt="hk_norm_givenHK", housekeepers = hk.names)
proteinData <- normalize(proteinData, norm_method = "neg", toElt="neg_norm")
proteinData <- normalize(proteinData, norm_method = "quant", toElt="q3_norm")
proteinData <- normalize(proteinData, norm_method = "subtractBackground", toElt="bgSub_norm")

# Spec 13: All normalization methods work on protein data.:------
test_that("Protein data is normalized",{
  expect_true(all(proteinData@assayData$exprs != proteinData@assayData$hk_norm))
  expect_true(all(proteinData@assayData$exprs != proteinData@assayData$neg_norm))
  expect_true(all(proteinData@assayData$exprs != proteinData@assayData$q3_norm))
  expect_true(all(proteinData@assayData$exprs != proteinData@assayData$bgSub_norm))
  
  expect_equal(proteinData@assayData$hk_norm, proteinData@assayData$hk_norm_givenHK)
  
  expect_true(all(c("hk_norm_hkFactors", "neg_norm_negFactors", "q3_norm_qFactors") %in% 
                    colnames(pData(proteinData))))
})


normfactors <- computeNormalizationFactors(object = proteinData)

normfactors_area <- computeNormalizationFactors(object = proteinData,
                                                  area = "AOI.Size.um2")

normfactors_area_nuc <- computeNormalizationFactors(object = proteinData,
                                                      area = "AOI.Size.um2",
                                                      nuclei = "Nuclei.Counts")

# Spec 3: computeNormalizationFactors and normalize calculations match.:------
test_that("computeNormalizationFactors vs normalization",{
  expect_true(all(pData(proteinData)$hk_norm_hkFactors == 
                    as.data.frame(normfactors)$`HK geomean`/exp(mean(log(as.data.frame(normfactors)$`HK geomean`)))))
  expect_true(all(pData(proteinData)$neg_norm_negFactors[,1] == 
                    as.data.frame(normfactors)$`Neg geomean`/exp(mean(log(as.data.frame(normfactors)$`Neg geomean`)))))

  expect_equal(as.numeric(dim(normfactors)), as.numeric(c(ncol(proteinData), 2)))
  expect_equal(as.numeric(dim(normfactors_area)), as.numeric(c(ncol(proteinData), 3)))
  expect_equal(as.numeric(dim(normfactors_area_nuc)), as.numeric(c(ncol(proteinData), 4)))
  
  expect_error(compute_normalization_factors(object = RNAData))
})



test_that("QC plots are plotted", {
  expect_error(fig <- plotNormFactorConcordance(object = proteinData, plotFactor = "Segment_Type", 
                                         normfactors = normfactors), NA)
  expect_true(class(fig)[1] == "gg")
  expect_error(fig, NA)
  
  expect_error(fig <- qcProteinSignal(object = proteinData,
                                 neg.names = igg.names), NA)
  expect_true(class(fig) == "function")
  expect_error(fig(), NA)
  
  proteinOrder <- qcProteinSignalNames(object = proteinData,
                                       neg.names = igg.names)
  
  expect_true(all(proteinOrder %in% rownames(proteinData)))
  expect_true(all(rownames(proteinData) %in% proteinOrder))
  
  expect_error(qcProteinSignal(object = RNAData))
})


