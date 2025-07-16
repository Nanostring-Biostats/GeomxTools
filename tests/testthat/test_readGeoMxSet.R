# Return colnames and rownames of testData@assayData$expr comparing to expected
# Return colnames and rownames of testData@phenoData@data comparing to expected
# Return colnames and rownames of testData@protocolData@data comparing to expected
# Return genes of testData@featureData@data comparing to expected
# Return testData@experimentData comparing to expected

library(GeomxTools)
library(testthat)
library(stringr)

datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)[1:10]
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
SampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

testData <- suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
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

pkcFile <- readPKCFile(PKCFiles)

DCCFiles <- DCCFiles[!basename(DCCFiles) %in% unique(sData(testData)$NTC_ID)]


# Spec 1: test that the column names and the rownames of testData@assayData$exprs match those in DCC files and PKC Files respectively:------
testthat::test_that("test that the column names and the rownames of testData@assayData$exprs match those in DCC files and PKC Files respectively", {
  expect_true(all(basename(DCCFiles) %in% colnames(testData@assayData$exprs)))
  expect_true(all(unique(pkcFile$RTS_ID) %in% rownames(testData@assayData$exprs)))
})




# Spec 2: test that the column names and the rownames of testData@phenoData$data match those in DCC files and PKC Files respectively:------ 
testthat::test_that("test that the column names and the rownames of testData@phenoData$exprs match those in DCC files and PKC Files respectively", {
  phenoDataDccColName <- "Sample_ID"
  protocolDataColNames <- c("aoi",
                            "cell_line",
                            "roi_rep",
                            "pool_rep",
                            "slide_rep")
  experimentDataColNames <- c("panel")
  pheno_tab <- readxl::read_xlsx(SampleAnnotationFile, sheet = 'CW005')
  colnames(pheno_tab) <- str_replace_all(colnames(pheno_tab),'\\.',' ')
  expect_true(all(basename(DCCFiles) %in% rownames(testData@phenoData@data)))
  expect_true(all(colnames(pheno_tab) %in% c(names(testData@phenoData@data), # what is pheno_tab?
                                             phenoDataDccColName,
                                             protocolDataColNames,
                                             experimentDataColNames)))
})




# Spec 3: test that the column names and the rownames of testData@protocolData$data match those in DCC files and PKC Files respectively:------
testthat::test_that("test that the column names and the rownames of testData@protocolData$exprs match those in DCC files and PKC Files respectively", {
  protocolDataColNames <- c("aoi",
                            "cell_line",
                            "roi_rep",
                            "pool_rep",
                            "slide_rep")
  expect_true(all(basename(DCCFiles) %in% rownames(testData@protocolData@data)))
  expect_true(all(protocolDataColNames %in% names(testData@protocolData@data)))
})




# Spec 4: test that the genes in testData@featureData@data match those in PKC Files:------
testthat::test_that("test that the genes in testData@featureData$exprs match those in PKC Files", {
  expect_true(dim(testData@featureData@data)[1] == length(unique(pkcFile$RTS_ID))) # QuickBase: length(unique(pkcFile$Gene))
  expect_true(all(unique(pkcFile$RTS_ID) %in% testData@featureData@data$RTS_ID))
})



# Spec 5: test that the names in testData@experimentData@other are in correct format:------
testthat::test_that("test that the names in testData@experimentData@other are in correct format", {
  experimentDataColNames <- c("panel")
  experimentDataColNames <- c(experimentDataColNames, 
                              "PKCFileName",
                              "PKCFileVersion",
                              "PKCFileDate",
                              "AnalyteType",
                              "MinArea",
                              "MinNuclei")
  expect_true(all(experimentDataColNames %in% names(testData@experimentData@other))) 
})


#Spec 6: test that the counts in testData@assayData$exprs match those in DCC files
testthat::test_that("test that the counts of testData@assayData$exprs match those in DCC files", {
  correct <- TRUE
  i <- 1
  while(correct == TRUE & i < length(DCCFiles)){
    dccFile <- readDccFile(DCCFiles[i])
    
    mtxCount <- testData@assayData$exprs[,basename(DCCFiles[i])]
    genes <- match(dccFile$Code_Summary$RTS_ID, names(mtxCount))
    
    correct <- all(mtxCount[genes] == dccFile$Code_Summary$Count) & all(mtxCount[!genes] == 0)
    
    i <- i+1
  }
  expect_true(correct)
})

n_aois <- 6

testData <-
  suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
                                          pkcFiles = PKCFiles,
                                          phenoDataFile = SampleAnnotationFile,
                                          phenoDataSheet = "CW005",
                                          phenoDataDccColName = "Sample_ID",
                                          protocolDataColNames = c("aoi",
                                                                   "cell_line",
                                                                   "roi_rep",
                                                                   "pool_rep",
                                                                   "slide_rep"),
                                          experimentDataColNames = c("panel"),
                                          col_types = c(rep("text", 4), "numeric", rep("text", 2),
                                                        "numeric", rep("text", 4)),
                                          n_max = n_aois + 1))

# Spec 7: test that ... gets translated to read_xlsx():------
testthat::test_that("test that ... gets translated to read_xlsx()", {
  expect_true(dim(testData)[2] == n_aois) 
  expect_true(class(pData(testData)$roi) == "numeric")
  expect_true(class(pData(testData)$area) == "numeric")
})

# Spec 8: test that only valid ... gets translated to read_xlsx():------
testthat::test_that("test that only valid ... gets translated to read_xlsx()", {
  expect_error(testData <-
                 suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
                                                         pkcFiles = PKCFiles,
                                                         phenoDataFile = SampleAnnotationFile,
                                                         phenoDataSheet = "CW005",
                                                         phenoDataDccColName = "Sample_ID",
                                                         protocolDataColNames = c("aoi",
                                                                                  "cell_line",
                                                                                  "roi_rep",
                                                                                  "pool_rep",
                                                                                  "slide_rep"),
                                                         experimentDataColNames = c("panel"),
                                                         trim_ws = c(rep("text", 4), "numeric", rep("text", 2),
                                                                       "numeric", rep("text", 4)))))
})

proteinDatadir <- system.file("extdata", "DSP_Proteogenomics_Example_Data",
                       package="GeomxTools")
proteinPKCs <- unzip(zipfile = file.path(proteinDatadir,  "/pkcs.zip"))

# Spec 10. PKCs are removed if they aren't in the provided config file.
testthat::test_that("PKCs are removed if they aren't in the config", {
  testDataWithExtraConfig <- 
    suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                            pkcFiles = c(PKCFiles, proteinPKCs[1]),
                                            phenoDataFile = SampleAnnotationFile,
                                            phenoDataSheet = "CW005",
                                            phenoDataDccColName = "Sample_ID",
                                            protocolDataColNames = c("aoi",
                                                                     "cell_line",
                                                                     "roi_rep",
                                                                     "pool_rep",
                                                                     "slide_rep"),
                                            experimentDataColNames = c("panel"),
                                            configFile = paste0(datadir, "/fakeTesting_config.ini")))
  
  expect_identical(testData@annotation, testDataWithExtraConfig@annotation)
  expect_identical(rownames(testData), rownames(testDataWithExtraConfig))
  
  # TargetGroup is not expected to be identical with an extra pkc
  targGroupCol <- which(colnames(fData(testData)) == "TargetGroup")
  expect_identical(fData(testData)[,-targGroupCol], fData(testDataWithExtraConfig)[,-targGroupCol])
})

DCCFiles <- unzip(zipfile = file.path(proteinDatadir,  "/DCCs.zip"))[1:10]
annots <- file.path(proteinDatadir, "Annotation.xlsx")

RNAData <- suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                                   pkcFiles = proteinPKCs,
                                                   phenoDataFile = annots,
                                                   phenoDataSheet = "Annotations",
                                                   phenoDataDccColName = "Sample_ID",
                                                   protocolDataColNames = c("Tissue", 
                                                                            "Segment_Type", 
                                                                            "ROI.Size")))


ProteinData <- suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                                       pkcFiles = proteinPKCs,
                                                       phenoDataFile = annots,
                                                       phenoDataSheet = "Annotations",
                                                       phenoDataDccColName = "Sample_ID",
                                                       protocolDataColNames = c("Tissue", 
                                                                                "Segment_Type", 
                                                                                "ROI.Size"),
                                                       analyte = "protein"))

# Spec 11: Only a single analyte is read in to a GeomMxSet object
testthat::test_that("only a single analyte is read in to a GeoMxSet object",{
  expect_false(dim(RNAData)["Features"] == dim(ProteinData)["Features"])
  expect_true(dim(RNAData)["Samples"] == dim(ProteinData)["Samples"])
  expect_false(analyte(RNAData) == analyte(ProteinData))
  expect_true(analyte(RNAData) == "RNA")
  expect_true(analyte(ProteinData) == "Protein")
  expect_true(featureType(RNAData) == "Probe")
  expect_true(featureType(ProteinData) == "Target")
  expect_false(all(annotation(RNAData) == annotation(ProteinData)))
  expect_false(any(rownames(RNAData) %in% rownames(ProteinData)))
  expect_true(all(colnames(RNAData) %in% colnames(ProteinData)))
  expect_true(all(sData(RNAData) == sData(ProteinData)))
})

# Spec 12: Object should be able to be loaded with Lab worksheet or data.frame
phenoData <- readxl::read_xlsx(SampleAnnotationFile, col_names = TRUE, sheet = "CW005")
sampleLW <- file.path(datadir, "GxT_testing_labWorksheet.txt")

DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)[1:10]
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
SampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

testData <- suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
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

dfData <- suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
                                                  pkcFiles = PKCFiles,
                                                  phenoData = phenoData,
                                                  phenoDataFile = SampleAnnotationFile,
                                                  phenoDataSheet = "CW005",
                                                  phenoDataDccColName = "Sample_ID",
                                                  protocolDataColNames = c("aoi",
                                                                           "cell_line",
                                                                           "roi_rep",
                                                                           "pool_rep",
                                                                           "slide_rep"),
                                                  experimentDataColNames = c("panel")))

lwData <- suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
                                                  pkcFiles = PKCFiles,
                                                  phenoDataFile = sampleLW,
                                                  phenoDataSheet = "CW005",
                                                  phenoDataDccColName = "Sample_ID",
                                                  protocolDataColNames = c("Aoi",
                                                                           "Cell_line",
                                                                           "Roi_rep",
                                                                           "Pool_rep",
                                                                           "Slide_rep"),
                                                  experimentDataColNames = c("Panel")))

testthat::test_that("GeoMxSet objects can be read in using lab worksheet and data.frames",{
  testthat::expect_equal(object = dfData, expected = testData)
  
  # change column names to be all lower case
  colnames(lwData@phenoData) <- tolower(colnames(lwData@phenoData))
  colnames(lwData@protocolData) <- tolower(colnames(lwData@protocolData))
  names(lwData@experimentData@other) <- tolower(names(lwData@experimentData@other))
  colnames(testData@phenoData) <- tolower(colnames(testData@phenoData))
  colnames(testData@protocolData) <- tolower(colnames(testData@protocolData))
  names(testData@experimentData@other) <- tolower(names(testData@experimentData@other))
  
  testthat::expect_equal(object = lwData, expected = testData)
})

testthat::test_that("Loading will error if expected columns are not available in phenoData",{
  testthat::expect_warning(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
                                                  pkcFiles = PKCFiles,
                                                  phenoData = phenoData[,-2], # remove "slide name" column
                                                  phenoDataFile = SampleAnnotationFile,
                                                  phenoDataSheet = "CW005",
                                                  phenoDataDccColName = "Sample_ID",
                                                  protocolDataColNames = c("aoi",
                                                                           "cell_line",
                                                                           "roi_rep",
                                                                           "pool_rep",
                                                                           "slide_rep"),
                                                  experimentDataColNames = c("panel")), 
                           regexp = "'slide name' is an expected column in phenoData.")
  
  testthat::expect_error(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
                                                pkcFiles = PKCFiles,
                                                phenoData = phenoData[,c(1,1,2)],
                                                phenoDataFile = SampleAnnotationFile,
                                                phenoDataSheet = "CW005",
                                                phenoDataDccColName = "Sample_ID",
                                                protocolDataColNames = c("aoi",
                                                                         "cell_line",
                                                                         "roi_rep",
                                                                         "pool_rep",
                                                                         "slide_rep"),
                                                experimentDataColNames = c("panel")),
                         regexp = "column names are duplicated in given")
  
  testthat::expect_error(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
                                                pkcFiles = PKCFiles,
                                                phenoData = phenoData,
                                                phenoDataFile = SampleAnnotationFile,
                                                phenoDataSheet = "CW005",
                                                phenoDataDccColName = "Sample_ID",
                                                protocolDataColNames = c("aoi",
                                                                         "cell_line",
                                                                         "roi_rep",
                                                                         "pool_rep",
                                                                         "slide_rep",
                                                                         "fake_column"), # fake column name
                                                experimentDataColNames = c("panel")),
                         regexp = "Columns specified in `protocolDataColNames` are not found in")
  
  testthat::expect_error(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
                                                pkcFiles = PKCFiles,
                                                phenoData = phenoData,
                                                phenoDataFile = SampleAnnotationFile,
                                                phenoDataSheet = "CW005",
                                                phenoDataDccColName = "Sample_ID",
                                                protocolDataColNames = c("aoi",
                                                                         "cell_line",
                                                                         "roi_rep",
                                                                         "pool_rep",
                                                                         "slide_rep"),
                                                experimentDataColNames = c("panel",
                                                                           "fake_column")), # fake column name
                         regexp = "Columns specified in `experimentDataColNames` are not found in")
  
  badPhenoData <- phenoData
  badPhenoData$`slide name`[1] <- "slide1"
  
  testthat::expect_warning(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
                                                pkcFiles = PKCFiles,
                                                phenoData = badPhenoData,
                                                phenoDataFile = SampleAnnotationFile,
                                                phenoDataSheet = "CW005",
                                                phenoDataDccColName = "Sample_ID",
                                                protocolDataColNames = c("aoi",
                                                                         "cell_line",
                                                                         "roi_rep",
                                                                         "pool_rep",
                                                                         "slide_rep"),
                                                experimentDataColNames = c("panel")), 
                         regexp = "No NTCs were found. These are determined by")
})

