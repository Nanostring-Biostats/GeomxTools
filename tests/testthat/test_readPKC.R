# Ensure that the package is installed in rstudio.
# Return colnames(pkcFile) comparing to expected
# Return names(metadata(pkcFile)) and compare to expected

library(GeomxTools)
library(testthat)

datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
pkcFile <- readPKCFile(PKCFiles)
lines <- suppressWarnings(c(readLines(file.path(PKCFiles[1])), readLines(file.path(PKCFiles[2]))))


# req 1: test that the column names of PKC files are in correct format:------
testthat::test_that("test that the column names of PKC files are in correct format", {
  expect_true(all(colnames(pkcFile) == c("RTS_ID", "Target", "Module", 
    "CodeClass", "ProbeID", "GeneID", "SystematicName", "Negative")))
})


# req 2: test that the names of metadata of PKC files are in correct format:------
testthat::test_that("test that the column names of PKC files are in correct format", {
  expect_true(all(names(metadata(pkcFile)) == 
    c("PKCFileName", "PKCModule", "PKCFileVersion", "PKCFileDate", 
    "AnalyteType", "MinArea","MinNuclei")))
})

# req3: test that the number of probes is correct:------
testthat::test_that("test that the number of probes is correct", {
  num_probes <- length(grep("\"ProbeID\":", lines))
  expect_true(dim(pkcFile)[1] == num_probes)
})

# Load a simulated second PKC version for the six gene module
verTestPKCFile <- unzip(zipfile = file.path(datadir, "/unittest_pkc.zip"))
multiPKCFiles <- c(PKCFiles[grepl("Six", PKCFiles)], verTestPKCFile)

# req4: check for warning if default PKC provided when no multiple versions:---
testthat::test_that("warning if default PKC given when not needed", {
  expect_warning(readPKCFile(PKCFiles, 
    default_pkc_vers="./Six-gene_test_v1_v1.1.pkc"))
})

# req5: check for error if default PKCs is not a valid pkc file:------
testthat::test_that("check for error if default file name wrong", {
  expect_error(readPKCFile(multiPKCFiles, default_pkc_vers=c("fake pkc name")))
})

# req6: check for error if multiple default PKCs for a module:------
testthat::test_that("check for error if multiple defaults per module", {
  expect_error(readPKCFile(multiPKCFiles, default_pkc_vers=c(
    "./Six-gene_test_v1_v1.1.pkc", "./Six-gene_test_v1_v2.5.pkc")))
})

# req7: check for warning when resolving multiple PKC versions:------
testthat::test_that("check for warning when resolving multiple PKC versions", {
  expect_warning(expect_warning(readPKCFile(multiPKCFiles), 
    "The following probes"), "The following PKC")
})

firstVer <- readPKCFile(multiPKCFiles[1L])
secondVer <- readPKCFile(multiPKCFiles[2L])
combineVer <- suppressWarnings(readPKCFile(multiPKCFiles))
combineReVer <- suppressWarnings(readPKCFile(multiPKCFiles, 
  default_pkc_vers=multiPKCFiles[1L]))
newProbesV1 <- setdiff(firstVer$RTS_ID, secondVer$RTS_ID)
newProbesV2 <- setdiff(secondVer$RTS_ID, firstVer$RTS_ID)
reassignedProbe <- which(combineVer$Target != combineReVer$Target)

# req8: check that only probes in all versions kept:------
testthat::test_that("check module probes in all PKC versions", {
  expect_equal(nrow(combineVer) * 2, (nrow(firstVer) + nrow(secondVer) - 
    length(newProbesV1) - length(newProbesV2)))
  expect_true(sum(c(newProbesV1, newProbesV2) %in% combineVer$RTS_ID) == 0)
  expect_equal(nrow(combineVer), nrow(combineReVer))
})

# req9: check that default PKC target assignments are used for probes:----
testthat::test_that("check module probes in all PKC versions", {
  expect_false(all(combineVer$Target == combineReVer$Target))
  expect_true(combineVer$RTS_ID[reassignedProbe] == 
    combineReVer$RTS_ID[reassignedProbe])
  expect_true(combineVer$Target[reassignedProbe] == 
    secondVer[combineVer$RTS_ID[reassignedProbe] == secondVer$RTS_ID, "Target"])
  expect_true(combineReVer$Target[reassignedProbe] == 
    firstVer[combineReVer$RTS_ID[reassignedProbe] == firstVer$RTS_ID, "Target"])
  expect_false(combineReVer$Module[reassignedProbe] %in% combineVer$Module)
  expect_false(combineVer$Module[reassignedProbe] %in% combineReVer$Module)
})
