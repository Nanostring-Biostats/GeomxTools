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

# req4: check for error if default PKCs is not a valid pkc file:------
testthat::test_that("test that the number of probes is correct", {
  expect_error(readPKCFile(PKCFiles, default_pkc_vers=c("fake pkc name")))
})


