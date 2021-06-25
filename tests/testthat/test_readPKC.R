# Ensure that the package is installed in rstudio.
# Return colnames(pkcFile) comparing to expected
# Return names(metadata(pkcFile)) and compare to expected

library(GeomxTools)
library(testthat)

datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
pkcFile <- readPKCFile(PKCFiles)


# req 1: test that the column names of PKC files are in correct format:------
testthat::test_that("test that the column names of PKC files are in correct format", {
  expect_true(all(colnames(pkcFile) == c("RTS_ID", "Target", "Module", "CodeClass", "ProbeID", "Negative")))
})
# This is different from what's in Quick Base: all(colnames(pkcFile) == c("RTS_ID", "Gene", "Module", "Negative")) 


# req 2: test that the names of metadata of PKC files are in correct format:------
testthat::test_that("test that the column names of PKC files are in correct format", {
  expect_true(all(names(metadata(pkcFile)) == c("PKCFileName", "PKCFileVersion", "PKCFileDate", "AnalyteType", "MinArea","MinNuclei")))
})



