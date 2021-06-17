### CONFIGURATION SECTION ###
# load data
demoData <- readRDS(file= system.file("extdata","DSP_NGS_Example_Data", "demoData.rds", package = "GeomxTools"))
demoData <- normalize(demoData , data_type="RNA", norm_method="quant",
                        desiredQuantile = .9, toElt = "q_norm")  ## test with user parameter
demoData <- normalize(demoData, norm_method="quant") ## test defaults quantile = .75
housekeepers <- c('C1orf43','GPI','OAZ1','POLR2A','PSMB2','RAB7A',
                               'SDHA','SNRPD3','TBC1D10B','TPM4','TUBB','UBB')
#run aggregateCounts function on the data
target_demoData <- aggregateCounts(demoData)
tolerance <- 0.000001 # accepted tolerance between normalized values
#############################

########### Quantile Normalization test
#### req 1 check that normfactors are in in pData of demoData
test_that("quantile norm factors are present", {
  expect_snapshot(demoData@phenoData@data[["q_norm_qFactors"]])
  expect_snapshot(demoData@phenoData@data[["normFactors"]])
})

#### req 2 verify calculation of q90 norm factors
# Compute normalization factors from demoData
# compute 90% quantile count for samples and divide by geomean
thresh <- assayDataApply(demoData, 2, quantile, probs = .9)
expectedOutputData <- thresh/geoMean(thresh)
expectedOutputData <- unname(expectedOutputData)

# Extract normalized data from object
actualOutputData <- sData(demoData)[,"q_norm_qFactors"]
test_that("quantile norm factors are correct", {
  expect_equal(expectedOutputData, actualOutputData)
})

#### req 2b verify calculation of q75 norm factors
# compute 75% quantile norm factors for samples and divide by geomean
thresh <- assayDataApply(demoData, 2, quantile, probs = .75)
expectedOutputData <- thresh/geoMean(thresh)
expectedOutputData <- unname(expectedOutputData)

# Extract normalized data from object
actualOutputData <- sData(demoData)[,"normFactors"]
test_that("quantile norm factors for 75% are correct", {
  expect_equal(expectedOutputData, actualOutputData)
})


#### req 3 verify calculation of quantile norm values
# compute normalize count for samples and divide by geomean
thresh <- assayDataApply(demoData, 2, quantile, probs = .75)
norm_quant <- function(x){
  x <- x/(thresh/geoMean(thresh))
}
expectedOutputData <- t(assayDataApply(demoData, 1, norm_quant))

# Extract normalized data from object
actualOutputData <- assayData(demoData)[["exprs_norm"]]
test_that("quantile norm values are correct", {
  expect_equal(expectedOutputData, actualOutputData)
})

########### Negative Normalization test
test_that("aggregateCounts", {
  expect_snapshot(featureData(target_demoData))
  expect_snapshot(exprs(target_demoData))
})

#### req 4 verify calculation of negative norm factors
#call negative normalization
target_demoData <- normalize(target_demoData , data_type="RNA", norm_method="neg",
                      toElt = "neg_norm")
# compute negative norm factors for samples and divide by geomean
negfactors <- assayDataApply(negativeControlSubset(target_demoData), 2, ngeoMean)
norm_neg <- function(x){
  x <- x/(negfactors/geoMean(negfactors))
}
expectedOutputData <- t(assayDataApply(target_demoData, 1, norm_neg))

# Extract normalized data from object
actualOutputData <- assayData(target_demoData)[["neg_norm"]]
test_that("negative norm values are correct", {
  expect_equal(expectedOutputData, actualOutputData)
})

########### Housekeeping Norm test
#### req 5 verify calculation of housekeeping norm factors
#call hk normalization
target_demoData <- normalize(target_demoData , data_type="RNA", norm_method="hk",
                             fromElt="exprs", toElt="hk_norm")
# compute hk norm factors for samples and divide by geomean
hksubset <- target_demoData[which(featureData(target_demoData)[["TargetName"]] %in% housekeepers),]
hkfactors <- assayDataApply(hksubset, 2, ngeoMean)
norm_hk <- function(x){
  x <- x/(hkfactors/geoMean(hkfactors))
}
expectedOutputData <- t(assayDataApply(target_demoData, 1, norm_hk))

# Extract normalized data from object
actualOutputData <- assayData(target_demoData)[["hk_norm"]]
test_that("hk norm values are correct", {
  expect_equal(expectedOutputData, actualOutputData)
})

# ########### Subtract Background Norm test
# #### req 6 verify calculation of subtract background norm factors
# #call subtract bg normalization
# target_demoData <- normalize(target_demoData , data_type="RNA",
#                              norm_method="subtractBackground", fromElt="exprs", toElt="bg_norm")
#
# # compute neg norm factors for samples and divide by geomean
# negfactors <- assayDataApply(negativeControlSubset(target_demoData), 2, ngeoMean)
# norm_bg <- function(x){
#   x <- (x- negfactors)
# }
# expectedOutputData <- t(assayDataApply(target_demoData, 1, norm_bg))
# head(exprs(target_demoData))
# head(negfactors)
# 3.776350                   -
# # Extract normalized data from object
# actualOutputData <- assayData(target_demoData)[["bg_norm"]]
# test_that("bg subtract norm values are correct", {
#   expect_equal(expectedOutputData, actualOutputData)
# })





