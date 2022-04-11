### CONFIGURATION SECTION ###
library(GeomxTools)
library(testthat)
library(EnvStats)

# load data
demoData <- readRDS(file= system.file("extdata","DSP_NGS_Example_Data", "demoData.rds", package = "GeomxTools"))
demoData <- shiftCountsOne(demoData)
demoData <- normalize(demoData ,  norm_method="quant",
                        desiredQuantile = .9, toElt = "q_norm")  ## test with user parameter
demoData <- normalize(demoData, norm_method="quant") ## test defaults quantile = .75
housekeepers <- c('C1orf43','GPI','OAZ1','POLR2A','PSMB2','RAB7A',
                               'SDHA','SNRPD3','TBC1D10B','TPM4','TUBB','UBB')
#run aggregateCounts function on the data
target_demoData <- aggregateCounts(demoData)
#############################

########### Quantile Normalization test
#### Spec 1 check that normfactors are in pData of demoData
test_that("check that normfactors are in pData of demoData", {
  expect_true(length(demoData@phenoData@data[["q_norm_qFactors"]]) == dim(demoData@assayData$exprs)[2])
  expect_true(length(demoData@phenoData@data[["normFactors"]]) == dim(demoData@assayData$exprs)[2])
})

#### Spec 2 verify calculation of q90 norm factors
# Compute normalization factors from demoData
# compute 90% quantile count for samples and divide by geomean
thresh <- assayDataApply(demoData, 2, quantile, probs = .9)
expectedOutputData <- thresh/ngeoMean(thresh)
expectedOutputData <- unname(expectedOutputData)

# Extract normalized data from object
actualOutputData <- sData(demoData)[,"q_norm_qFactors"]
test_that("quantile norm factors are correct", {
  expect_equal(expectedOutputData, actualOutputData)
})

#### Spec 2b verify calculation of q75 norm factors
# compute 75% quantile norm factors for samples and divide by geomean
thresh <- assayDataApply(demoData, 2, quantile, probs = .75)
expectedOutputData <- thresh/ngeoMean(thresh)
expectedOutputData <- unname(expectedOutputData)

# Extract normalized data from object
actualOutputData <- sData(demoData)[,"normFactors"]
test_that("quantile norm factors for 75% are correct", {
  expect_equal(expectedOutputData, actualOutputData)
})


#### Spec 3 verify calculation of quantile norm values
# compute normalize count for samples and divide by geomean
thresh <- assayDataApply(demoData, 2, quantile, probs = .75)
norm_quant <- function(x){
  x <- x/(thresh/ngeoMean(thresh))
}
expectedOutputData <- t(assayDataApply(demoData, 1, norm_quant))

# Extract normalized data from object
actualOutputData <- assayData(demoData)[["exprs_norm"]]
test_that("quantile norm values are correct", {
  expect_equal(expectedOutputData, actualOutputData)
})


#### Spec 4 verify calculation of negative norm factors
#subset dataset for single panel
sub_target_demoData <- subset(target_demoData, subset = Module == "VnV_GeoMx_Hs_CTA_v1.2")
#call negative normalization
sub_target_demoData <- normalize(sub_target_demoData ,  norm_method="neg",
                      toElt = "neg_norm")
# compute negative norm factors for samples and divide by geomean
negfactors <- assayDataApply(negativeControlSubset(sub_target_demoData), 2, ngeoMean)
norm_neg <- function(x){
  x <- x/(negfactors/ngeoMean(negfactors))}
expectedOutputData <- t(assayDataApply(sub_target_demoData, 1, norm_neg))

# Extract normalized data from object
actualOutputData <- assayData(sub_target_demoData)[["neg_norm"]]
test_that("negative norm values are correct for single panel", {
  expect_equal(expectedOutputData, actualOutputData)
})

#### Spec 4a verify calculation of negative norm factors for multipanel
# compute neg_norm factors for multipanel
# call negative normalization
target_demoData <- normalize(target_demoData ,  norm_method="neg",
                                 toElt = "neg_norm")
expectedOutputData_1 <- pData(target_demoData)[["NegGeoMean_VnV_GeoMx_Hs_CTA_v1.2"]]/
  ngeoMean(pData(target_demoData)[["NegGeoMean_VnV_GeoMx_Hs_CTA_v1.2"]])
expectedOutputData_2 <- pData(target_demoData)[["NegGeoMean_Six-gene_test_v1_v1.1"]]/
ngeoMean(pData(target_demoData)[["NegGeoMean_Six-gene_test_v1_v1.1"]])
expectedOutputData <- (cbind("VnV_GeoMx_Hs_CTA_v1.2" = expectedOutputData_1,
                             "Six-gene_test_v1_v1.1" = expectedOutputData_2))
actualOutputData <- pData(target_demoData)[["neg_norm_negFactors"]]
dimnames(expectedOutputData) <- NULL
dimnames(actualOutputData) <- NULL

test_that("negative norm factors for multipanel are correct", {
  expect_equal(expectedOutputData, actualOutputData)
})

# check the normalized values for multipanel
if (length(unique(fData(target_demoData)[["Module"]])) > 1){
  # compute negative norm factors for samples and divide by geomean
  # subset data per module
  pool <- as.list(unique(fData(target_demoData)[["Module"]]))
  expectedNormMat <- lapply(pool, function (x) {
    poolSubSet <- subset(target_demoData, subset = Module == x)
    negfactors <- assayDataApply(negativeControlSubset(poolSubSet), 2, ngeoMean)
    norm_neg <- function(x){
      x <- x/(negfactors/ngeoMean(negfactors))
    }
    expectedNormMat_i <- t(assayDataApply(poolSubSet, 1, norm_neg))
  })
  expectedOutputData <- do.call(rbind, expectedNormMat)
}

# Extract normalized data from object
actualOutputData <- assayData(target_demoData)[["neg_norm"]]
test_that("negative norm values are correct", {
  expect_equal(expectedOutputData, actualOutputData)
})

negs <- which(fData(target_demoData)$CodeClass == "Negative")
test_that("Error is given if no negatives are in dataset", {
  expect_error(normalize(target_demoData[-negs[1],],  norm_method="neg",
                         toElt = "neg_norm"))
  expect_error(normalize(target_demoData[-negs[2],],  norm_method="neg",
                         toElt = "neg_norm"))
})

test_that("Error is given if dataset is not collapsed", {
  expect_error(normalize(demoData ,  norm_method="neg",
                         toElt = "neg_norm"))
})

########### Housekeeping Norm test
#### Spec 5 verify calculation of housekeeping norm factors
#call hk normalization
target_demoData <- normalize(target_demoData ,  norm_method="hk",
                             fromElt="exprs", toElt="hk_norm")
# compute hk norm factors for samples and divide by geomean
hksubset <- target_demoData[which(featureData(target_demoData)[["TargetName"]] %in% housekeepers),]
hkfactors <- assayDataApply(hksubset, 2, ngeoMean)
norm_hk <- function(x){
  x <- x/(hkfactors/ngeoMean(hkfactors))
}
expectedOutputData <- t(assayDataApply(target_demoData, 1, norm_hk))

# Extract normalized data from object
actualOutputData <- assayData(target_demoData)[["hk_norm"]]
test_that("hk norm values are correct", {
  expect_equal(expectedOutputData, actualOutputData)
})


# ########### Subtract ground Norm test
# #### Spec 6 verify calculation of subtract background norm factors
# #call subtract bg normalization
test_that("background subtraction throws warning on probe level data", {
  expect_warning(normalize(demoData , 
                           norm_method="subtractBackground", fromElt="exprs", toElt="bg_norm"))
})


target_demoData <- normalize(object = target_demoData , 
                             norm_method="subtractBackground", fromElt="exprs", 
                             toElt="bg_norm_byPanel", byPanel = TRUE)
target_demoData <- normalize(target_demoData , 
                             norm_method="subtractBackground", fromElt="exprs", 
                             toElt="bg_norm", byPanel = FALSE)

panels <- unique(fData(target_demoData)$Module)
# compute neg norm factors for samples and divide by geomean
bkcounts <- NULL
for(i in panels){
  negs <- assayDataApply(subset(negativeControlSubset(target_demoData), subset=Module == i), 2, ngeoMean)
  counts <- subset(target_demoData, subset=Module == i)
  bkcounts <- rbind(bkcounts, sweep(counts@assayData$exprs, 2, negs, "-"))
}

bkcounts[bkcounts < 0] <- 0

expectedOutputData <- bkcounts[featureNames(target_demoData), sampleNames(target_demoData)]

# Extract normalized data from object
actualOutputData <- assayData(target_demoData)[["bg_norm_byPanel"]]
test_that("bg subtract byPanel norm values are correct", {
  expect_equal(actualOutputData, expectedOutputData)
})

negs <- assayDataApply(negativeControlSubset(target_demoData), 2, ngeoMean)
counts <- target_demoData
bkcounts <- sweep(counts@assayData$exprs, 2, negs, "-")

bkcounts[bkcounts < 0] <- 0

expectedOutputData <- bkcounts[featureNames(target_demoData), sampleNames(target_demoData)]

# Extract normalized data from object
actualOutputData <- assayData(target_demoData)[["bg_norm"]]
test_that("bg subtract norm values are correct", {
  expect_equal(actualOutputData, expectedOutputData)
})

negs <- which(fData(target_demoData)$CodeClass == "Negative")
test_that("bg subtract errors with no negatives", {
  expect_error(normalize(target_demoData[-negs,] , 
                         norm_method="subtractBackground", fromElt="exprs", 
                         toElt="bg_norm", byPanel = FALSE))
})

