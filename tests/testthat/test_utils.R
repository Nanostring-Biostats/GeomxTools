library(GeomxTools)
library(testthat)
library(EnvStats)

fakeCounts <- 1:10
fakeCountsWith0 <- c(0, fakeCounts)
fakeCountsWithNA <- c(fakeCounts, NA)
fakeCountsWith0_NA <- c(0, fakeCountsWithNA)
fakeCountwith0 <- 0

#copied non-exported function from utils.R
thresholdValues <- function(x, thresh=0.5) {
  if (thresh <= 0) {
    warning("Parameter, thresh, cannot be set to less than or equal to 0. 
              The default threshold of 0.5 was used instead.")
    thresh <- 0.5
  }
  if (min(x, na.rm = TRUE) < thresh) {
    x <- x[!is.na(x)] + thresh
  }
  return(x)
}

# Spec 0: test thresholding behaves as expected:------
testthat::test_that("test thresholding behaves as expected", {
  expect_identical(fakeCounts, thresholdValues(fakeCounts))
  expect_identical(fakeCountsWith0+0.5, thresholdValues(fakeCountsWith0))
  expect_identical(fakeCountsWithNA, thresholdValues(fakeCountsWithNA))
  expect_identical(fakeCountsWith0+0.5, thresholdValues(fakeCountsWith0_NA))
  expect_identical(0.5, thresholdValues(fakeCountwith0))
})

# Spec 1: test that ngeoMean is calculated as expected:------
testthat::test_that("test that ngeoMean is calculated as expected", {
  expect_identical(ngeoMean(fakeCounts), EnvStats::geoMean(thresholdValues(fakeCounts), na.rm = TRUE))
  expect_identical(ngeoMean(fakeCountsWith0), EnvStats::geoMean(thresholdValues(fakeCountsWith0), na.rm = TRUE))
  expect_identical(ngeoMean(fakeCountsWithNA), EnvStats::geoMean(thresholdValues(fakeCountsWithNA), na.rm = TRUE))
  expect_identical(ngeoMean(fakeCountsWith0_NA), EnvStats::geoMean(thresholdValues(fakeCountsWith0_NA), na.rm = TRUE))
  expect_identical(ngeoMean(fakeCountwith0), EnvStats::geoMean(thresholdValues(fakeCountwith0), na.rm = TRUE))
})

# Spec 2: test that ngeoSD is calculated as expected:------
testthat::test_that("test that ngeoSD is calculated as expected", {
  expect_identical(ngeoSD(fakeCounts), EnvStats::geoSD(thresholdValues(fakeCounts), na.rm = TRUE))
  expect_identical(ngeoSD(fakeCountsWith0), EnvStats::geoSD(thresholdValues(fakeCountsWith0), na.rm = TRUE))
  expect_identical(ngeoSD(fakeCountsWithNA), EnvStats::geoSD(thresholdValues(fakeCountsWithNA), na.rm = TRUE))
  expect_identical(ngeoSD(fakeCountsWith0_NA), EnvStats::geoSD(thresholdValues(fakeCountsWith0_NA), na.rm = TRUE))
  expect_identical(ngeoSD(fakeCountwith0), EnvStats::geoSD(thresholdValues(fakeCountwith0), na.rm = TRUE))
})


# Spec 3: test that logtBase is calculated as expected:------
testthat::test_that("test that logtBase is calculated as expected", {
  expect_identical(logtBase(fakeCounts), log(thresholdValues(fakeCounts), base = 2))
  expect_identical(logtBase(fakeCountsWith0), log(thresholdValues(fakeCountsWith0), base = 2))
  expect_identical(logtBase(fakeCountsWithNA), log(thresholdValues(fakeCountsWithNA), base = 2))
  expect_identical(logtBase(fakeCountsWith0_NA), c(log(thresholdValues(fakeCountsWith0_NA), base = 2), NA))
  expect_identical(logtBase(fakeCountwith0), log(thresholdValues(fakeCountwith0), base = 2))
})


#These functions are tested in test_shiftCountsOne.R:
#shiftCountsOne()

