### CONFIGURATION SECTION ###
library(GeomxTools)
library(testthat)

datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
DCCFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
PKCFiles <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
SampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

demoData <-
  suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles, # QuickBase: readNanoStringGeomxSet, need to change it.
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

#Shift counts to one to mimic how DSPDA handles zero counts
demoData <- shiftCountsOne(demoData, elt="exprs", useDALogic=TRUE) 


#QC flags
demoData <- setSegmentQCFlags(demoData)
demoData <- setBioProbeQCFlags(demoData)

#aggregate counts
demoData <- aggregateCounts(demoData)

#normalize
demoData <- normalize(demoData , data_type="RNA", norm_method="quant",
                      desiredQuantile = .9, toElt = "q_norm")

writeNanoStringGeoMxSet(demoData, dir = "./")
