library(GeomxTools)

datadir <- "~/NAS_data/mgriswold/GeomxTools_protein/"
DCCFiles <- dir(paste0(datadir, "/shilah_DCCs_100621"), pattern=".dcc$", full.names=TRUE)
PKCFiles <- dir(paste0(datadir, "/pkcs"), pattern=".pkc$", full.names=TRUE)
SampleAnnotationFile <- file.path(datadir, "ProteoGen_Seg_Annotation.xlsx")

testData <-
  suppressWarnings(readNanoStringGeoMxSet(dccFiles = DCCFiles, 
                                          pkcFiles = PKCFiles,
                                          phenoDataFile = SampleAnnotationFile,
                                          phenoDataSheet = "Annotations",
                                          phenoDataDccColName = "Sample_ID",
                                          protocolDataColNames = c("Tissue",
                                                                   "Segment_Type")))
aggTestData <- aggregateCounts(testData)

igg.names <- igg_names(aggTestData)
hk.names <- hk_names(aggTestData)

plot_concordance(igg.names, aggTestData, "Segment_Type")
plot_concordance(igg.names, aggTestData, "Tissue")

normfactors <- compute_normalization_factors(object = aggTestData, 
                                             # area = "AOI.Size.um2", 
                                             nuclei = "Nuclei.Counts")

plot_normFactor_concordance(aggTestData, "Segment_Type", normfactors)

qc_protein_signal(object = aggTestData,
                  neg.names = igg.names,
                  targetAnnotations = fData(aggTestData))


