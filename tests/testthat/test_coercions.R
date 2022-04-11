### CONFIGURATION SECTION ###
library(GeomxTools)
library(testthat)

# load data
demoData <- readRDS(file= system.file("extdata","DSP_NGS_Example_Data", "demoData.rds", package = "GeomxTools"))
demoData <- shiftCountsOne(demoData)

noQC <- demoData

demoData <- setSegmentQCFlags(demoData, qcCutoffs = list(percentSaturation = 45))
demoData <- setBioProbeQCFlags(demoData)

#run aggregateCounts function on the data
target_demoData <- aggregateCounts(demoData)
noQC <- aggregateCounts(noQC)

# Spec 1: The coercion of a GeoMxSet object shall warn users when coercing non-normalized data but will coerce when forced.:------
# advised to only work on normalized data
test_that("GeomxSet object has been normalized - Seurat",{
    expect_error(to.Seurat(object = target_demoData, normData = "exprs"))
})
test_that("GeomxSet object has been normalized - SpatialExperiment",{
    expect_error(to.SpatialExperiment(object = target_demoData, normData = "exprs"))
})

# but will work on non-normalized data
test_that("if GeomxSet object hasn't been normalized, it can be forced to convert - Seurat",{
    expect_visible(to.Seurat(object = target_demoData, normData = "exprs", forceRaw = TRUE))
})
test_that("if GeomxSet object hasn't been normalized, it can be forced to convert - SpatialExperiment",{
    expect_visible(to.SpatialExperiment(object = target_demoData, normData = "exprs", forceRaw = TRUE))
})

target_demoData <- normalize(target_demoData, norm_method="quant")
noQC <- normalize(noQC, norm_method="quant")
#############################

# Spec 1: The coercion of a GeoMxSet object shall only occur on target level data.:------
# only works on target level data
test_that("coercion is only run on target level data - Seurat",{
    expect_error(to.Seurat(object = demoData, normData = "exprs_norm"))
})
test_that("coercion is only run on target level data - SpatialExperiment",{
    expect_error(to.SpatialExperiment(object = demoData, normData = "exprs_norm"))
})

# Spec 3: The coercion of a GeoMxSet object shall only occur when a valid norm 
#           data count matrix is provided by the user. :------
# needs normData
test_that("normData is required - Seurat",{
    expect_error(to.Seurat(object = target_demoData))
})
test_that("normData is required - SpatialExperiment",{
    expect_error(to.SpatialExperiment(object = target_demoData))
})


# normData must be valid
test_that("normData is a valid assayData name - Seurat",{
    expect_error(to.Seurat(object = target_demoData, normData = "invalid"))
})
test_that("normData is a valid assayData name - SpatialExperiment",{
    expect_error(to.SpatialExperiment(object = target_demoData, normData = "invalid"))
})

library(Seurat)
library(SpatialExperiment)

seurat_object <- to.Seurat(object = target_demoData, normData = "exprs_norm", ident = "cell_line")
noQC_seurat_object <- to.Seurat(object = noQC, normData = "exprs_norm", ident = "cell_line")

spe_object <- to.SpatialExperiment(object = target_demoData, normData = "exprs_norm")
noQC_spe_object <- to.SpatialExperiment(object = noQC, normData = "exprs_norm")


# Spec 4: The coercion of a GeoMxSet object shall copy the wanted data from the 
#           GeoMxSet object to the correct location in the coerced object. :------
# ident is equal to column given 
test_that("Seurat Ident is equal to given column",{
    expect_true(all(Seurat::Idents(seurat_object) == sData(target_demoData)$cell_line))
})

# count matrix in correct location
test_that("Count matrix is in the correct location - Seurat", {
    expect_true(all(GetAssayData(seurat_object, "counts") == 
                      assayDataElement(target_demoData, "exprs_norm")))
    expect_identical(colnames(GetAssayData(seurat_object, "counts")), 
                     colnames(assayDataElement(target_demoData, "exprs_norm")))
    expect_identical(rownames(GetAssayData(seurat_object, "counts")), 
                     rownames(assayDataElement(target_demoData, "exprs_norm")))
})
test_that("Count matrix is in the correct location - SpatialExperiment", {
    expect_true(all(assay(spe_object, i = "GeoMx") == 
                      assayDataElement(target_demoData, "exprs_norm")))
    expect_identical(colnames(assay(spe_object, i = "GeoMx")), 
                     colnames(assayDataElement(target_demoData, "exprs_norm")))
    expect_identical(rownames(assay(spe_object, i = "GeoMx")), 
                     rownames(assayDataElement(target_demoData, "exprs_norm")))
})

sequencingMetrics <- c("FileVersion", "SoftwareVersion", "Date", "Plate_ID", 
                       "Well", "SeqSetId", "Raw", "Trimmed", "Stitched", 
                       "Aligned", "umiQ30", "rtsQ30", "DeduplicatedReads", 
                       "NTC_ID", "NTC", "Trimmed (%)", "Stitched (%)", 
                       "Aligned (%)", "Saturated (%)")

QCMetrics <- "QCFlags"

identMetrics <- colnames(sData(target_demoData))[!colnames(sData(target_demoData)) %in% 
                                                   c(sequencingMetrics, QCMetrics)]

# pheno data in correct spot
test_that("pheno data is in the correct location - Seurat", {
    expect_true(all(seurat_object[[]][,colnames(seurat_object[[]]) %in% 
                                              gsub("\\W", ".", identMetrics)] == 
                        sData(target_demoData)[,identMetrics]))
    expect_true(all(noQC_seurat_object[[]][,colnames(noQC_seurat_object[[]]) %in% 
                                                   gsub("\\W", ".", identMetrics)] == 
                        sData(noQC)[,identMetrics]))
})
test_that("pheno data is in the correct location - SpatialExperiment", {
    expect_true(all(all(colData(spe_object)[,colnames(colData(spe_object)) %in% 
                                              identMetrics] == 
                            sData(target_demoData)[,identMetrics])))
    expect_true(all(all(colData(noQC_spe_object)[,colnames(colData(noQC_spe_object)) %in% 
                                                   identMetrics] == 
                            sData(noQC)[,identMetrics])))
})

# sequencing metrics in correct spot
test_that("sequencing metrics are in the correct location - Seurat", {
    expect_true(all(Misc(seurat_object)[["sequencingMetrics"]] == 
                      sData(target_demoData)[,colnames(sData(target_demoData)) %in% 
                                               sequencingMetrics]))
    expect_true(all(Misc(noQC_seurat_object)[["sequencingMetrics"]] == 
                      sData(noQC)[,colnames(sData(noQC)) %in% sequencingMetrics]))
})
test_that("sequencing metrics are in the correct location - SpatialExperiment", {
    expect_true(all(all(spe_object@metadata$sequencingMetrics == 
                          sData(target_demoData)[,colnames(sData(target_demoData)) %in% 
                                                   sequencingMetrics])))
    expect_true(all(all(noQC_spe_object@metadata$sequencingMetrics == 
                          sData(noQC)[,colnames(sData(noQC)) %in% sequencingMetrics])))
})

# sequencing metrics in correct spot
test_that("QC metrics are in the correct location - Seurat", {
    expect_true(all(Misc(seurat_object)[["QCMetrics"]][["QCFlags"]] == 
                      sData(target_demoData)[,colnames(sData(target_demoData)) %in% 
                                               QCMetrics]))
    expect_true(all(Misc(noQC_seurat_object)[["QCMetrics"]] == 
                      sData(noQC)[,colnames(sData(noQC)) %in% QCMetrics]))
})
test_that("QC metrics are in the correct location - SpatialExperiment", {
    expect_true(all(all(spe_object@metadata$QCMetrics$QCFlags == 
                          sData(target_demoData)[,colnames(sData(target_demoData)) %in% 
                                                   QCMetrics])))
    expect_true(all(all(noQC_spe_object@metadata$QCMetrics == 
                          sData(noQC)[,colnames(sData(noQC)) %in% QCMetrics])))
})


# experiment data in correct spot
test_that("experiment data are in the correct location - Seurat", {
    expect_identical(Misc(seurat_object)[which(!names(Misc(seurat_object)) %in% 
                                                c("sequencingMetrics", "QCMetrics"))], 
                     otherInfo(experimentData(target_demoData)))
    expect_identical(Misc(noQC_seurat_object)[which(!names(Misc(noQC_seurat_object)) %in% 
                                                     c("sequencingMetrics", "QCMetrics"))], 
                     otherInfo(experimentData(noQC)))
})
test_that("experiment data are in the correct location - SpatialExperiment", {
    expect_identical(spe_object@metadata[which(!names(spe_object@metadata)  %in% 
                                                 c("sequencingMetrics", "QCMetrics"))], 
                     otherInfo(experimentData(target_demoData)))
    expect_identical(noQC_spe_object@metadata[which(!names(noQC_spe_object@metadata)  %in% 
                                                      c("sequencingMetrics", "QCMetrics"))], 
                     otherInfo(experimentData(noQC)))
})

# feature metadata in correct spot
test_that("feature data is in the correct location - Seurat", {
    expect_true(all(seurat_object@assays$GeoMx@meta.features == 
                      fData(target_demoData)))
})
test_that("feature data is in the correct location - SpatialExperiment", {
    expect_true(all(all(rowData(spe_object) == fData(target_demoData))))
})

pData(target_demoData)$Xcoord <- runif(nrow(sData(target_demoData)), 0, 100)
pData(target_demoData)$Ycoord <- runif(nrow(sData(target_demoData)), 0, 100)

seurat_object <- to.Seurat(object = target_demoData, normData = "exprs_norm", 
                           ident = "cell_line", coordinates = c("Xcoord", 
                                                                "Ycoord"))
spe_object <- to.SpatialExperiment(object = target_demoData, 
                                   normData = "exprs_norm", 
                                   coordinates = c("Xcoord", "Ycoord"))

# coordinates there if given
test_that("coordinates, when given, are in the correct location - Seurat", {
    expect_identical(seurat_object@images$image@coordinates,
                     sData(target_demoData)[,c("Xcoord", "Ycoord")])
})
test_that("coordinates, when given, are in the correct location - SpatialExperiment", {
    expect_true(all(all(spatialCoords(spe_object) == 
                          sData(target_demoData)[,c("Xcoord", "Ycoord")])))
})

# Spec 4: The coercion of a GeoMxSet object shall warn users when the coordinate 
#           column names are not valid. :------

# errors if coordinates columns are not valid
test_that("invalid coordinates give an error - Seurat", {
    expect_error(to.Seurat(object = target_demoData, normData = "exprs_norm", 
                           ident = "cell_line", coordinates = c("Invalid", 
                                                                "Ycoord")))
    expect_error(to.Seurat(object = target_demoData, normData = "exprs_norm", 
                           ident = "cell_line", coordinates = c("Xcoord", 
                                                                "Invalid")))
    expect_error(to.Seurat(object = target_demoData, normData = "exprs_norm", 
                           ident = "cell_line", coordinates = c("Ycoord")))
})
test_that("invalid coordinates give an error - SpatialExperiment", {
    expect_error(to.SpatialExperiment(object = target_demoData, 
                                      normData = "exprs_norm", 
                                      coordinates = c("Invalid", "Ycoord")))
    expect_error(to.SpatialExperiment(object = target_demoData, 
                                      normData = "exprs_norm", 
                                      coordinates = c("Xcoord", "Invalid")))
    expect_error(to.SpatialExperiment(object = target_demoData, 
                                      normData = "exprs_norm", 
                                      coordinates = c("Ycoord")))
})

