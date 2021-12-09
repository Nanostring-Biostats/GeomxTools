### CONFIGURATION SECTION ###
library(GeomxTools)
library(testthat)

# load data
demoData <- readRDS(file= system.file("extdata","DSP_NGS_Example_Data", "demoData.rds", package = "GeomxTools"))
demoData <- shiftCountsOne(demoData)

#run aggregateCounts function on the data
target_demoData <- aggregateCounts(demoData)

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

demoData <- normalize(demoData, norm_method="quant") 
target_demoData <- normalize(target_demoData, norm_method="quant") 
#############################

# only works on target level data
test_that("coercion is only run on target level data - Seurat",{
    expect_error(to.Seurat(object = demoData, normData = "exprs_norm"))
})
test_that("coercion is only run on target level data - SpatialExperiment",{
    expect_error(to.SpatialExperiment(object = demoData, normData = "exprs_norm"))
})

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
spe_object <- to.SpatialExperiment(object = target_demoData, normData = "exprs_norm")

# ident is equal to column given 
test_that("Seurat Ident is equal to given column",{
    expect_true(all(Seurat::Idents(seurat_object) == sData(target_demoData)$cell_line))
})

# count matrix in correct location
test_that("Count matrix is in the correct location - Seurat", {
    expect_true(all(seurat_object@assays$GeoMx@counts == target_demoData@assayData$exprs_norm))
    expect_identical(colnames(seurat_object@assays$GeoMx@counts), colnames(target_demoData@assayData$exprs_norm))
    expect_identical(rownames(seurat_object@assays$GeoMx@counts), rownames(target_demoData@assayData$exprs_norm))
})
test_that("Count matrix is in the correct location - SpatialExperiment", {
    expect_true(all(assay(spe_object, i = "GeoMx") == target_demoData@assayData$exprs_norm))
    expect_identical(colnames(assay(spe_object, i = "GeoMx")), colnames(target_demoData@assayData$exprs_norm))
    expect_identical(rownames(assay(spe_object, i = "GeoMx")), rownames(target_demoData@assayData$exprs_norm))
})

sequencingMetrics <- c("FileVersion", "SoftwareVersion", "Date", "Plate_ID", "Well", "SeqSetId", "Raw", "Trimmed", 
                       "Stitched", "Aligned", "umiQ30", "rtsQ30", "DeduplicatedReads", "NTC_ID", "NTC")

# pheno data in correct spot
test_that("pheno data is in the correct location - Seurat", {
    expect_true(all(seurat_object@meta.data[,-c(1:3)] == sData(target_demoData)[,!colnames(sData(target_demoData)) %in% sequencingMetrics]))
})
test_that("pheno data is in the correct location - SpatialExperiment", {
    expect_true(all(all(colData(spe_object) == sData(target_demoData)[,!colnames(sData(target_demoData)) %in% sequencingMetrics])))
})

# sequencing metrics in correct spot
test_that("sequencing metrics are in the correct location - Seurat", {
    expect_true(all(seurat_object@misc$sequencingMetrics == sData(target_demoData)[,colnames(sData(target_demoData)) %in% sequencingMetrics]))
})
test_that("sequencing metrics are in the correct location - SpatialExperiment", {
    expect_true(all(all(spe_object@metadata$sequencingMetrics == sData(target_demoData)[,colnames(sData(target_demoData)) %in% sequencingMetrics])))
})

# experiment data in correct spot
test_that("sequencing metrics are in the correct location - Seurat", {
    expect_identical(seurat_object@misc[which(names(seurat_object@misc) != "sequencingMetrics")],target_demoData@experimentData@other)
})
test_that("sequencing metrics are in the correct location - SpatialExperiment", {
    expect_identical(spe_object@metadata[which(names(spe_object@metadata) != "sequencingMetrics")],target_demoData@experimentData@other)
})

# feature metadata in correct spot
test_that("feature data is in the correct location - Seurat", {
    expect_true(all(seurat_object@assays$GeoMx@meta.features == fData(target_demoData)))
})
test_that("feature data is in the correct location - SpatialExperiment", {
    expect_true(all(all(rowData(spe_object) == fData(target_demoData))))
})

target_demoData@phenoData$Xcoord <- runif(nrow(sData(target_demoData)), 0, 100)
target_demoData@phenoData$Ycoord <- runif(nrow(sData(target_demoData)), 0, 100)

seurat_object <- to.Seurat(object = target_demoData, normData = "exprs_norm", ident = "cell_line", coordinates = c("Xcoord", "Ycoord"))
spe_object <- to.SpatialExperiment(object = target_demoData, normData = "exprs_norm", coordinates = c("Xcoord", "Ycoord"))

# coordinates there if given
test_that("coordinates, when given, are in the correct location - Seurat", {
    expect_identical(seurat_object@images$image@coordinates,sData(target_demoData)[,c("Xcoord", "Ycoord")])
})
test_that("coordinates, when given, are in the correct location - SpatialExperiment", {
    expect_true(all(all(spatialCoords(spe_object) == sData(target_demoData)[,c("Xcoord", "Ycoord")])))
})

# errors if coordinates columns are not valid
test_that("coordinates, when given, are in the correct location - Seurat", {
    expect_error(to.Seurat(object = target_demoData, normData = "exprs_norm", ident = "cell_line", coordinates = c("Invalid", "Ycoord")))
    expect_error(to.Seurat(object = target_demoData, normData = "exprs_norm", ident = "cell_line", coordinates = c("Xcoord", "Invalid")))
    expect_error(to.Seurat(object = target_demoData, normData = "exprs_norm", ident = "cell_line", coordinates = c("Ycoord")))
})
test_that("coordinates, when given, are in the correct location - SpatialExperiment", {
    expect_error(to.SpatialExperiment(object = target_demoData, normData = "exprs_norm", coordinates = c("Invalid", "Ycoord")))
    expect_error(to.SpatialExperiment(object = target_demoData, normData = "exprs_norm", coordinates = c("Xcoord", "Invalid")))
    expect_error(to.SpatialExperiment(object = target_demoData, normData = "exprs_norm", coordinates = c("Ycoord")))
})

