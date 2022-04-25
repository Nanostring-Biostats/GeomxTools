# Class definition
.NanoStringGeoMxSet <- setClass("NanoStringGeoMxSet",
    contains = "NanoStringExperiment",
    slots = c(featureType = "character",
        analyte = "character"),
    prototype = prototype(
        new("NanoStringExperiment"),
            .__classVersion__ = 
                c(NanoStringExperiment = 
                    paste(packageVersion("NanoStringExperiment"), 
                        collapse="."),
                GeomxTools = 
                    paste(packageVersion("GeomxTools"), 
                        collapse=".")),
            featureType = "Probe",
            analyte = "RNA"))

# Show method
setMethod("show", signature = "NanoStringGeoMxSet",
function(object) {
    methods::callNextMethod(object)
    cat("feature: ")
    cat(featureType(object))
    cat("\n")
    cat("analyte: ")
    cat(analyte(object))
    cat("\n")
})

setMethod("updateObject", signature = "NanoStringGeoMxSet",
    function(object){
        return(NanoStringGeoMxSet(object))
    })

# Constructors
setGeneric("NanoStringGeoMxSet",
function(assayData,
         phenoData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = Biobase::annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = Biobase::MIAME(),
         annotation = character(),
         protocolData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         dimLabels = c("TargetName", "SampleID"),
         signatures = SignatureSet(),
         design = NULL,
         featureType = "Probe",
         analyte = "RNA",
         ...){
    standardGeneric("NanoStringGeoMxSet")},
signature = "assayData")

setMethod("NanoStringGeoMxSet", "matrix",
function(assayData,
         phenoData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = Biobase::annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = Biobase::MIAME(),
         annotation = character(),
         protocolData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         dimLabels = c("TargetName", "SampleID"),
         signatures = SignatureSet(),
         design = NULL,
         featureType = "Probe",
         analyte = "RNA",
         ...)
{
    nse <- NanoStringExperiment(
        assayData = assayData,
        phenoData = phenoData,
        featureData = featureData,
        experimentData = experimentData,
        annotation = annotation,
        protocolData = protocolData,
        dimLabels = dimLabels,
        signatures = signatures,
        design = design, ...)
    .NanoStringGeoMxSet(nse,
        featureType = featureType,
        analyte = analyte)
})

setMethod("NanoStringGeoMxSet", "missing",
function(assayData,
         phenoData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = Biobase::annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = Biobase::MIAME(),
         annotation = character(),
         protocolData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         dimLabels = c("TargetName", "SampleID"),
         signatures = SignatureSet(),
         design = NULL,
         featureType = "Probe",
         analyte = "RNA",
         ...)
{
    methods::callGeneric(
        assayData = matrix(integer(), nrow = 0L, ncol = 0L), 
        phenoData = phenoData,
        featureData = featureData,
        experimentData = experimentData,
        annotation = annotation,
        protocolData = protocolData,
        dimLabels = dimLabels,
        signatures = signatures,
        design = design,
        featureType = featureType,
        analyte = analyte,
        ...)
})

setMethod("NanoStringGeoMxSet", "environment",
function(assayData,
         phenoData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = Biobase::annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = Biobase::MIAME(),
         annotation = character(),
         protocolData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         dimLabels = c("TargetName", "SampleID"),
         signatures = SignatureSet(),
         design = NULL,
         featureType = "Probe",
         analyte = "RNA",
         ...)
{
    nse <- NanoStringExperiment(
        assayData = assayData,
        phenoData = phenoData,
        featureData = featureData,
        experimentData = experimentData,
        annotation = annotation,
        protocolData = protocolData,
        dimLabels = dimLabels,
        signatures = signatures,
        design = design,
        ...)
    .NanoStringGeoMxSet(nse,
        featureType = featureType,
        analyte = analyte)
})

setMethod("NanoStringGeoMxSet", "ExpressionSet",
function(assayData,
         dimLabels = c("TargetName", "SampleID"),
         signatures = SignatureSet(),
         design = NULL,
         featureType = "Probe",
         analyte = "RNA",
         ...)
{
    se <- SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(
        assayData)
    methods::callGeneric(
        assayData = se,
        dimLabels = dimLabels,
        signatures = signatures,
        design = design,
        featureType = featureType,
        analyte = analyte,
        ...)
})

setMethod("NanoStringGeoMxSet", "NanoStringExperiment",
function(assayData,
         featureType = "Probe",
         analyte = "RNA",
         ...)
{
    .NanoStringGeoMxSet(assayData,
        featureType = featureType,
        analyte = analyte)
})

setMethod("NanoStringGeoMxSet", "SummarizedExperiment",
function(assayData,
         phenoData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = Biobase::annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = Biobase::MIAME(),
         annotation = character(),
         protocolData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         dimLabels = c("TargetName", "SampleID"),
         signatures = SignatureSet(),
         design = NULL,
         featureType = "Probe",
         analyte = "RNA",
         ...)
{
    nse <- NanoStringExperiment(
        assayData = assayData,
        phenoData = phenoData,
        featureData = featureData,
        experimentData = experimentData,
        annotation = annotation,
        protocolData = protocolData,
        dimLabels = dimLabels,
        signatures = signatures,
        design = design,
        ...)
    .NanoStringGeoMxSet(nse,
        featureType = featureType,
        analyte = analyte)
})

# Copy constructor
setMethod("NanoStringGeoMxSet", "NanoStringGeoMxSet",
function(assayData,
         featureType = "Probe",
         analyte = "RNA",
         ...)
{
    if (featureType %in% slotNames(object)) {
        featureType <- featureType(object)
    }
    if (analyte %in% slotNames(object)) {
        analyte <- analyte(object)
    }
    methods::callGeneric(
        assayData = assayData(assayData),
        phenoData = Biobase::phenoData(assayData),
        featureData = Biobase::featureData(assayData),
        experimentData = Biobase::experimentData(assayData),
        annotation = Biobase::annotation(assayData),
        protocolData = Biobase::protocolData(assayData),
        dimLabels = dimLabels(assayData),
        signatures = signatures(assayData),
        design = design(assayData),
        featureType = featureType,
        analyte = analyte,
        ...)
})

# Coercion
setAs("ExpressionSet", "NanoStringGeoMxSet",
    function(from) NanoStringGeoMxSet(from))

setAs("SummarizedExperiment", "NanoStringGeoMxSet",
    function(from) NanoStringGeoMxSet(from))

setAs("NanoStringExperiment", "NanoStringGeoMxSet",
    function(from) NanoStringGeoMxSet(from))
