setClassUnion("formulaOrNULL", c("formula", "NULL"))

# Class definition
setClass("NanoStringGeoMxSet",
         contains = "NanoStringRccSet",
         slots = c(dimLabels = "character",
                   signatures = "SignatureSet",
                   design = "formulaOrNULL",
                   featureType = "character"),
         prototype = prototype(
             new("VersionedBiobase",
                 versions = c(classVersion("ExpressionSet"),
                              NanoStringGeoMxSet = "2.0.0")),
             signatures = SignatureSet(),
             design = NULL,
             featureType = "Probe"))

# Show method
setMethod("show", signature = "NanoStringGeoMxSet",
function(object) {
    methods::callNextMethod(object)
    cat("feature: ")
    cat(featureType(object))
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
         ...)
    standardGeneric("NanoStringGeoMxSet"),
signature = "assayData")

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
         ...)
{
  assayData <- assayDataNew(exprs = matrix(integer(), nrow = 0L, ncol = 0L))
  methods::callGeneric(assayData = assayData, phenoData = phenoData,
              featureData = featureData, experimentData = experimentData,
              annotation = annotation, protocolData = protocolData,
              dimLabels = dimLabels, signatures = signatures, design = design,
              featureType = featureType,
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
         ...)
{
  new2("NanoStringGeoMxSet",
       assayData = assayData,
       phenoData = phenoData,
       featureData = featureData,
       experimentData = experimentData,
       annotation = annotation,
       protocolData = protocolData,
       dimLabels = dimLabels,
       signatures = signatures,
       design = design,
       featureType = featureType,
       ...)
})

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
         ...)
{
  assayData <- assayDataNew(exprs = assayData)
  methods::callGeneric(assayData = assayData, phenoData = phenoData,
              featureData = featureData, experimentData = experimentData,
              annotation = annotation, protocolData = protocolData,
              dimLabels = dimLabels, signatures = signatures, design = design,
              featureType = featureType,
              ...)
})

setMethod("NanoStringGeoMxSet", "ExpressionSet",
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
         ...)
{
  methods::callGeneric(assayData = copyEnv(assayData(assayData)),
              phenoData = Biobase::phenoData(assayData),
              featureData = Biobase::featureData(assayData),
              experimentData = Biobase::experimentData(assayData),
              annotation = Biobase::annotation(assayData),
              protocolData = Biobase::protocolData(assayData),
              dimLabels = dimLabels,
              signatures = signatures,
              design = design,
              featureType = featureType,
              ...)
})

setMethod("NanoStringGeoMxSet", "NanoStringGeoMxSet",
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
         ...)
{
  methods::callGeneric(assayData = copyEnv(assayData(assayData)),
              phenoData = Biobase::phenoData(assayData),
              featureData = Biobase::featureData(assayData),
              experimentData = Biobase::experimentData(assayData),
              annotation = Biobase::annotation(assayData),
              protocolData = Biobase::protocolData(assayData),
              dimLabels = dimLabels(assayData),
              signatures = signatures(assayData),
              design = design(assayData),
              featureType = featureType,
              ...)
})

# Coersion
setAs("ExpressionSet", "NanoStringGeoMxSet",
      function(from) NanoStringGeoMxSet(from))
