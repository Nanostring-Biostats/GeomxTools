setClassUnion("formulaOrNULL", c("formula", "NULL"))

# Class definition
setClass("NanoStringGeomxSet",
         contains = "ExpressionSet",
         slots = c(dimLabels = "character",
                   signatures = "SignatureSet",
                   design = "formulaOrNULL",
                   featureType = "character"),
         prototype = prototype(
             new("VersionedBiobase",
                 versions = c(classVersion("ExpressionSet"),
                              NanoStringGeomxSet = "1.0.0")),
             signatures = SignatureSet(),
             design = NULL,
             featureType = "Probe"))

# Show method
setMethod("show", signature = "NanoStringGeomxSet",
function(object) {
    methods::callNextMethod(object)
    cat("feature: ")
    cat(featureType(object))
    cat("signature: ")
    if (length(signatures(object)) == 0L)
        cat("none\n")
    else
        cat("use 'signatures(object)'")
})

# Constructors
setGeneric("NanoStringGeomxSet",
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
    standardGeneric("NanoStringGeomxSet"),
signature = "assayData")

setMethod("NanoStringGeomxSet", "missing",
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

setMethod("NanoStringGeomxSet", "environment",
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
  new2("NanoStringGeomxSet",
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

setMethod("NanoStringGeomxSet", "matrix",
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

setMethod("NanoStringGeomxSet", "ExpressionSet",
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

setMethod("NanoStringGeomxSet", "NanoStringGeomxSet",
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
setAs("ExpressionSet", "NanoStringGeomxSet",
      function(from) NanoStringGeomxSet(from))
