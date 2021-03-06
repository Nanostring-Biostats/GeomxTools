setClassUnion("formulaOrNULL", c("formula", "NULL"))

# Class definition
setClass("NanoStringGeomxSet",
         contains = "ExpressionSet",
         slots = c(dimLabels = "character",
                   signatures = "SignatureSet",
                   design = "formulaOrNULL"),
         prototype = prototype(
           new("VersionedBiobase",
               versions = c(classVersion("ExpressionSet"),
                            NanoStringGeomxSet = "1.0.0")),
           signatures = SignatureSet(),
           design = NULL))

# Show method
setMethod("show", signature = "NanoStringGeomxSet",
function(object) {
  methods::callNextMethod(object)
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
         dimLabels = c("GeneName", "SampleID"),
         signatures = SignatureSet(),
         design = NULL,
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
         dimLabels = c("GeneName", "SampleID"),
         signatures = SignatureSet(),
         design = NULL,
         ...)
{
  assayData <- assayDataNew(exprs = matrix(integer(), nrow = 0L, ncol = 0L))
  methods::callGeneric(assayData = assayData, phenoData = phenoData,
              featureData = featureData, experimentData = experimentData,
              annotation = annotation, protocolData = protocolData,
              dimLabels = dimLabels, signatures = signatures, design = design,
              ...)
})

setMethod("NanoStringGeomxSet", "environment",
function(assayData,
         phenoData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = Biobase::annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = Biobase::MIAME(),
         annotation = character(),
         protocolData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         dimLabels = c("GeneName", "SampleID"),
         signatures = SignatureSet(),
         design = NULL,
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
       ...)
})

setMethod("NanoStringGeomxSet", "matrix",
function(assayData,
         phenoData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = Biobase::annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = Biobase::MIAME(),
         annotation = character(),
         protocolData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         dimLabels = c("GeneName", "SampleID"),
         signatures = SignatureSet(),
         design = NULL,
         ...)
{
  assayData <- assayDataNew(exprs = assayData)
  methods::callGeneric(assayData = assayData, phenoData = phenoData,
              featureData = featureData, experimentData = experimentData,
              annotation = annotation, protocolData = protocolData,
              dimLabels = dimLabels, signatures = signatures, design = design,
              ...)
})

setMethod("NanoStringGeomxSet", "ExpressionSet",
function(assayData,
         phenoData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = Biobase::annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = Biobase::MIAME(),
         annotation = character(),
         protocolData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         dimLabels = c("GeneName", "SampleID"),
         signatures = SignatureSet(),
         design = NULL,
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
              ...)
})

setMethod("NanoStringGeomxSet", "NanoStringGeomxSet",
function(assayData,
         phenoData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         featureData = Biobase::annotatedDataFrameFrom(assayData, byrow = TRUE),
         experimentData = Biobase::MIAME(),
         annotation = character(),
         protocolData = Biobase::annotatedDataFrameFrom(assayData, byrow = FALSE),
         dimLabels = c("GeneName", "SampleID"),
         signatures = SignatureSet(),
         design = NULL,
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
              ...)
})

# Coersion
setAs("ExpressionSet", "NanoStringGeomxSet",
      function(from) NanoStringGeomxSet(from))
