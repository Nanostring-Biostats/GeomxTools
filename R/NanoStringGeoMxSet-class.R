setClassUnion("formulaOrNULL", c("formula", "NULL"))

# Class definition
setClass("NanoStringGeoMxSet",
         contains = "NanoStringRccSet",
         slots = c(dimLabels = "character",
                   signatures = "SignatureSet",
                   design = "formulaOrNULL",
                   featureType = "character",
                   analyte = "character"),
         prototype = prototype(
             new("VersionedBiobase",
                 versions = c(classVersion("ExpressionSet"),
                              NanoStringGeoMxSet = "3.0.0")),
             signatures = SignatureSet(),
             design = NULL,
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

# for when not inheriting from NCtools
# setMethod("updateObject", signature = "NanoStringGeoMxSet",
#     function(object){
#         if(!"analyte" %in% names(getObjectSlots(object))){
#           object@analyte <- "RNA"
#           
#           object@.__classVersion__$NanoStringGeoMxSet <- "2.1.6"
#         }
# 
#         return(object)
#     })

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
         analyte = "RNA",
         ...)
{
  assayData <- assayDataNew(exprs = matrix(integer(), nrow = 0L, ncol = 0L))
  methods::callGeneric(assayData = assayData, phenoData = phenoData,
              featureData = featureData, experimentData = experimentData,
              annotation = annotation, protocolData = protocolData,
              dimLabels = dimLabels, signatures = signatures, design = design,
              featureType = featureType, analyte = analyte,
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
       analyte = analyte,
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
         analyte = "RNA",
         ...)
{
  assayData <- assayDataNew(exprs = assayData)
  methods::callGeneric(assayData = assayData, phenoData = phenoData,
              featureData = featureData, experimentData = experimentData,
              annotation = annotation, protocolData = protocolData,
              dimLabels = dimLabels, signatures = signatures, design = design,
              featureType = featureType, analyte = analyte,
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
         analyte = "RNA",
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
              analyte = analyte, 
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
         analyte = "RNA",
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
              analyte = analyte,
              ...)
})

# Coersion
setAs("ExpressionSet", "NanoStringGeoMxSet",
      function(from) NanoStringGeoMxSet(from))
