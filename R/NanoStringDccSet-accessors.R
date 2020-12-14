assayDataElement2 <- function(object, elt)
{
  if (elt %in% assayDataElementNames(object))
    assayDataElement(object, elt)
  else
    stop("'elt' not present in assayData(object)")
}

# sData Accessor
setGeneric("sData", signature = "object",
           function(object) standardGeneric("sData"))
setMethod("sData", "NanoStringDccSet",
          function(object) cbind(pData(object), pData(protocolData(object))))

# svarLabels Accessor
setGeneric("svarLabels", signature = "object",
           function(object) standardGeneric("svarLabels"))
setMethod("svarLabels", "NanoStringDccSet",
          function(object) c(varLabels(object), varLabels(protocolData(object))))

# dimLabels Accessor and Replacer
setMethod("dimLabels", "NanoStringDccSet", function(object) object@dimLabels)
setReplaceMethod("dimLabels", c("NanoStringDccSet", "character"),
                 function(object, value) {
                   object@dimLabels <- value
                   object
                 })

# design Accessor and Replacer
setMethod("design", "NanoStringDccSet", function(object) object@design)
setReplaceMethod("design", c("NanoStringDccSet", "formula"),
                 function(object, value) {
                   object@design <- value
                   object
                 })
setReplaceMethod("design", c("NanoStringDccSet", "ANY"),
                 function(object, value) {
                   object@design <- as.formula(value)
                   object
                 })
setReplaceMethod("design", c("NanoStringDccSet", "NULL"),
                 function(object, value) {
                   object@design <- NULL
                   object
                 })
