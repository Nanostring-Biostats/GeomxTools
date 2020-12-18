assayDataElement2 <- function(object, elt)
{
  if (elt %in% Biobase::assayDataElementNames(object))
    Biobase::assayDataElement(object, elt)
  else
    stop("'elt' not present in assayData(object)")
}

# sData Accessor
setGeneric("sData", signature = "object",
           function(object) standardGeneric("sData"))
setMethod("sData", "NanoStringGeomxSet",
          function(object) cbind(pData(object), pData(protocolData(object))))

# svarLabels Accessor
setGeneric("svarLabels", signature = "object",
           function(object) standardGeneric("svarLabels"))
setMethod("svarLabels", "NanoStringGeomxSet",
          function(object) c(varLabels(object), varLabels(protocolData(object))))

# dimLabels Accessor and Replacer
setMethod("dimLabels", "NanoStringGeomxSet", function(object) object@dimLabels)
setReplaceMethod("dimLabels", c("NanoStringGeomxSet", "character"),
                 function(object, value) {
                   object@dimLabels <- value
                   object
                 })

# design Accessor and Replacer
setMethod("design", "NanoStringGeomxSet", function(object) object@design)
setReplaceMethod("design", c("NanoStringGeomxSet", "formula"),
                 function(object, value) {
                   object@design <- value
                   object
                 })
setReplaceMethod("design", c("NanoStringGeomxSet", "ANY"),
                 function(object, value) {
                   object@design <- stats::as.formula(value)
                   object
                 })
setReplaceMethod("design", c("NanoStringGeomxSet", "NULL"),
                 function(object, value) {
                   object@design <- NULL
                   object
                 })


