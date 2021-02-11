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
                   return(object)
                 })

# design Accessor and Replacer
setMethod("design", "NanoStringGeomxSet", function(object) object@design)
setReplaceMethod("design", c("NanoStringGeomxSet", "formula"),
                 function(object, value) {
                   object@design <- value
                   return(object)
                 })
setReplaceMethod("design", c("NanoStringGeomxSet", "ANY"),
                 function(object, value) {
                   object@design <- stats::as.formula(value)
                   return(object)
                 })
setReplaceMethod("design", c("NanoStringGeomxSet", "NULL"),
                 function(object, value) {
                   object@design <- NULL
                   return(object)
                 })

setGeneric("featureType", signature = "object",
           function(object, ...) standardGeneric("featureType"))
setMethod("featureType", "NanoStringGeomxSet", function(object) object@featureType)
setReplaceMethod("featureType", c("NanoStringGeomxSet", "character"),
                 function(object, value) {
                   if (value %in% c("Probe", "Target")) {
                       object@featureType <- value
                   } else {
                       stop("featureType must be 'Probe' or 'Target'")
                   }
                   return(object)
                 })

