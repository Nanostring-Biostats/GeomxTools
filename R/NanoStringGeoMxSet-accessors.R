# sData Accessor
setGeneric("sData", signature = "object",
           function(object) standardGeneric("sData"))
setMethod("sData", "NanoStringGeoMxSet",
          function(object) cbind(pData(object), pData(protocolData(object))))

# svarLabels Accessor
setGeneric("svarLabels", signature = "object",
           function(object) standardGeneric("svarLabels"))
setMethod("svarLabels", "NanoStringGeoMxSet",
          function(object) c(varLabels(object), varLabels(protocolData(object))))

# dimLabels Accessor and Replacer
setMethod("dimLabels", "NanoStringGeoMxSet", function(object) object@dimLabels)
setReplaceMethod("dimLabels", c("NanoStringGeoMxSet", "character"),
                 function(object, value) {
                   object@dimLabels <- value
                   return(object)
                 })

# design Accessor and Replacer
setMethod("design", "NanoStringGeoMxSet", function(object) object@design)
setReplaceMethod("design", c("NanoStringGeoMxSet", "formula"),
                 function(object, value) {
                   object@design <- value
                   return(object)
                 })
setReplaceMethod("design", c("NanoStringGeoMxSet", "ANY"),
                 function(object, value) {
                   object@design <- stats::as.formula(value)
                   return(object)
                 })
setReplaceMethod("design", c("NanoStringGeoMxSet", "NULL"),
                 function(object, value) {
                   object@design <- NULL
                   return(object)
                 })

setGeneric("featureType", signature = "object",
           function(object) standardGeneric("featureType"))
setMethod("featureType", "NanoStringGeoMxSet", function(object) object@featureType)
setGeneric("featureType<-", signature = c("object", "value"), 
           function(object, value) standardGeneric("featureType<-"))
setReplaceMethod("featureType", c("NanoStringGeoMxSet", "character"),
                 function(object, value) {
                   if (value %in% c("Probe", "Target")) {
                       object@featureType <- value
                   } else {
                       stop("featureType must be 'Probe' or 'Target'")
                   }
                   return(object)
                 })

