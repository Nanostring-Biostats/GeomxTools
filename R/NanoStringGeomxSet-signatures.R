# signatures Accessor and Replacer
setGeneric("signatures", signature = "object",
           function(object) standardGeneric("signatures"))
setMethod("signatures", "NanoStringGeomxSet",
          function(object) object@signatures)

setGeneric("signatures<-", signature = c("object", "value"),
           function(object, value) standardGeneric("signatures<-"))
setReplaceMethod("signatures", c("NanoStringGeomxSet", "SignatureSet"),
                 function(object, value) {
                   object@signatures <- value
                   return(object)
                 })

setGeneric("signatureScores", signature = "object",
           function(object, ...) standardGeneric("signatureScores"))
setMethod("signatureScores", "NanoStringGeomxSet",
          function(object, elt = "exprs") {
            if (length(signatures(object)) == 0L) {
              return(matrix(numeric(), nrow = 0L, ncol = ncol(object),
                            dimnames = list(NULL, colnames(object))))
            }
            exprs <- t(assayDataElement2(object, elt))
            colnames(exprs) <- featureData(object)[["TargetName"]]
            sigFuncList <- signatureFuncs( object )
            linWeights <- stats::weights( signatures( object ) )[names( sigFuncList )[which( sigFuncList %in% "default" )]]
            nonLinFuncs <- sigFuncList[which( !( sigFuncList %in% "default" ) )]
            scores <- .sigCalc(exprs, linWeights)
            while (length(idx <- which(rowSums(is.na(scores)) > 0L))) {
              subscores <-
                .sigCalc(cbind(exprs, t(scores[-idx, , drop = FALSE])),
                         stats::weights(signatures(object))[idx])
              if (all(is.na(subscores)))
                break
              else
                scores[idx, ] <- subscores
            }
            nonLinScores <- t( sapply( nonLinFuncs , function( x , elt ) eval( parse( text = paste( x , "( object , fromElt = \"" , elt , "\" )" , sep = "" ) ) ) , elt ) )
            if (ncol(nonLinScores) > 0) {
                scores <- rbind( scores , nonLinScores )
            }
            return( scores[names( stats::weights( signatures( object ) ) ),] )
          })

setGeneric( "signatureGroups" , signature = "object" ,
            function (object , ... ) standardGeneric( "signatureGroups" ) )
setMethod("signatureGroups", "NanoStringGeomxSet",
          function( object ) {
            return( groups( object@signatures ) )
          } )

setGeneric("setSignatureGroups<-", signature = c("object", "value"),
           function(object, value) standardGeneric("setSignatureGroups<-"))
setReplaceMethod("setSignatureGroups", c("NanoStringGeomxSet", "factor"),
                 function(object, value) {
                   groups( object@signatures ) <- value
                   return( object )
                 })
setReplaceMethod("setSignatureGroups", c("NanoStringGeomxSet", "character"),
                 function(object, value) {
                   groups( object@signatures ) <- value
                   return( object )
                 })

setGeneric( "signatureFuncs" , signature = "object" ,
            function (object , ... ) standardGeneric( "signatureFuncs" ) )
setMethod("signatureFuncs", "NanoStringGeomxSet",
          function( object ) {
            return( getSigFuncs( object@signatures ) )
          } )

setGeneric("setSignatureFuncs<-", signature = c("object", "value"),
           function(object, value) standardGeneric("setSignatureFuncs<-"))
setReplaceMethod("setSignatureFuncs", c("NanoStringGeomxSet", "character"),
                 function(object, value) {
                   setSigFuncs( object@signatures ) <- value
                   return( object )
                 })
