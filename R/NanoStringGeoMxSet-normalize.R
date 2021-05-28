#' normalize
#' @description normalize GeoMxSet using different normalization methods
#' @param object name of the object class to perform normalization on
#' @param norm_method the normalization method to be applied on the object
#' @param data_type the data type of the object. Values maybe RNA, protein. 
#' @param fromElt name of the assayDataElement to normalize
#' @param toElt name of the assayDataElement to store normalized values
#' @return a NanoStringGeoMxSet object with normalized counts and normalized factors
#' @export
#' @examples
#' norm_object <- normalize(demoData)
 
setMethod("normalize", "NanoStringGeoMxSet",
    function(object, norm_method=c("quant", "neg", "hk"), data_type=c("RNA", "protein"), 
             fromElt = "exprs", toElt = "exprs_norm", ...) 
        {
        norm_method <- match.arg(norm_method)
            switch(norm_method, "quant" = {quantileNorm(object, data_type=data_type)}, 
                "neg" = {negNorm(object, data_type=data_type)}, 
                "hk" = {hkNorm(object, data_type=data_type)})
        })

quantileNorm <- function(object, data_type, fromElt = "exprs" , toElt = "q_norm", desiredQuantile = .75) 
{ 
   #object<-checkQCFlags(object)
   ## Get quantile of counts for each sample
   qs <- apply(exprs(object), 2, function(x) quantile(x,desiredQuantile))
   ## Save the normfactors for desired quantile
   if (toElt != "q_norm")
       pData(object)[[paste(toElt, "qFactors", sep= "_")]] <- qs/ngeoMean(qs) 
   else
       pData(object)[["qnormFactors"]] <- qs/ngeoMean(qs) 
   assayDataElement(object, toElt, validate=TRUE) <- sweep(assayDataElement(object, fromElt), 2L, qs/ngeoMean(qs), FUN = "/")
   return(object)
} 

negNorm <- function(object, data_type, fromElt = "exprs" , toElt = "neg_norm") 
{ 
     if (!featureType(object) == "Target")
     {
         stop ("Negative Background normalization is for collapsed data set.  
            Run function aggregateCounts() to collapse the probes to targets.\n")
     } else
     {  
         negsubset <- subset(object, subset = CodeClass %in% c("Negative01", "Negative"))
         negs <- apply(exprs(negsubset), 2, function(x) ngeoMean(x))
         pData(object)[["negnormFactors"]] <- negs/ngeoMean(negs)
         assayDataElement(object, toElt) <- sweep(assayDataElement(object, fromElt), 2L, negs/ngeoMean(negs), FUN = "/")
         return(object)
     }
} 

hkNorm <- function(object, data_type, fromElt = "exprs" , toElt = "hk_norm") 
{ 
    housekeepers <- c('C1orf43','GPI','OAZ1','POLR2A','PSMB2','RAB7A',
                    'SDHA','SNRPD3','TBC1D10B','TPM4','TUBB','UBB')
    if (!featureType(object) == "Target")
    {
    stop ("Housekeeping normalization is for collapsed data set.  
            Run function aggregateCounts() to collapse the probes to targets.\n")
    } else
    {  
    hksubset <- subset(object, subset = TargetName %in% housekeepers)
    hks <- apply(exprs(hksubset), 2, function(x) ngeoMean(x))
    pData(object)[["hknormFactors"]] <- hks/ngeoMean(hks)
    assayDataElement(object, toElt) <- sweep(assayDataElement(object, fromElt), 2L, hks/ngeoMean(hks), FUN = "/")
    return(object)
    }
} 

#' Check QC Flags in the GeoMxSet and removes the probe or sample from the object
#' @rdname checkQCFlags
#' @param object name of the NanoStringGeoMxSet object to check the QC Flags
#' @param ...  for other arguments 
#' @return a NanoStringGeoMxSet object probes and samples failing QC removed
#' @export
#' @examples 
#' DemoData <- checkQCFlags(Demodata)
setGeneric("checkQCFlags", signature = c("object"), 
           function(object, ...)
               standardGeneric("checkQCFlags"))


#' checkQCFlags
#' @param NanoStringGeoMxSet 
#' @param removeLowLocalOutliers logical, if TRUE it sets outlier counts to zero,  default is FALSE,
#' @return NanoStringGeoMxSet
#' @export
#'
#' @examples
#' QCobject <- checkQCFlags(demoData)
setMethod("checkQCFlags", "NanoStringGeoMxSet",
          function(object, removeLowLocalOutliers = FALSE,  ...)
          {
              ## Remove all samples that failed AOI QC should it have been run
              AOIQCFlags <- protocolData(object)[["QCFlags"]]
              if (is.null(AOIQCFlags))
              {
                  warning("AOI QC has not been run on this data set.  Proceed with caution.\n")
              } else 
              {
                  QCResultsIndex <- which(apply(AOIQCFlags , 1L , function(x) sum(x) == 0L))
                  object <- object[, QCResultsIndex]
              }
              ## Remove all low probe count and ratio probes that failed QC  
              ProbeQCFlags <- fData(object)[["QCFlags"]]
              if (is.null(ProbeQCFlags))
              {
                  warning("Probe QC has not been run on this data set.  Proceed with caution.\n")
              } else 
              {
                  ProbeQCFlags <- ProbeQCFlags[, c("LowProbeCount", "LowProbeRatio", "GlobalOutlier")]
                  probeQCResultsIndex <- which(apply(ProbeQCFlags , 1L , function(x) sum(x) == 0L))
                  object <- object[probeQCResultsIndex, ]
              }
              ## Check if option to remove local outliers is set to TRUE
              if (removeLowLocalOutliers == TRUE)
              {
                  ProbeQCFlags <- fData(object)[["QCFlags"]]
                  ProbeQCFlags <- ProbeQCFlags[ , grepl( "HighLocalOutlier|LowLocalOutlier" , names(ProbeQCFlags))]
                  ## RV: This will remove all probes that has flags. Need to modify this to replace only the sample.
                  probeQCResultsIndex <- which(apply(ProbeQCFlags , 1L , function(x) x == TRUE))
                  object <- object[probeQCResultsIndex, ]
              }
            return(object)
          })