setMethod("normalize", "NanoStringGeomxSet",
    function(object, norm_method=c("quantile", "neg", "hk"), data_type="RNA", 
             fromElt = "exprs", toElt = "exprs_norm", ...) 
        {
        norm_method <- match.arg(norm_method)
            switch(norm_method, "quantile" = {quantileNorm(object, data_type=data_type, quantile)}, 
                "neg" = {neg_norm(object, data_type=data_type)}, 
                "hk" = {hkNorm(object, data_type=data_type)})
        })

quantileNorm <- function(object, data_type, fromElt = "exprs" , toElt = "q_norm") 
{ 
   #object<-checkQCFlags(object)
   ## Get quantile of counts for each sample
   quantiles <- apply(exprs(object), 2, function(x) quantile(x,.75))
   pData(object)[["qnormFactors"]] <- quantiles/ngeoMean(quantiles) 
   assayDataElement(object, toElt) <- sweep(assayDataElement(object, fromElt), 2L, quantiles/ngeoMean(quantiles), FUN = "*")
   return(object)
} 
 
 negNorm <- function(object, data_type, fromElt = "exprs" , toElt = "exprs_norm") 
{ 
} 

hkNorm <- function(object, data_type) 
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
    assayDataElement(object, toElt) <- sweep(assayDataElement(object, fromElt), 2L, hks/ngeoMean(hks), FUN = "*")
    return(object)
    }
} 

setGeneric("checkQCFlags", signature = c("object"), 
           function(object, ...)
               standardGeneric("checkQCFlags"))

setMethod("checkQCFlags", "NanoStringGeomxSet",
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


