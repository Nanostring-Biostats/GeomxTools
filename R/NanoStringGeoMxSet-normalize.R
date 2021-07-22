HOUSEKEEPERS <- c('C1orf43','GPI','OAZ1','POLR2A','PSMB2','RAB7A',
                                    'SDHA','SNRPD3','TBC1D10B','TPM4','TUBB','UBB')
#' normalize
#' @description normalize GeoMxSet using different normalization methods
#' @param object name of the object class to perform normalization on
#' @param norm_method the normalization method to be applied on the object
#' @param data_type the data type of the object. Values maybe RNA, protein.
#' @param fromElt name of the assayDataElement to normalize
#' @param toElt name of the assayDataElement to store normalized values
#' @return a NanoStringGeoMxSet object with normalized counts and normalized factors
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'                        package="GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' norm_object <- normalize(demoData)
#' @export

setMethod("normalize", "NanoStringGeoMxSet",
    function(object, norm_method=c("quant", "neg", "hk", "subtractBackground"),
             data_type=c("RNA", "protein"), fromElt = "exprs", toElt = "exprs_norm",
             housekeepers = HOUSEKEEPERS, ...)
        {
        norm_method <- match.arg(norm_method)
            switch(norm_method,
                "quant" = {quantileNorm(object, data_type=data_type,
                                        toElt = toElt , fromElt = fromElt, ...)},
                "neg" = {negNorm(object, data_type=data_type,
                                 toElt = toElt, fromElt = fromElt, ...)},
                "hk" = {hkNorm(object, data_type=data_type,
                               toElt = toElt, fromElt = fromElt,
                               housekeepers = housekeepers, ...)},
                "subtractBackground" = {subtractBackground(object, data_type=data_type,
                                                           toElt = toElt,
                                                           fromElt = fromElt,...)}
                )
        })

quantileNorm <- function(object, data_type, desiredQuantile = .75, toElt, fromElt)
{
   ## Get quantile of counts for each sample
   qs <- apply(exprs(object), 2, function(x) stats::quantile(x,desiredQuantile))
   ## Save the normfactors for desired quantile
   if (toElt != "exprs_norm")
       pData(object)[[paste(toElt, "qFactors", sep= "_")]] <- qs/ngeoMean(qs)
   else
       pData(object)[["normFactors"]] <- qs/ngeoMean(qs)
   assayDataElement(object, toElt, validate=TRUE) <- sweep(assayDataElement(object, fromElt), 2L, qs/ngeoMean(qs), FUN = "/")
   return(object)
}

negNorm <- function(object, data_type, toElt, fromElt)
{
    if (!featureType(object) == "Target")
    {
        stop ("Error: Negative Background normalization is for collapsed data set.
        Run function aggregateCounts() to collapse the probes to targets.\n")
    }
    if (is.null(fData(object)[["Module"]]))
    {
        stop ("Error: Module is not specified in the object. Check your GeoMxSet object. \n")
    }

    # check if single panel
    pools <- as.list(unique(fData(object)[["Module"]]))
    if (length(pools) == 1){
        # compute negative norm factors for samples and divide by geomean
        negsubset <- subset(object, subset = CodeClass == "Negative")
        negs <- apply(exprs(negsubset), 2, function(x) ngeoMean(x))

        ## Save the normfactors in desired pData element
        if (toElt != "exprs_norm")
            pData(object)[[paste(toElt, "negFactors", sep= "_")]] <- negs/ngeoMean(negs)
        else
            pData(object)[["normFactors"]] <- negs/ngeoMean(negs)
        assayDataElement(object, toElt) <- sweep(assayDataElement(object, fromElt), 2L, negs/ngeoMean(negs), FUN = "/")
        return(object)
    }
    # for multipanel
    if (length(pools) > 1){
        lapply(pools, function (x){
            pool_neg <- subset(object, subset = CodeClass == "Negative" & Module == x)
            # check if there is only 1 negative probe per panel
            if (length(fData(pool_neg)[["TargetName"]]) < 1) {
                stop(paste0("Error: No negative probes could be located for panel ",
                            x, "."))}
            if (length(fData(pool_neg)[["TargetName"]]) > 1) {
            stop(paste0("Error: More than one negative was located for panel",
                        x, "."))}
        })
        ## Compute neg factors for each subset by module and use for normalization
        multiNegNorm <- lapply(pools, function (x){
            subsetModule <- subset(object, subset = Module == x)
            # compute negative norm factors for samples and divide by geomean
            negs <- as.numeric(exprs(negativeControlSubset(subsetModule)))
            normFactors <- negs/ngeoMean(negs)
            norm_neg <- function(x){
                x <- x/(negs/geoMean(negs))
            }
            norm_exprs <- t(assayDataApply(subsetModule, 1, norm_neg))
            norm_list <- list(normFactors = normFactors, norm_exprs = norm_exprs)
            return(norm_list)
        })
        names(multiNegNorm) <- pools
        ## Save the normfactors in desired pData element
        if (toElt != "exprs_norm")
            pData(object)[[paste(toElt, "negFactors", sep= "_")]] <- do.call(cbind, lapply(multiNegNorm, '[[', 1))
        else
            pData(object)[["normFactors"]] <- do.call(cbind, lapply(multiNegNorm, '[[', 1))
        ## Save the exprs_norm in desired pData element
        assayDataElement(object, toElt) <- do.call(rbind, lapply(multiNegNorm, '[[', 2))
        return(object)
    }
}

hkNorm <- function(object, data_type, toElt, fromElt, housekeepers)
{
    if (!featureType(object) == "Target")
    {
    stop ("Housekeeping normalization is for collapsed data set.
            Run function aggregateCounts() to collapse the probes to targets.\n")
    } else
    {
    hksubset <- subset(object, subset = TargetName %in% housekeepers)
    hks <- apply(exprs(hksubset), 2, function(x) ngeoMean(x))
    ## Save the normfactors in desired pData element
    if (toElt != "exprs_norm")
        pData(object)[[paste(toElt, "hkFactors", sep= "_")]] <- hks/ngeoMean(hks)
    else
        pData(object)[["hknormFactors"]] <- hks/ngeoMean(hks)
    assayDataElement(object, toElt) <- sweep(assayDataElement(object, fromElt), 2L, hks/ngeoMean(hks), FUN = "/")
    return(object)
    }
}


#subtract background
subtractBackground <- function(object, data_type, toElt, fromElt)
{
    if (!featureType(object) == "Target")
    {
        negsubset <- subset(object, subset = Codeclass %in% c("Negative01", "Negative"))
        negs <- apply(exprs(negsubset), 2, function(x) ngeoMean(x))
        assayDataElement( object, toElt ) <-
            t(assayDataApply( object, MARGIN=1L, FUN=`-`, t(negs) , elt=fromElt ))
    } else
    assayDataElement( object, toElt ) <-
    t(assayDataApply( object, MARGIN=1L, FUN=`-`, t(exprs(object)["Negative Probe", ]) , elt=fromElt ))
    return(object)
}

#' Check QC Flags in the GeoMxSet and removes the probe or sample from the object
#' @rdname checkQCFlags
#' @param object name of the NanoStringGeoMxSet object to check the QC Flags
#' @param ...  for other arguments
#' @return a NanoStringGeoMxSet object probes and samples failing QC removed
#' @export
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'                        package="GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
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
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'                        package="GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
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
