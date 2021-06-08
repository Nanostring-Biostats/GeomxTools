 ngeoMean <- function(v) {
    v[v == 0] <- 1
    return(EnvStats::geoMean(v, na.rm = TRUE))
}

ngeoSD <- function(v) {
  v[v == 0] <- 1
  return(geoSD(v, na.rm=T))
}

log10t <- function(x, thresh = 0.5) {
    if (min(x, na.rm = TRUE) < thresh) {
        x[!is.na(x) & x >= 0 & x < thresh] <- thresh
    }
    log2(x)
}

#' Add one to all counts in an expression matrix
#' 
#' @param object name of the NanoStringGeoMxSet object to perform QC on
#' @param elt expression matrix element in /code{assayDataElement}
#'        to shift all counts by
#' 
#' @examples
#' shiftOne(demoData)
#' 
#' @export
#' 
shiftOne <- function(object) {
    assay
}

#NEO rewrite in data.tables for speed
collapseCounts <- function(object) {
    probeCounts <- cbind(fData(object)[, c("TargetName", "Module")], 
        assayDataElement(object, elt="exprs"))
    collapsedCounts <- aggregate(formula=. ~ TargetName + Module, 
        data=probeCounts, FUN=ngeoMean)
    return(collapsedCounts)
}