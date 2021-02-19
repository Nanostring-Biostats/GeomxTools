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

#NEO rewrite in data.tables for speed
collapseCounts <- function(object) {
    probeCounts <- cbind(fData(object)[, c("TargetName", "Module")], 
        assayDataElement(object, elt="exprs"))
    collapsedCounts <- aggregate(formula=. ~ TargetName + Module, 
        data=probeCounts, FUN=ngeoMean)
    return(collapsedCounts)
}