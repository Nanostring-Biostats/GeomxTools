ngeoMean <- function(v) {
    v[v == 0] <- 1
    return(EnvStats::geoMean(v, na.rm = TRUE))
}

#NEO rewrite in data.tables for speed
collapseCounts <- function(object) {
    probeCounts <- cbind(fData(object)[, c("Target", "Module")], 
        assayDataElement(object, elt="exprs"))
    collapsedCounts <- aggregate(formula=. ~ Target + Module, 
        data=probeCounts, FUN=ngeoMean)
    return(collapsedCounts)
}