ngeoMean <- function(v) {
    v[v == 0] <- 1
    return(EnvStats::geoMean(v, na.rm = TRUE))
}

collapseCounts <- function(object) {
    probeCounts <- cbind(fData(object)[, c("Target", "Module")], 
        assayDataElement(object, elt="exprs"))
    collapsedCounts <- aggregate(formula=. ~ Target + Module, 
        data=probeCounts, FUN=ngeoMean)
    return(collapsedCounts)
}