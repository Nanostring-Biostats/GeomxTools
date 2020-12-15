ngeoMean <- function(v) {
  v[v == 0] <- 1
  return(geoMean(v, na.rm = TRUE))
}
