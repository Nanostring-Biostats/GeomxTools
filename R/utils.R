ngeoMean <- function(v) {
  v[v == 0] <- 1
  return(EnvStats::geoMean(v, na.rm = TRUE))
}
