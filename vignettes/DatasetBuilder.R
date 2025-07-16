library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)

datasetBuilder <- function(exprs_matrix, class) {
  exprs_transposed <- t(exprs_matrix) # transposing expression matrix, so that samples are rows and genes are columns
  combo_data <- cbind(exprs_transposed, class) # binding class column onto end of expression matrix
  combo_data <- data.frame(combo_data)
  return(combo_data)
}

