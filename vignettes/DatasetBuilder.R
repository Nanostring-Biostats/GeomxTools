library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(data.table)

datasetBuilder <- function(exprs_matrix, class) {
  exprs_transposed <- t(exprs_matrix)
  combo_data <- cbind(exprs_transposed, class)
  combo_data <- data.frame(combo_data)
  
  return(combo_data)
}

