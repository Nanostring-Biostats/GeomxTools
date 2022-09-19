library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(dplyr)
library(ggforce)
library(scales)
library(cowplot)
library(caTools)
library(data.table)
library(caret)
library(e1071)

splitgeomx <- function(geomx_set_object, split_ratio) {
  
  dimslice <- dim(pData(geomx_set_object))[1] * split_ratio
  sampled <- slice_sample(pData(protocolData(geomx_set_object)), n = dimslice)
  sampled_ids <- rownames(sampled)

  split_train <- geomx_set_object[, c(sampled_ids)]
  split_test <- geomx_set_object[, -which(colnames(geomx_set_object) %in% sampled_ids)]
  
  return_list <- c(split_train, split_test)
  
  return(return_list)
}