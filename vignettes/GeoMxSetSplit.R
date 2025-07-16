library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(dplyr)

splitgeomx <- function(geomx_set_object, split_ratio) {
  
  dimslice <- dim(pData(geomx_set_object))[1] * split_ratio # determines how many samples will be in the training data
  sampled <- slice_sample(pData(protocolData(geomx_set_object)), n = dimslice) # sampling without replacement
  sampled_ids <- rownames(sampled) #getting rownames (which are sample ids)

  split_train <- geomx_set_object[, c(sampled_ids)] # everything in sampled_ids is in the training set
  split_test <- geomx_set_object[, -which(colnames(geomx_set_object) %in% sampled_ids)] # everything no included in sampled_ids is in the testing set
  
  return_list <- c(split_train, split_test) # list containing the train and test GeoMxSet objects
  
  return(return_list)
}
