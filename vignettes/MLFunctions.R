# Loading package
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(data.table)
library(ModelMetrics)
library(caret)

ml_function <- function(geomx_set_object, classifier_column, model_type, split_ratio=0.7, ...) {
  
  source("~/Downloads/GeomxTools-nayana_machine_learning/vignettes/GeoMxSetSplit.R")
  
  train_test_list <- splitgeomx(geomx_set_object, split_ratio)
  train_split <- train_test_list[[1]]
  test_split <- train_test_list[[2]]
  
  source("~/Downloads/GeomxTools-nayana_machine_learning/vignettes/DatasetBuilder.R")
  
  train_data <- datasetBuilder(exprs(train_split), pData(train_split)[[classifier_column]])
  test_data <- datasetBuilder(exprs(test_split), pData(test_split)[[classifier_column]])
  
  switch(model_type,
         "nb" = custom_nb(train_split, test_split, classifier_column, ...),
         "knn" = custom_knn(train_split, test_split, classifier_column, ...),
         "rf" = custom_rf(train_split, test_split, classifier_column, ...),
         "svm" = custom_svm(train_split, test_split, classifier_column, ...))
}

# Naive Bayes Input Parameters:
# tuneLength: amount of granularity in the tuning parameter grid
# trControl: specifies training resampling method
# tuneGrid: data frame with tuning values
custom_nb <- function(train_split, test_split, classifier_column, tuneLength=2, trControl=caret::trainControl(method="cv"), tuneGrid=NULL, ...) {

  library(naivebayes)
  
  nb_model <- caret::train(t(exprs(train_split)),
                    as.factor(pData(train_split)[[classifier_column]]),
                    method = "naive_bayes",
                    tuneLength=tuneLength,
                    trControl=trControl,
                    tuneGrid=tuneGrid, ...)
  
  print(nb_model)
  
  imp_nb <- caret::varImp(nb_model, useModel=FALSE, nonpara=FALSE, scale=FALSE)
  print(imp_nb)
  
  nb_predict <- predict(nb_model, newdata = t(exprs(test_split)))
  nb_predict_prob <- predict(nb_model, newdata = t(exprs(test_split)), type = "prob")[,2]
  
  cm_nb <- caret::confusionMatrix(nb_predict, as.factor(pData(test_split)[[classifier_column]]))
  print(cm_nb)
  
  library(ROCR)
  pred <- prediction(nb_predict_prob, as.factor(pData(test_split)[[classifier_column]]))
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
  
  auc <- performance(pred, "auc")
  auc <- auc@y.values[[1]]
  print("AUC")
  print(auc)
}

# KNN Input Parameters:
# tuneLength: amount of granularity in the tuning parameter grid
# trControl: specifies training resampling method
# tuneGrid: data frame with tuning values
          # k: number of neighbors
custom_knn <- function(train_split, test_split, classifier_column, tuneLength=8, trControl=caret::trainControl(method="cv"), tuneGrid=expand.grid(k=3:10), ...) {
  
  library(caret)
  
  knn_model <- caret::train(t(exprs(train_split)), 
                     as.factor(pData(train_split)[[classifier_column]]), 
                     method = "knn",
                     tuneLength=tuneLength,
                     trControl=trControl, 
                     tuneGrid=tuneGrid, ...)
  
  print(knn_model)
  
  imp_knn <- caret::varImp(knn_model, useModel=FALSE, nonpara=FALSE, scale=FALSE)
  print(imp_knn)
  
  knn_predict <- predict(knn_model, newdata = t(exprs(test_split)))
  knn_predict_prob <- predict(knn_model, newdata = t(exprs(test_split)), type = "prob")[,2]
  
  cm_knn <- caret::confusionMatrix(knn_predict, as.factor(pData(test_split)[[classifier_column]]))
  print(cm_knn)
  
  library(ROCR)
  pred <- prediction(knn_predict_prob, as.factor(pData(test_split)[[classifier_column]]))
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
  
  auc <- performance(pred, "auc")
  auc <- auc@y.values[[1]]
  print("AUC")
  print(auc)
}

# Random Forest Input Parameters:
# tuneLength: amount of granularity in the tuning parameter grid
# trControl: specifies training resampling method
# tuneGrid: data frame with tuning values
custom_rf <- function(train_split, test_split, classifier_column, tuneLength=5, trControl=caret::trainControl(method="cv"), tuneGrid=NULL, ...) {
  
  library(randomForest)
  
  rf_model <- caret::train(t(exprs(train_split)), 
                            as.factor(pData(train_split)[[classifier_column]]), 
                            method = "rf",
                            tuneLength=tuneLength,
                            trControl=trControl, 
                            tuneGrid=tuneGrid, ...)
  
  print(rf_model)
  
  imp_rf <- caret::varImp(rf_model, useModel=FALSE, nonpara=FALSE, scale=FALSE)
  print(imp_rf)
  
  rf_predict <- predict(rf_model, newdata = t(exprs(test_split)))
  rf_predict_prob <- predict(rf_model, newdata = t(exprs(test_split)), type = "prob")[,2]
  
  cm_rf <- caret::confusionMatrix(rf_predict, as.factor(pData(test_split)[[classifier_column]]))
  print(cm_rf)
  
  library(ROCR)
  pred <- prediction(rf_predict_prob, as.factor(pData(test_split)[[classifier_column]]))
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
  
  auc <- performance(pred, "auc")
  auc <- auc@y.values[[1]]
  print("AUC")
  print(auc)
}

# SVM Input Parameters:
# tuneLength: amount of granularity in the tuning parameter grid
# trControl: specifies training resampling method
# tuneGrid: data frame with tuning values
custom_svm <- function(train_split, test_split, classifier_column, tuneLength=3, trControl=caret::trainControl(method="cv"), tuneGrid=NULL, ...) {
  
  library(LiblineaR)
  
  
  svm_model <- caret::train(t(exprs(train_split)),
                     as.factor(pData(train_split)[[classifier_column]]),
                     method = "svmLinearWeights2",
                     tuneLength=tuneLength,
                     trControl=trControl,
                     tuneGrid=tuneGrid)
  
  print(svm_model)
  
  imp_svm <- caret::varImp(svm_model, useModel=FALSE, nonpara=FALSE, scale=FALSE)
  print(imp_svm)
  
  svm_predict <- predict(svm_model, newdata = t(exprs(test_split)))
  
  cm_svm <- caret::confusionMatrix(svm_predict, as.factor(pData(test_split)[[classifier_column]]))
  print(cm_svm)
}
