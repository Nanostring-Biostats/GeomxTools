# Loading package
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(data.table)
library(ModelMetrics)

ml_function <- function(geomx_set_object, classifier_column, model_type, split_ratio=0.7, ...) {
  
  source("~/hades/nyeshlur/GeoMxSetSplit.R")
  
  train_test_list <- splitgeomx(geomx_set_object, split_ratio)
  train_split <- train_test_list[[1]]
  test_split <- train_test_list[[2]]
  
  source("~/hades/nyeshlur/DatasetBuilder.R")
  
  train_data <- datasetBuilder(exprs(train_split), pData(train_split)[[classifier_column]])
  test_data <- datasetBuilder(exprs(test_split), pData(test_split)[[classifier_column]])
  
  switch(model_type,
         "nb" = custom_nb(train_split, test_split, classifier_column, ...),
         "knn" = custom_knn(train_split, test_split, classifier_column, ...),
         "xgb" = custom_xgb(train_split, test_split, classifier_column, ...),
         "rf" = custom_rf(train_split, test_split, classifier_column, ...),
         "svm" = custom_svm(train_split, test_split, classifier_column, ...))
}

custom_nb <- function(train_split, test_split, classifier_column, seed=50, tuneLength=2, trControl=trainControl(method="repeatedcv", repeats=3), tuneGrid=NULL, ...) {

  library(naivebayes)
  
  set.seed(seed)
  
  nb_model <- train(t(exprs(train_split)),
                    as.factor(pData(train_split)[[classifier_column]]),
                    method = "naive_bayes",
                    tuneLength=tuneLength,
                    trControl=trControl,
                    tuneGrid=tuneGrid)
  
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

custom_knn <- function(train_split, test_split, classifier_column, seed=50, tuneLength=NULL, trControl=trainControl(method="repeatedcv", repeats=3), tuneGrid=expand.grid(k=3:10), ...) {
  
  library(caret)
  
  set.seed(seed)
  
  knn_model <- train(t(exprs(train_split)), 
                     as.factor(pData(train_split)[[classifier_column]]), 
                     method = "knn",
                     tuneLength=tuneLength,
                     trControl=trControl, 
                     tuneGrid=tuneGrid)
  
  print(knn_model)
  plot(knn_model)
  
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

custom_xgb <- function(train_split, test_split, classifier_column, seed=50, tuneLength=1, trControl=trainControl(method="none"), tuneGrid=expand.grid(nrounds=100, max_depth=6, eta=0.3, gamma=0, subsample=1, colsample_bytree=1, rate_drop=0, skip_drop=0, min_child_weight=1), ...) {
  
  library(xgboost)
  
  set.seed(seed)
  
  xgb_model <- train(t(exprs(train_split)), 
                     as.factor(pData(train_split)[[classifier_column]]), 
                     method = "xgbDART", 
                     tuneLength=tuneLength,
                     trControl=trControl,
                     tuneGrid=tuneGrid)
  
  print(xgb_model)
  
  imp_xgb <- caret::varImp(xgb_model, useModel=FALSE, nonpara=FALSE, scale=FALSE)
  print(imp_xgb)
  
  xgb_predict <- predict(xgb_model, newdata = t(exprs(test_split)))
  xgb_predict_prob <- predict(xgb_model, newdata = t(exprs(test_split)), type = "prob")[,2]
  
  cm_xgb <- caret::confusionMatrix(xgb_predict, as.factor(pData(test_split)[[classifier_column]]))
  print(cm_xgb)
  
  library(ROCR)
  pred <- prediction(xgb_predict_prob, as.factor(pData(test_split)[[classifier_column]]))
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
  
  auc <- performance(pred, "auc")
  auc <- auc@y.values[[1]]
  print("AUC")
  print(auc)
}

custom_rf <- function(train_split, test_split, classifier_column, seed=50, ntreeTry=100, stepFactor=5, improve=0.05, trace=FALSE, plot=TRUE, doBest=TRUE, ...) {
  
  library(randomForest)
  
  set.seed(seed)
  
  random_forest_model <- tuneRF(t(exprs(train_split)), as.factor(pData(train_split)[[classifier_column]]), ntreeTry=ntreeTry, stepFactor=stepFactor, improve=improve, trace=trace, plot=plot, doBest=doBest, ...)
  
  plot(random_forest_model)
  
  print(random_forest_model)
  
  importance(random_forest_model)
  varImpPlot(random_forest_model, sort = TRUE, n.var = 10, scale = TRUE)
  
  rf_prediction <- predict(random_forest_model, t(exprs(test_split)), type = "response", norm.votes = FALSE)
  rf_prediction_prob <- predict(random_forest_model, t(exprs(test_split)), type = "prob", norm.votes = FALSE)[,2]
  
  cm_rf <- caret::confusionMatrix(rf_prediction, as.factor(pData(test_split)[[classifier_column]]))
  print(cm_rf)
  
  library(ROCR)
  pred <- prediction(rf_prediction_prob, as.factor(pData(test_split)[[classifier_column]]))
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
  
  auc <- performance(pred, "auc")
  auc <- auc@y.values[[1]]
  print("AUC")
  print(auc)
}

custom_svm <- function(train_split, test_split, classifier_column, seed=100, tuneLength=3, trControl=trainControl(method="repeatedcv", repeats=3), tuneGrid=NULL, ...) {
  
  library(LiblineaR)
  
  set.seed(seed)
  
  svm_model <- train(t(exprs(train_split)),
                     as.factor(pData(train_split)[[classifier_column]]),
                     method = "svmLinearWeights2",
                     tuneLength=tuneLength,
                     trControl=trControl,
                     tuneGrid=tuneGrid)
  
  print(svm_model)
  
  plot(svm_model)
  
  imp_svm <- caret::varImp(svm_model, useModel=FALSE, nonpara=FALSE, scale=FALSE)
  print(imp_svm)
  
  svm_predict <- predict(svm_model, newdata = t(exprs(test_split)))
  
  cm_svm <- caret::confusionMatrix(svm_predict, as.factor(pData(test_split)[[classifier_column]]))
  print(cm_svm)
}

