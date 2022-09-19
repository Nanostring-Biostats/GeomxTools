classifier_metrics <- function(confusion_matrix_table) {
  TP <- confusion_matrix_table[1:1, 1:1]
  FP <- confusion_matrix_table[1:1, 2:2]
  TN <- confusion_matrix_table[2:2, 2:2]
  FN <- confusion_matrix_table[2:2, 1:1]
  
  accuracy <- (TP+TN)/(TP+TN+FP+FN)
  accuracy
  
  precision <- TP/(TP+FP)
  precision
  
  recall <- TP/(TP+FN)
  recall
  
  specificity <- TN/(TN+FP)
  specificity
  
  metrics_list <- c("accuarcy" = accuracy, "precision" = precision, "recall" = recall, "specificity" = specificity)
  
  return(metrics_list)
}