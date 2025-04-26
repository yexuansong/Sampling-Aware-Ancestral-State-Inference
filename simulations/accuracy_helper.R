library(caret)  # For confusion matrix and other metrics
library(irr)

######################
# This function returns 4 accuracy measurements 

accuracy_helper <- function(asr,true_states,tip_states){
  
  reconstructed_states <- apply(asr, 1, which.max)
  comparison_data <- data.frame(
    true = true_states,
    reconstructed = c(tip_states, reconstructed_states)
  )
  comparison_data <- tail(comparison_data,length(tip_states)+1)
  
  # Overall accuracy
  accuracy <- mean(comparison_data$true == comparison_data$reconstructed)
  
  # Confusion matrix
  #conf_matrix <- confusionMatrix(factor(comparison_data$reconstructed), factor(comparison_data$true))
  
  # Cohen's Kappa
  #kappa <- kappa2(data.frame(true = comparison_data$true, reconstructed = comparison_data$reconstructed))$value
  
  prob_accuracy <- mean(sapply(1:nrow(asr), function(i) {
    asr[i, true_states[length(tip_states) + i]]
  }))
  return(list(
    acc=accuracy,
    #conf=conf_matrix,
    #ka=kappa,
    prob_acc=prob_accuracy
  ))
}