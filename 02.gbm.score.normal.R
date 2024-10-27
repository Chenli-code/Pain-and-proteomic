#=====  Model protein =====
library(caret)
library(foreach)
library(doParallel)
library(dplyr)
library(rms)
library(ROCit)
library(cutpointr)
library(xgboost)
library(reticulate)
# library(tensorflow)
library(caret)
library(pROC)
##### protein##### 
library(glmnet)
library(tidyverse)
library(tidymodels)
library(recipes)
library(parsnip)
# library(workflows)
library(tune)
library(gbm)
#install.packages("shapviz")
library(shapviz)
#install.packages("xgboost")
library(xgboost)


# load data
protein_selected = readRDS('name_protein.rds')
load('Data_analysis_R2.RData')
data = data %>% filter(!is.na(missing_per_protein))

##### normalization for protein data #####
for (pro in protein_selected) {
  cat(pro,'-')
  data[[pro]] = scale(data[[pro]])
}

outcomes = c('chronic_pain_status',
             'Head_and_Face_Pain',
             'Abdominal_Pain','Generalized_Pain',
             'Neck_Shoulder_and_Back_Hip_and_Knee_Pain')

# Initialize an empty dataframe to store variable importance results
auc_df = data.frame(outcome = NA,test=NA,train=NA,exposure =NA)
importance_df = data.frame(outcome = NA, pro = NA, importance = NA, exposure = NA)
model_list_p = list()


Protein_association = read.csv('Results_wide_association/Protein_association-chronic.csv')
Protein_association = Protein_association %>% filter(p_value < 0.05/2923)


##### Protein #####
for (outcome in outcomes) {
  
  cat(outcome, ' - ')
  ##### all ProRS #####
  # Select the proteins associated with the current outcome
  lasso_protein = Protein_association[which(Protein_association$outcome == outcome), 'exposure']
  data_protein = data
  data_protein[, outcome] = ifelse(data_protein[, outcome] == 1,'Y','N')
  data_protein[, outcome] = as.factor(data_protein[, outcome])
  set.seed(998)
  inTraining <- createDataPartition(data_protein[, outcome], p = 0.8, list = FALSE)
  train <- data_protein[ inTraining,]
  test  <- data_protein[-inTraining,]
  model_list_p[[paste0(outcome,'_train')]] = train
  model_list_p[[paste0(outcome,'_test')]] = test
  
  # Train a ranger model with hyperparameter tuning (excluding num.trees)
  rf_model <- caret::train(
    as.formula(paste0(outcome, '~.')),
    data = na.omit(train[, c(outcome, lasso_protein)]),
    method = 'gbm',
    trControl = trainControl(method = "cv",
                             number = 10,  # 10-fold cross-validation
                             verboseIter = FALSE,
                             classProbs = TRUE)
  )
  # rf_model = model_list_p[[outcome]]
  model_list_p[[outcome]] = rf_model
  # Calculate variable importance
  importance <- varImp(rf_model, scale = TRUE)$importance
  importance = as.data.frame(importance)
  colnames(importance) = 'importance'
  importance$pro = rownames(importance)
  importance$outcome = outcome
  importance$exposure = 'Protein'
  
  # Append the importance results to the dataframe
  importance_df = rbind(importance_df, importance)
  
  # Predict risk scores in test data
  matching_ind <- which(complete.cases(test[, c(lasso_protein)]))
  ProRS_test = predict(rf_model, test[matching_ind, c(lasso_protein)], type = 'prob')[,'Y']
  matching_indices <- which(complete.cases(train[, c(lasso_protein)]))
  ProRS = predict(rf_model, train[matching_indices, c(lasso_protein)], type = 'prob')
  
  auc = data.frame(outcome = outcome,
                   test = rocit(ProRS_test,test[matching_ind,outcome])$AUC,
                   train = rocit(ProRS[, 'Y'],train[matching_indices,outcome])$AUC,
                   exposure = 'all')
  auc_df = rbind(auc_df,auc)
  # Predict risk scores
  ProRS = predict(rf_model, data[, c(lasso_protein)], type = 'prob')
  
  # Get indices of rows without NAs for consistent row numbers
  matching_indices <- which(complete.cases(data[, c(lasso_protein)]))
  
  # Update the original dataframe with the risk scores
  data[matching_indices, paste0(outcome, '_ProRS')] = ProRS[, 'Y']
  
  
  
  ##### simple ProRS #####
  top10 = importance %>% arrange(-importance) 
  top10 = top10[1:10,'pro']
  # Train a ranger model with hyperparameter tuning (excluding num.trees)
  rf_model_s <- caret::train(
    as.formula(paste0(outcome, '~.')), 
    data = na.omit(train[, c(outcome, top10)]),
    method = 'gbm',
    trControl = trainControl(method = "cv", 
                             number = 10,  # 10-fold cross-validation
                             verboseIter = FALSE,
                             classProbs = TRUE)
  )
  # Predict risk scores in test data
  matching_ind <- which(complete.cases(test[, c(top10)]))
  ProRS_test = predict(rf_model_s, test[matching_ind, c(top10)], type = 'prob')[,'Y']
  matching_indices <- which(complete.cases(train[, c(top10)]))
  ProRS = predict(rf_model_s, train[matching_indices, c(top10)], type = 'prob')
  
  
  auc = data.frame(outcome = outcome,
                   test = rocit(ProRS_test,test[matching_ind,outcome])$AUC,
                   train = rocit(ProRS[, 'Y'],train[matching_indices,outcome])$AUC,
                   exposure = 'simple')
  auc_df = rbind(auc_df,auc)
  # Predict risk scores
  ProRS = predict(rf_model_s, data[, c(top10)], type = 'prob')
  
  # Get indices of rows without NAs for consistent row numbers
  matching_indices <- which(complete.cases(data[, c(top10)]))
  
  # Update the original dataframe with the risk scores
  data[matching_indices, paste0(outcome, '_simpleRS')] = ProRS[, 'Y']
  
}


ProRS = data[, c('f.eid',paste0(outcomes, '_ProRS'),paste0(outcomes, '_simpleRS'))]


# Save the variable importance results to a CSV file
write.csv(importance_df, file = 'Results_wide_association/relative_df_normalization_chronic.csv', row.names = FALSE)

save(ProRS, file = 'ProRS_gbm.RData')
