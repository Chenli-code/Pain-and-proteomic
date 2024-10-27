#===== xgboost Model protein =====
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
library(ggpubr)

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

####
Protein_association = read.csv('Results_wide_association/Protein_association-chronic.csv')
Protein_association = Protein_association %>% filter(p_value < 0.05/2923)

assay = readRDS('protein_assay.rds')
assay = assay %>% mutate(pro = paste0('P',coding))



auc_df = data.frame()
importance_df = data.frame()
model_list_p = list()

df = data

##### Protein #####
for (outcome in outcomes) {
  
  cat(outcome, ' - ')
  # Select the proteins associated with the current outcome
  lasso_protein <- Protein_association[which(Protein_association$outcome == outcome), 'exposure']

  
  X <- df[, lasso_protein]
  y <- df[, outcome]
  non_na_indices <- complete.cases(X, y)
  X <- X[non_na_indices, ]
  y <- y[non_na_indices]
  
  ctrl <- trainControl(method = "cv",  
                       number = 10,   
                       verboseIter = TRUE)
  
  param_grid <- expand.grid(nrounds = c(50, 100),  
                            max_depth = c(3, 5),  
                            eta = c(0.01, 0.05), 
                            gamma = c(0, 0.1), 
                            colsample_bytree = c(0.6, 0.8),
                            min_child_weight = c(3, 5),  
                            subsample = c(0.7))  
  
  xgb_model <- train(x = X, 
                     y = y, 
                     method = "xgbTree", 
                     trControl = ctrl, 
                     tuneGrid = param_grid)
  
  
  best_params <- xgb_model$bestTune
  print(best_params)
  
  params <- list(objective = "binary:logistic", 
                 eval_metric = "auc", 
                 eta = best_params$eta, 
                 max_depth = best_params$max_depth, 
                 gamma = best_params$gamma, 
                 colsample_bytree = best_params$colsample_bytree, 
                 min_child_weight = best_params$min_child_weight, 
                 subsample = best_params$subsample)
  
  dtrain <- xgb.DMatrix(data = as.matrix(X), label = y)
  
  nrounds <- best_params$nrounds
  
  cv_results <- xgb.cv(params = params,
                       data = dtrain,
                       nrounds = 500,   
                       nfold = 10,     
                       showsd = TRUE,              
                       stratified = TRUE,  
                       early_stopping_rounds = 10,  
                       verbose = TRUE)            
  
  best_nrounds <- cv_results$best_iteration
  
  xgb_model_final <- xgb.train(params = params, 
                               data = dtrain, 
                               nrounds = best_nrounds)
  
  shap_xgboost <- shapviz(xgb_model_final, X_pred = as.matrix(X))
  
  
  new_names <- colnames(X)
  names <- assay[which(assay$pro %in% new_names), c('pro', 'Assay')]
  names <- names[match(new_names, names$pro), ]
  rownames(names) <- 1:nrow(names)
  
  p <- sv_importance(shap_xgboost, kind = "beeswarm") + theme_bw() +
    scale_y_discrete(labels = setNames(names$Assay, new_names))
  p2 <- sv_importance(shap_xgboost) + theme_bw() +
    scale_y_discrete(labels = setNames(names$Assay, new_names))
  
  ggsave(filename = paste0('results/02_shap_', outcome, '.png'), 
         plot = ggarrange(
           p, p2, 
           nrow = 1, ncol = 2, common.legend = FALSE, widths = c(0.6, 0.4)
         ), 
         width = 10, height = 6, dpi = 300)
  
  importance <- sv_importance(shap_xgboost)
  importance <- as.data.frame(importance$data)
  importance$outcome <- outcome
  importance$exposure <- 'Protein'
  importance_df <- rbind(importance_df, importance)
  
  data[, paste0(outcome, '_ProRS')] <- predict(xgb_model, data[, lasso_protein])
  
  auc <- data.frame(outcome = outcome,
                    test = rocit(data[, paste0(outcome, '_ProRS')], data[, outcome])$AUC,
                    exposure = 'all')
  auc_df <- rbind(auc_df, auc)
  
  top10_protein <- as.character(importance[1:10, 'feature'])
  
  X <- df[, top10_protein]
  y <- df[, outcome]
  
  non_na_indices <- complete.cases(X, y)
  X <- X[non_na_indices, ]
  y <- y[non_na_indices]
  dtrain_simple <- xgb.DMatrix(data = as.matrix(X), label = y)
  
  param_grid <- expand.grid(nrounds = c(50, 100), 
                            max_depth = c(3, 5), 
                            eta = c(0.01, 0.05),
                            gamma = c(0, 0.1),   
                            colsample_bytree = c(0.6, 0.8), 
                            min_child_weight = c(3, 5),    
                            subsample = c(0.7))          
  
  xgb_model_simple_cv <- train(x = X, 
                            y = y, 
                            method = "xgbTree", 
                            trControl = ctrl, 
                            tuneGrid = param_grid)
  
  
  best_params_simple <- xgb_model_simple_cv$bestTune
  print(best_params_simple)

  
  params_simple <- list(objective = "binary:logistic", 
                        eval_metric = "auc", 
                        eta = best_params_simple$eta, 
                        max_depth = best_params_simple$max_depth, 
                        gamma = best_params_simple$gamma, 
                        colsample_bytree = best_params_simple$colsample_bytree, 
                        min_child_weight = best_params_simple$min_child_weight, 
                        subsample = best_params_simple$subsample)
  
  cv_results_simple <- xgb.cv(params = params_simple,
                              data = dtrain_simple,
                              nrounds = 500,  
                              nfold = 10, 
                              early_stopping_rounds = 10, 
                              metrics = "auc",
                              maximize = TRUE,
                              stratified = TRUE,
                              verbose = TRUE)
  
  best_nrounds_simple <- cv_results_simple$best_iteration
  
  xgb_model_simple <- xgb.train(params = params_simple,
                                data = dtrain_simple,
                                nrounds = best_nrounds_simple,
                                verbose = TRUE)
  
  data[, paste0(outcome, '_simpleRS')] <- predict(xgb_model_simple_cv, data[, top10_protein])
  
  auc <- data.frame(outcome = outcome,
                    test = rocit(data[, paste0(outcome, '_simpleRS')], data[, outcome])$AUC,
                    exposure = 'simple')
  auc_df <- rbind(auc_df, auc)
  
}

ProRS = data[, c('f.eid',paste0(outcomes, '_ProRS'),paste0(outcomes, '_simpleRS'))]


# Save the variable importance results to a CSV file
write.csv(importance_df, file = 'Results_wide_association/shap_cv_df_normalization_chronic.csv', row.names = FALSE)
write.csv(auc_df, file = 'Results_wide_association/auc_cv_xgboost.csv', row.names = FALSE)

save(ProRS, file = 'ProRS_xgboost_cv.RData')
