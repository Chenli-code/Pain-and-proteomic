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
protein_selected = readRDS('D:/Desktop/PhD/10-学术项目/34-chronic_pain/Prediction_model/DataPreparation/name_protein.rds')
load('D:/Desktop/PhD/10-学术项目/34-chronic_pain/Temp_folder2/Data_analysis_R2.RData')
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

assay = readRDS('D:/Desktop/PhD/10-学术项目/34-chronic_pain/Prediction_model/DataPreparation/PRO_Data_Preprocessing/protein_assay.rds')
assay = assay %>% mutate(pro = paste0('P',coding))



auc_df = data.frame()
importance_df = data.frame()
model_list_p = list()

df = data

##### Protein #####
for (outcome in outcomes) {
  
  cat(outcome, ' - ')
  # Select the proteins associated with the current outcome
  # Use Lasso (L1 Regularization) to further refine important features
  # set.seed(123)
  # df <- data[, c(outcome, protein_selected)]
  # df <- df[!is.na(df[, outcome]), ]
  # y <- df[, outcome]
  # X_selected_df <- df[, protein_selected]
  # x_matrix <- model.matrix(~., data = X_selected_df)[, -1]
  # lasso_model <- cv.glmnet(x_matrix, y, alpha = 1, family = "binomial")
  # coef_matrix <- coef(lasso_model, s = "lambda.min")
  # coef_df <- as.data.frame(as.matrix(coef_matrix))
  # lasso_protein <- rownames(coef_df)[which(coef_df != 0)[-1]]
  
  # Step 1: 选择相关蛋白质 ----
  # 从Protein_association表中选择与当前outcome相关的蛋白
  lasso_protein <- Protein_association[which(Protein_association$outcome == outcome), 'exposure']
  
  # Step 2: 构建XGBoost模型 ----
  # 定义特征和目标变量
  
  X <- df[, lasso_protein]
  y <- df[, outcome]
  non_na_indices <- complete.cases(X, y)
  X <- X[non_na_indices, ]
  y <- y[non_na_indices]
  # # Step 3: 使用caret包进行Bootstrapping调参 ----
  # # 将数据集转换为trainControl对象，使用Bootstrapping
  # cat('使用Bootstrapping进行XGBoost模型调参...\n')
  # ctrl <- trainControl(method = "boot",   # 使用 Bootstrapping
  #                      number = 100,      # Bootstrap 重采样次数
  #                      verboseIter = TRUE)
  # # 设置参数网格
  # param_grid <- expand.grid(nrounds = c(50, 100),  # 迭代轮数（nrounds），降低运行时间
  #                           max_depth = c(3, 5),   # 最大树深度（max_depth）
  #                           eta = c(0.01, 0.05),   # 学习率（eta）
  #                           gamma = c(0, 0.1),     # 树分裂所需的最小损失减少值
  #                           colsample_bytree = c(0.6, 0.8),  # 特征子采样比例（colsample_bytree）
  #                           min_child_weight = c(3, 5),      # 叶子节点的最小权重和（min_child_weight）
  #                           subsample = c(0.7))              # 样本子采样比例（subsample）
  # # 使用train()函数进行参数调优
  # xgb_model <- train(x = X, 
  #                    y = y, 
  #                    method = "xgbTree", 
  #                    trControl = ctrl, 
  #                    tuneGrid = param_grid)
  # cat('调参完成。\n')
  
  # Step 3: 使用caret包进行10折交叉验证调参 ----
  # 将数据集转换为trainControl对象，使用10折交叉验证
  cat('使用10折交叉验证进行XGBoost模型调参...\n')
  ctrl <- trainControl(method = "cv",   # 使用 10 折交叉验证
                       number = 10,      # 10 折交叉验证次数
                       verboseIter = TRUE)
  
  # 设置参数网格
  param_grid <- expand.grid(nrounds = c(50, 100),  # 迭代轮数（nrounds），降低运行时间
                            max_depth = c(3, 5),   # 最大树深度（max_depth）
                            eta = c(0.01, 0.05),   # 学习率（eta）
                            gamma = c(0, 0.1),     # 树分裂所需的最小损失减少值
                            colsample_bytree = c(0.6, 0.8),  # 特征子采样比例（colsample_bytree）
                            min_child_weight = c(3, 5),      # 叶子节点的最小权重和（min_child_weight）
                            subsample = c(0.7))              # 样本子采样比例（subsample）
  
  # 使用train()函数进行参数调优
  xgb_model <- train(x = X, 
                     y = y, 
                     method = "xgbTree", 
                     trControl = ctrl, 
                     tuneGrid = param_grid)
  cat('10折交叉验证调参完成。\n')
  
  
  
  
  # 输出最佳参数配置
  best_params <- xgb_model$bestTune
  print(best_params)
  
  # Step 4: 使用最佳参数训练最终模型 ----
  cat('使用最佳参数训练最终模型并增加早停法...\n')
  params <- list(objective = "binary:logistic", 
                 eval_metric = "auc", 
                 eta = best_params$eta, 
                 max_depth = best_params$max_depth, 
                 gamma = best_params$gamma, 
                 colsample_bytree = best_params$colsample_bytree, 
                 min_child_weight = best_params$min_child_weight, 
                 subsample = best_params$subsample)
  
  # 重新将整个数据集用于训练
  dtrain <- xgb.DMatrix(data = as.matrix(X), label = y)
  
  nrounds <- best_params$nrounds
  
  # 使用最佳参数训练最终模型，增加 early stopping
  # 使用 xgb.cv 进行10折交叉验证
  cv_results <- xgb.cv(params = params,
                       data = dtrain,
                       nrounds = 500,                   # 最大迭代次数
                       nfold = 10,                      # 10 折交叉验证
                       showsd = TRUE,                   # 显示标准差
                       stratified = TRUE,               # 是否分层抽样
                       early_stopping_rounds = 10,      # 早停法，连续10轮无改善则停止
                       verbose = TRUE)                  # 显示详细输出
  
  # 输出最佳迭代次数
  best_nrounds <- cv_results$best_iteration
  cat('交叉验证完成，最佳迭代次数为:', best_nrounds, '\n')
  
  # 使用最佳参数和迭代次数训练最终模型
  xgb_model_final <- xgb.train(params = params, 
                               data = dtrain, 
                               nrounds = best_nrounds)
  
  cat('最终模型训练完成。\n')
  
  
  # Step 4: 计算SHAP值并绘制重要性蜂群图 ----
  cat('计算SHAP值...\n')
  shap_xgboost <- shapviz(xgb_model_final, X_pred = as.matrix(X))
  cat('SHAP值计算完成。\n')
  
  # Step 5: 绘制变量重要性蜂群图 ----
  new_names <- colnames(X)
  names <- assay[which(assay$pro %in% new_names), c('pro', 'Assay')]
  names <- names[match(new_names, names$pro), ]
  rownames(names) <- 1:nrow(names)
  
  # 绘制变量重要性蜂群图
  cat('绘制并保存SHAP重要性图...\n')
  p <- sv_importance(shap_xgboost, kind = "beeswarm") + theme_bw() +
    scale_y_discrete(labels = setNames(names$Assay, new_names))
  p2 <- sv_importance(shap_xgboost) + theme_bw() +
    scale_y_discrete(labels = setNames(names$Assay, new_names))
  
  # print(ggarrange(
  #   p, p2, 
  #   labels = c('Distribution of SHAP', 'Mean SHAP'),
  #   nrow = 1, ncol = 2, common.legend = FALSE, widths = c(0.5, 0.5)
  # ))
  # export::graph2ppt(file = "results/02_shap.pptx", width = 10, height = 13.5/2, append = TRUE)
  ggsave(filename = paste0('results/02_shap_', outcome, '.png'), 
         plot = ggarrange(
           p, p2, 
           # labels = c('Distribution of SHAP', 'Mean SHAP'),
           nrow = 1, ncol = 2, common.legend = FALSE, widths = c(0.6, 0.4)
         ), 
         width = 10, height = 6, dpi = 300)
  cat('SHAP图保存完成。\n')
  
  # Step 6: 保存SHAP重要性数据 ----
  importance <- sv_importance(shap_xgboost)
  importance <- as.data.frame(importance$data)
  importance$outcome <- outcome
  importance$exposure <- 'Protein'
  importance_df <- rbind(importance_df, importance)
  
  # Step 7: 计算风险评分 ----
  cat('计算风险评分并预测...\n')
  data[, paste0(outcome, '_ProRS')] <- predict(xgb_model, data[, lasso_protein])
  cat('风险评分计算完成。\n')
  
  # Step 8: 计算AUC值 ----
  cat('计算AUC值...\n')
  auc <- data.frame(outcome = outcome,
                    test = rocit(data[, paste0(outcome, '_ProRS')], data[, outcome])$AUC,
                    exposure = 'all')
  auc_df <- rbind(auc_df, auc)
  
  # Step 9: 使用Top 10蛋白质构建简单模型 ----
  cat('使用Top 10蛋白质构建简单模型...')
  top10_protein <- as.character(importance[1:10, 'feature'])
  
  # 定义训练集特征和目标变量
  X <- df[, top10_protein]
  y <- df[, outcome]
  
  non_na_indices <- complete.cases(X, y)
  X <- X[non_na_indices, ]
  y <- y[non_na_indices]
  dtrain_simple <- xgb.DMatrix(data = as.matrix(X), label = y)
  
  # 设置参数网格
  param_grid <- expand.grid(nrounds = c(50, 100),  # 迭代轮数（nrounds）
                            max_depth = c(3, 5),   # 最大树深度（max_depth）
                            eta = c(0.01, 0.05),   # 学习率（eta）
                            gamma = c(0, 0.1),     # 树分裂所需的最小损失减少值
                            colsample_bytree = c(0.6, 0.8),  # 特征子采样比例（colsample_bytree）
                            min_child_weight = c(3, 5),      # 叶子节点的最小权重和（min_child_weight）
                            subsample = c(0.7))              # 样本子采样比例（subsample）
  
  # 使用train()函数进行参数调优
  xgb_model_simple_cv <- train(x = X, 
                            y = y, 
                            method = "xgbTree", 
                            trControl = ctrl, 
                            tuneGrid = param_grid)
  cat('简单模型调参完成。')
  
  # 获取最佳参数配置
  best_params_simple <- xgb_model_simple_cv$bestTune
  print(best_params_simple)
  
  # 使用最佳参数和早停法训练最终模型
  cat('使用最佳参数和整个数据集训练简单模型...')
  
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
                              nrounds = 500,  # 设置较大轮数，以允许早停法确定最佳迭代次数
                              nfold = 10,      # 10折交叉验证
                              early_stopping_rounds = 10,  # 早停法
                              metrics = "auc",
                              maximize = TRUE,
                              stratified = TRUE,
                              verbose = TRUE)
  
  # 获取最佳迭代次数
  best_nrounds_simple <- cv_results_simple$best_iteration
  cat('交叉验证最佳迭代次数为:', best_nrounds_simple, '\n')
  
  # 使用最佳参数和最佳迭代次数训练简单模型
  cat('使用最佳参数和整个数据集训练简单模型...')
  xgb_model_simple <- xgb.train(params = params_simple,
                                data = dtrain_simple,
                                nrounds = best_nrounds_simple,
                                verbose = TRUE)
  cat('简单模型训练完成。\n')
  
  # Step 10: 预测并计算简单模型的AUC值 ----
  cat('计算简单模型的风险评分和AUC值...')
  data[, paste0(outcome, '_simpleRS')] <- predict(xgb_model_simple_cv, data[, top10_protein])
  # 计算AUC值
  auc <- data.frame(outcome = outcome,
                    test = rocit(data[, paste0(outcome, '_simpleRS')], data[, outcome])$AUC,
                    exposure = 'simple')
  auc_df <- rbind(auc_df, auc)
  cat('简单模型的AUC值计算完成。\n')
}

ProRS = data[, c('f.eid',paste0(outcomes, '_ProRS'),paste0(outcomes, '_simpleRS'))]


# Save the variable importance results to a CSV file
write.csv(importance_df, file = 'Results_wide_association/shap_cv_df_normalization_chronic.csv', row.names = FALSE)
write.csv(auc_df, file = 'Results_wide_association/auc_cv_xgboost.csv', row.names = FALSE)

save(ProRS, file = 'D:/Desktop/PhD/10-学术项目/34-chronic_pain/Temp_folder2/ProRS_xgboost_cv.RData')
