#===== Model =====
library(caret)
library(foreach)
library(doParallel)
library(dplyr)
library(rms)
library(ROCit)
library(cutpointr)
library(boot)
library(isotone)

# load data
load('Data_analysis_R2.RData')
data = data %>% filter(!is.na(missing_per_protein))

load('ProRS_xgboost_cv.RData')

df = merge(data,ProRS,by='f.eid',all.y=T)

life_style = c('Sleeplessness','Fed_up','tiredness','nerves','life_stress','BMI_cat')
for (variable in life_style) {
  df[[variable]] = as.numeric(df[[variable]])
}

outcomes = c('chronic_pain_status',
             'Head_and_Face_Pain',
             'Abdominal_Pain','Neck_Shoulder_and_Back_Hip_and_Knee_Pain',
             'Generalized_Pain')


results <- list()

for (out in outcomes) {
  # simple_pro = pro_im[which(pro_im$outcome == out),]$pro
  
  df2 = df
  df2[, out] <- as.factor(df2[, out])
  df2 = df2[!is.na(df2[, out]),]
  df2 = df2[,c(out, life_style,paste0(out,'_ProRS'),paste0(out,'_simpleRS'))]
  df2 = na.omit(df2)
  set.seed(998)
  inTraining <- createDataPartition(df2[, out], p = 0.8, list = FALSE)
  train <- df2[inTraining, ]
  test <- df2[-inTraining, ]
  
  cat(out, '====', '\n')
  results[[paste0(out, '_test')]] <- test
  # results[[paste0(out, '_train')]] <- train
  
  trControl <- trainControl(method = "cv", 
                            number = 10,
                            summaryFunction = defaultSummary,  # Default summary function for classification
                            # classProbs = TRUE, 
                            verboseIter = F)
  
  
  # Compute AUC function
  compute_auc <- function(data, indices) {
    data_boot <- data[indices, ]
    roc_curve <- rocit(score = data_boot$predicted, class = data_boot$actual)
    auc_value <- roc_curve$AUC
    return(auc_value)
  }
  
  # Bootstrap confidence interval
  n_bootstrap <- 1000
  
  # ===== bench model + lifestyle =====
  cat('bench model + lifestyle', '\n')
  life_model <- caret::train(
    x = train[, c( life_style)],
    y = train[, out],
    method = 'gbm',
    trControl = trControl
  )
  results[[paste0(out, '_m.life')]] <- life_model
  
  # Prediction in train data
  life_prob1 <- predict(life_model, newdata = train[, c( life_style)], type = "prob")[, '1']
  results[[paste0(out, '_life.prob1')]] <- life_prob1
  
  # ROC curve
  life_auc1 <- rocit(score = life_prob1, class = as.factor(train[, out]))$AUC
  results[[paste0(out, '_life.auc1')]] <- life_auc1
  
  # Prediction in test data
  life_prob2 <- predict(life_model, newdata = test[, c( life_style)], type = "prob")[, '1']
  # Calibration using Platt Scaling
  platt_model <- glm(test[, out] ~ life_prob2, family = binomial)
  life_prob2_cali <- predict(platt_model, newdata = data.frame(predicted_probs = life_prob2), type = "response")
  
  results[[paste0(out, '_life.prob2')]] <- life_prob2_cali
  
  # ROC curve
  life_auc2 <- rocit(score = life_prob2, class = as.factor(test[, out]))$AUC
  results[[paste0(out, '_life.auc2')]] <- life_auc2
  
  # Bootstrap confidence interval
  life_ci_auc2 <- boot.ci(boot(data = data.frame(actual = as.numeric(test[, out]) - 1, predicted = life_prob2),
                               statistic = compute_auc, R = n_bootstrap),
                          type = "perc", conf = 0.95)$percent[4:5]
  
  results[[paste0(out, '_life.ci_auc2')]] <- life_ci_auc2
  cat('\n', life_auc2, '\n')
  
  
  # predict all data
  life_prob_all <- predict(life_model, newdata = df[, c(life_style)], type = "prob")[, '1']
  df[,paste0(out,'_cs')] = life_prob_all
  # ===== ProRS =====
  cat('ProRS', '\n')
  
  # Prediction in train data
  # ProRS_prob1 <- predict(ProRS_model, newdata = train[, c(paste0(out,'_ProRS'))], type = "prob")[, '1']
  ProRS_prob1 =  train[, c(paste0(out,'_ProRS'))]
  platt_model <- glm(train[, out] ~ ProRS_prob1, family = binomial)
  ProRS_prob1 <- predict(platt_model, newdata = data.frame(predicted_probs = ProRS_prob1), type = "response")
  results[[paste0(out, '_ProRS.prob1')]] <- ProRS_prob1
  # ROC curve
  ProRS_auc1 <- rocit(score = ProRS_prob1, class = as.factor(train[, out]))$AUC
  results[[paste0(out, '_ProRS.auc1')]] <- ProRS_auc1
  
  # Prediction in test data
  # ProRS_prob2 <- predict(ProRS_model, newdata = test[, c(paste0(out,'_ProRS'))], type = "prob")[, '1']
  ProRS_prob2 = test[, c(paste0(out,'_ProRS'))]
  # Calibration using Platt Scaling
  
  
  platt_model <- glm(test[, out] ~ ProRS_prob2, family = binomial)
  ProRS_prob2 <- predict(platt_model, newdata = data.frame(predicted_probs = ProRS_prob2), type = "response")
  results[[paste0(out, '_ProRS.prob2')]] <- ProRS_prob2
  
  # ROC curve
  ProRS_auc2 <- rocit(score = ProRS_prob2, class = as.factor(test[, out]))$AUC
  results[[paste0(out, '_ProRS.auc2')]] <- ProRS_auc2
  
  # Bootstrap confidence interval
  ProRS_ci_auc2 <- boot.ci(boot(data = data.frame(actual = as.numeric(test[, out]) - 1, predicted = ProRS_prob2),
                                statistic = compute_auc, R = n_bootstrap),
                           type = "perc", conf = 0.95)$percent[4:5]
  
  results[[paste0(out, '_ProRS.ci_auc2')]] <- ProRS_ci_auc2
  cat('\n', ProRS_auc2, '\n')
  
  
  # ===== bench model + lifestyle + ProRS =====
  cat('bench model + lifestyle + ProRS', '\n')
  
  
  pro_model <- caret::train(
    x = train[, c( life_style, paste0(out,'_ProRS'))],
    y = train[, out],
    method = 'gbm',
    trControl = trControl
  )
  
  results[[paste0(out, '_m.pro')]] <- pro_model
  
  
  # Prediction in train data
  pro_prob1 <- predict(pro_model, newdata = train[, c(life_style, paste0(out,'_ProRS'))], type = "prob")[,'1']
  results[[paste0(out, '_pro.prob1')]] <- pro_prob1
  
  # ROC curve
  pro_auc1 <- rocit(score = pro_prob1, class = as.factor(train[, out]))$AUC
  results[[paste0(out, '_pro.auc1')]] <- pro_auc1
  
  # Prediction in test data
  pro_prob2 <- predict(pro_model, newdata = test[, c( life_style, paste0(out,'_ProRS'))], type = "prob")[, '1']
  # Calibration using Platt Scaling
  platt_model <- glm(test[, out] ~ pro_prob2, family = binomial)
  pro_prob2_cali <- predict(platt_model, newdata = data.frame(predicted_probs = pro_prob2), type = "response")
  
  results[[paste0(out, '_pro.prob2')]] <- pro_prob2
  
  # ROC curve
  pro_auc2 <- rocit(score = pro_prob2, class = as.factor(test[, out]))$AUC
  results[[paste0(out, '_pro.auc2')]] <- pro_auc2
  
  # Bootstrap confidence interval
  pro_ci_auc2 <- boot.ci(boot(data = data.frame(actual = as.numeric(test[, out]) - 1, predicted = pro_prob2),
                              statistic = compute_auc, R = n_bootstrap),
                         type = "perc", conf = 0.95)$percent[4:5]
  
  results[[paste0(out, '_pro.ci_auc2')]] <- pro_ci_auc2
  cat('\n', pro_auc2, '\n')
  
  
  
  # ===== simple =====
  cat('simple', '\n')
  
  # Prediction in train data
  simple_prob1 = train[, paste0(out, '_simpleRS')]
  results[[paste0(out, '_simple.prob1')]] <- simple_prob1
  
  # ROC curve
  simple_auc1 <- rocit(score = simple_prob1, class = as.factor(train[, out]))$AUC
  results[[paste0(out, '_simple.auc1')]] <- simple_auc1
  
  # Prediction in test data
  # simple_prob2 <- predict(simple_model, newdata = test[, c( simple_pro)], type = "prob")[, '1']
  simple_prob2 <- test[, paste0(out, '_simpleRS')]
  # Calibration using Platt Scaling
  platt_model <- glm(test[, out] ~ simple_prob2, family = binomial)
  simple_prob2_cali <- predict(platt_model, newdata = data.frame(predicted_probs = simple_prob2), type = "response")
  
  results[[paste0(out, '_simple.prob2')]] <- simple_prob2_cali
  
  # ROC curve
  simple_auc2 <- rocit(score = simple_prob2, class = as.factor(test[, out]))$AUC
  results[[paste0(out, '_simple.auc2')]] <- simple_auc2
  
  # Bootstrap confidence interval
  simple_ci_auc2 <- boot.ci(boot(data = data.frame(actual = as.numeric(test[, out]) - 1, predicted = simple_prob2),
                                 statistic = compute_auc, R = n_bootstrap),
                            type = "perc", conf = 0.95)$percent[4:5]
  
  results[[paste0(out, '_simple.ci_auc2')]] <- simple_ci_auc2
  cat('\n', simple_auc2, '\n')
  
  
  # ===== simple_life =====
  cat('simple_life', '\n')
  simple_life_model <- caret::train(
    x = train[, c(life_style, paste0(out, '_simpleRS'))],
    y = train[, out],
    method = 'gbm',
    trControl = trControl
  )
  results[[paste0(out, '_m.simple_life')]] <- simple_life_model
  
  # Prediction in train data
  simple_life_prob1 <- predict(simple_life_model, newdata = train[, c(life_style, paste0(out, '_simpleRS'))], type = "prob")[, '1']
  results[[paste0(out, '_simple_life.prob1')]] <- simple_life_prob1
  
  # ROC curve
  simple_life_auc1 <- rocit(score = simple_life_prob1, class = as.factor(train[, out]))$AUC
  results[[paste0(out, '_simple_life.auc1')]] <- simple_life_auc1
  
  # Prediction in test data
  simple_life_prob2 <- predict(simple_life_model, newdata = test[, c(life_style, paste0(out, '_simpleRS'))], type = "prob")[, '1']
  # Calibration using Platt Scaling
  platt_model <- glm(test[, out] ~ simple_life_prob2, family = binomial)
  simple_life_prob2_cali <- predict(platt_model, newdata = data.frame(predicted_probs = simple_life_prob2), type = "response")
  
  results[[paste0(out, '_simple_life.prob2')]] <- simple_life_prob2_cali
  
  # ROC curve
  simple_life_auc2 <- rocit(score = simple_life_prob2, class = as.factor(test[, out]))$AUC
  results[[paste0(out, '_simple_life.auc2')]] <- simple_life_auc2
  
  # Bootstrap confidence interval
  simple_life_ci_auc2 <- boot.ci(boot(data = data.frame(actual = as.numeric(test[, out]) - 1, predicted = simple_life_prob2),
                                      statistic = compute_auc, R = n_bootstrap),
                                 type = "perc", conf = 0.95)$percent[4:5]
  
  results[[paste0(out, '_simple_life.ci_auc2')]] <- simple_life_ci_auc2
  cat('\n', simple_life_auc2, '\n')
  
  
}





save(results, file = 'Model.results_PROS_chronic.RData')
