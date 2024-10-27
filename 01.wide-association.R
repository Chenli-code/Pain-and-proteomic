library(dplyr)
library(ROCit)
library(foreach)
library(parallel)
library(doParallel)
# rm(list = ls())

# load data
setwd("Prediction_model_R2")
protein_selected = readRDS('name_protein.rds')

##### chronic pain #####
load('Data_analysis_R2.RData')

data_protein = data %>% filter(!is.na(missing_per_protein))

outcomes = c('chronic_pain_status',
             'Head_and_Face_Pain','Neck_Shoulder_and_Back_Pain',
             'Abdominal_Pain','Hip_and_Knee_Pain',
             'Generalized_Pain','Neck_Shoulder_and_Back_Hip_and_Knee_Pain')

# Setup parallel backend
# Create a cluster using all available cores minus one
cl <- makePSOCKcluster(5)
# Register the cluster as the default
registerDoParallel(cl)
##### Protein-wide association #####
# Protein_association = run_logistic_regression(data,protein_selected,outcomes)
##### Protein-wide association #####
Protein_association <- data.frame(
  exposure = NA,
  outcome = NA,
  estimate = NA,
  se = NA,
  p_value = NA,
  stringsAsFactors = FALSE
)
for (outcome in outcomes) {
  cat(outcome, '\n')
  results <- foreach(pro = protein_selected, .combine = 'rbind', .packages = c("stats")) %dopar% {
    cat(pro, '-', outcome, '\n')
    # Building the model
    formula_string <- paste(outcome, " ~ Age + Gender + ", pro)
    model <- tryCatch({
      glm(as.formula(formula_string), data = data_protein, family = binomial())
    }, error = function(e) {
      print(paste("Error fitting model for", pro, "and", outcome))
      return(NULL)  # Return NULL if the model fails to fit
    })
    
    # Check if the model was successfully fitted
    if (!is.null(model) && pro %in% names(coef(model))) {
      est <- coef(summary(model))[pro, "Estimate"]
      stande  <- coef(summary(model))[pro, "Std. Error"]
      p_val <- coef(summary(model))[pro, "Pr(>|z|)"]
    } else {
      est <- NA
      stande <- NA
      p_val <- NA
    }
    
    # Construct result data frame
    result <- data.frame(
      exposure = pro,
      outcome = outcome,
      estimate = est,
      se = stande,
      p_value = p_val,
      stringsAsFactors = FALSE
    )
    
    return(result)
  }
  # # Combine the output from foreach
  Protein_association = rbind(Protein_association,results)
}
# stopCluster(cl)
write.csv(Protein_association,file = 'Results_wide_association/Protein_association-chronic.csv',row.names = F)



##### acute pain #####
load('Data_analysis_acute_R2.RData')

data_protein = data %>% filter(!is.na(missing_per_protein))

outcomes = c('acute_pain_status',
             'Head_and_Face_Pain','Neck_Shoulder_and_Back_Pain',
             'Abdominal_Pain','Hip_and_Knee_Pain',
             'Generalized_Pain','Neck_Shoulder_and_Back_Hip_and_Knee_Pain')

# Setup parallel backend
# Create a cluster using all available cores minus one
cl <- makePSOCKcluster(5)
# Register the cluster as the default
registerDoParallel(cl)
##### Protein-wide association #####
# Protein_association = run_logistic_regression(data,protein_selected,outcomes)
##### Protein-wide association #####
Protein_association <- data.frame(
  exposure = NA,
  outcome = NA,
  estimate = NA,
  se = NA,
  p_value = NA,
  stringsAsFactors = FALSE
)
for (outcome in outcomes) {
  cat(outcome, '\n')
  results <- foreach(pro = protein_selected, .combine = 'rbind', .packages = c("stats")) %dopar% {
    cat(pro, '-', outcome, '\n')
    # Building the model
    formula_string <- paste(outcome, " ~ Age + Gender + ", pro)
    model <- tryCatch({
      glm(as.formula(formula_string), data = data_protein, family = binomial())
    }, error = function(e) {
      print(paste("Error fitting model for", pro, "and", outcome))
      return(NULL)  # Return NULL if the model fails to fit
    })
    
    # Check if the model was successfully fitted
    if (!is.null(model) && pro %in% names(coef(model))) {
      est <- coef(summary(model))[pro, "Estimate"]
      stande  <- coef(summary(model))[pro, "Std. Error"]
      p_val <- coef(summary(model))[pro, "Pr(>|z|)"]
    } else {
      est <- NA
      stande <- NA
      p_val <- NA
    }
    
    # Construct result data frame
    result <- data.frame(
      exposure = pro,
      outcome = outcome,
      estimate = est,
      se = stande,
      p_value = p_val,
      stringsAsFactors = FALSE
    )
    
    return(result)
  }
  # # Combine the output from foreach
  Protein_association = rbind(Protein_association,results)
}
# stopCluster(cl)
write.csv(Protein_association,file = 'Results_wide_association/Protein_association-acute.csv',row.names = F)


stopImplicitCluster()


