library(dplyr)
library(ROCit)
library(foreach)
library(parallel)
library(doParallel)
# rm(list = ls())

# load data
setwd("D:/Desktop/PhD/10-学术项目/34-chronic_pain/Prediction_model_R2")
protein_selected = readRDS('D:/Desktop/PhD/10-学术项目/34-chronic_pain/Prediction_model/DataPreparation/name_protein.rds')

##### chronic pain #####
load('D:/Desktop/PhD/10-学术项目/34-chronic_pain/Temp_folder2/Data_analysis_R2.RData')

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
load('D:/Desktop/PhD/10-学术项目/34-chronic_pain/Temp_folder2/Data_analysis_acute_R2.RData')

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




# Plot of the Volcano -----------------------------------------------------
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggrepel)
coding_assay = readRDS('D:/Desktop/PhD/10-学术项目/34-chronic_pain/Prediction_model/DataPreparation/PRO_Data_Preprocessing/protein_assay.rds')
coding_assay = coding_assay %>% mutate(pro = paste0('P',coding))
##### Chronic pain #####

# volcano plot of chronic pain --------------------------------------------
dfp = read.csv('Results_wide_association/Protein_association-chronic.csv')
dfp = merge(dfp,coding_assay,by.y='pro',by.x='exposure',all.x=T)

outcome_s = c(#'overall_pain','acute_pain_status',
  'chronic_pain_status',
  'Head_and_Face_Pain',
  'Abdominal_Pain',
  'Neck_Shoulder_and_Back_Hip_and_Knee_Pain',
  'Generalized_Pain')

lable =  c(#'Overall pain','Acute pain',
  'Overall chronic pain',
  'Head and face pain',
  'Abdominal pain',
  'Musculoskeletal pain',
  'Widespread pain')

for (i in 1:5) {
  
  outs = outcome_s[i]
  out_label = lable[i]
  dfp1 = dfp[which(dfp$outcome == outs),]
  dfp1$log_p_value <- -log10(dfp1$p_value)
  
  dfp1 <- dfp1 %>%
    mutate(Panel2 = case_when(
      estimate > 0 & p_value < 0.05/2920 ~ "Positive association",
      estimate < 0 & p_value < 0.05/2920 ~ "Negative association",
      TRUE ~ "Not Significant"
    ))
  # 找出负相关和正相关的最大 p 值的前五个
  if(outs == 'Generalized_Pain'){
    pro10 = dfp1 %>% filter(outcome == outs) %>% arrange(desc(log_p_value)) %>% 
      mutate(dir = ifelse(estimate >0,'positive','negative'))%>%group_by(dir) %>% 
      slice_head(n = 10)
    top = dfp1 %>% filter(exposure %in% pro10$exposure)
    outcomes = dfp1 %>% filter(outcome == outs) %>% 
      mutate(dir = ifelse(estimate >0,'positive','negative'))%>%
      group_by(dir) %>% 
      arrange(desc(log_p_value)) %>%
      slice_head(n = 10)
  }else{
    pro10 = dfp1 %>% filter(outcome == outs) %>% arrange(desc(log_p_value)) %>% 
      mutate(dir = ifelse(estimate >0,'positive','negative'))%>%group_by(dir) %>% 
      slice_head(n = 3)
    top = dfp1 %>% filter(exposure %in% pro10$exposure)
    outcomes = dfp1 %>% filter(outcome == outs) %>% 
      mutate(dir = ifelse(estimate >0,'positive','negative'))%>%
      group_by(dir) %>% 
      arrange(desc(log_p_value)) %>%
      slice_head(n = 3)
  }
  
  # 定义8种颜色，分别对应8个不同的Panel名称
  colors <- c("Cardiometabolic" = "#374E5599",  # 色系中的第一种颜色
              "Cardiometabolic II" = "#374E554C",  # 第一种颜色的 70% 透明度
              "Inflammation" = "#DF8F4499",
              "Inflammation II" = "#DF8F444C",  # 第二种颜色的 70% 透明度
              "Neurology" = "#00A1D599",
              "Neurology II" = "#00A1D54C",  # 第三种颜色的 70% 透明度
              "Oncology" = "#B2474599",
              "Oncology II" = "#B247454C",
              "Not Significant" = "grey")
  
  colors <- c(
    
    "Cardiometabolic" = "orange","Cardiometabolic II" = "green",
    "Inflammation II" = "red", "Inflammation" = "cyan",
    
    
    "Neurology" = "brown","Neurology II" = "pink",
    "Oncology" = "purple", "Oncology II" = "blue",
    "Not Significant" = "grey")
  
  threshold <- -log10(0.05/2920)
  assign(paste0('p_',outs),
         # 创建火山图
         ggplot(dfp1, aes(x = estimate, y = log_p_value)) +
           geom_point(aes(fill = Panel2,color=Panel2),shape=21,size=1,alpha = 0.6) +
           scale_fill_manual(values = c("Positive association" = "red", "Negative association" = "blue", "Not Significant" = "gray")) +
           scale_color_manual(values = c("Positive association" = "red", "Negative association" = "blue", "Not Significant" = "gray")) +
           labs(x = "Estimate",
                y = "-log10(p-value)",
                # title = unique(dfp1$outcome_label),
                title = out_label
                ) +
           xlim(-2.5,2.5)+
           theme_bw()  +
           geom_hline(yintercept = -log10(0.05/2920),size=1, linetype = "dashed", color = "orange") +
           annotate("text", x = max(dfp1$estimate), y = -log10(0.05/2920), label = paste0("Bonferroni correction"),
                    vjust = +1, hjust = 1, color = "black") +
           geom_vline(xintercept = 0,size=1, linetype = "dashed") +
           geom_text_repel(data = top, aes(label = Assay),
                           size = 3, 
                           box.padding = unit(0.35, "lines"),
                           point.padding = unit(0.3, "lines"),
                           segment.color = 'black', 
                           nudge_x = -0.2) +
           geom_point(data = top, aes(x = estimate, y = log_p_value,fill = Panel2),shape=21,size=2,color='black',alpha = 1) +
           theme(legend.title = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position = 'top')
  )
}
# Arrange the first row of plots
volcano_chronic <- ggarrange(
  p_Generalized_Pain, 
  ggarrange(p_chronic_pain_status, p_Head_and_Face_Pain,
            p_Abdominal_Pain, p_Neck_Shoulder_and_Back_Hip_and_Knee_Pain,
            nrow = 2, ncol = 2,common.legend = TRUE, legend = 'none'),
  # Labels removed as they are commented out in the original code
  nrow = 1, ncol = 2, common.legend = TRUE, legend = 'top'
)


# volcano plot of acute pain --------------------------------------------
dfp = read.csv('Results_wide_association/Protein_association-acute.csv')
dfp = merge(dfp,coding_assay,by.y='pro',by.x='exposure',all.x=T)

outcome_s = c(#'overall_pain','acute_pain_status',
  'acute_pain_status',
  'Head_and_Face_Pain',
  'Abdominal_Pain',
  'Neck_Shoulder_and_Back_Hip_and_Knee_Pain',
  'Generalized_Pain')

lable =  c(#'Overall pain','Acute pain',
  'Overall acute pain',
  'Head and face pain',
  'Abdominal pain',
  'Musculoskeletal pain',
  'Generalized pain')

for (i in 1:5) {
  
  outs = outcome_s[i]
  out_label = lable[i]
  dfp1 = dfp[which(dfp$outcome == outs),]
  dfp1$log_p_value <- -log10(dfp1$p_value)
  
  dfp1 <- dfp1 %>%
    mutate(Panel2 = case_when(
      estimate > 0 & p_value < 0.05/2920 ~ "Positive association",
      estimate < 0 & p_value < 0.05/2920 ~ "Negative association",
      TRUE ~ "Not Significant"
    ))
  # 找出负相关和正相关的最大 p 值的前五个
  pro10 = dfp1 %>% filter(outcome == outs) %>% arrange(desc(log_p_value)) %>% 
    mutate(dir = ifelse(estimate >0,'positive','negative'))%>%group_by(dir) %>% 
    slice_head(n = 10)
  top = dfp1 %>% filter(exposure %in% pro10$exposure)
  outcomes = dfp1 %>% filter(outcome == outs) %>% 
    mutate(dir = ifelse(estimate >0,'positive','negative'))%>%
    group_by(dir) %>% 
    arrange(desc(log_p_value)) %>%
    slice_head(n = 10)
  # 定义8种颜色，分别对应8个不同的Panel名称
  colors <- c("Cardiometabolic" = "#374E5599",  # 色系中的第一种颜色
              "Cardiometabolic II" = "#374E554C",  # 第一种颜色的 70% 透明度
              "Inflammation" = "#DF8F4499",
              "Inflammation II" = "#DF8F444C",  # 第二种颜色的 70% 透明度
              "Neurology" = "#00A1D599",
              "Neurology II" = "#00A1D54C",  # 第三种颜色的 70% 透明度
              "Oncology" = "#B2474599",
              "Oncology II" = "#B247454C",
              "Not Significant" = "grey")
  
  colors <- c(
    
    "Cardiometabolic" = "orange","Cardiometabolic II" = "green",
    "Inflammation II" = "red", "Inflammation" = "cyan",
    
    
    "Neurology" = "brown","Neurology II" = "pink",
    "Oncology" = "purple", "Oncology II" = "blue",
    "Not Significant" = "grey")
  
  threshold <- -log10(0.05/2920)
  assign(paste0('p_',outs),
         # 创建火山图
         ggplot(dfp1, aes(x = estimate, y = log_p_value)) +
           geom_point(aes(fill = Panel2,color=Panel2),shape=21,size=1,alpha = 0.6) +
           scale_fill_manual(values = c("Positive association" = "red", "Negative association" = "blue", "Not Significant" = "gray")) +
           scale_color_manual(values = c("Positive association" = "red", "Negative association" = "blue", "Not Significant" = "gray")) +
           labs(x = "Estimate",
                y = "-log10(p-value)",
                # title = unique(dfp1$outcome_label),
                title = out_label
           ) +
           theme_bw()  +
           geom_hline(yintercept = -log10(0.05/2920),size=1, linetype = "dashed", color = "orange") +
           annotate("text", x = max(dfp1$estimate), y = -log10(0.05/2920), label = paste0("Bonferroni correction"),
                    vjust = +1, hjust = 1, color = "black") +
           geom_vline(xintercept = 0,size=1, linetype = "dashed") +
           geom_text_repel(data = top, aes(label = Assay),
                           size = 3, 
                           box.padding = unit(0.35, "lines"),
                           point.padding = unit(0.3, "lines"),
                           segment.color = 'black', 
                           nudge_x = -0.2) +
           geom_point(data = top, aes(x = estimate, y = log_p_value,fill = Panel2),shape=21,size=2,color='black',alpha = 1) +
           theme(legend.title = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position = 'top')
  )
}
# Arrange the first row of plots
volcano_acute <- ggarrange(
  p_Generalized_Pain, p_acute_pain_status, p_Head_and_Face_Pain,
  p_Abdominal_Pain, p_Neck_Shoulder_and_Back_Hip_and_Knee_Pain,
  # Labels removed as they are commented out in the original code
  nrow = 2, ncol = 3, common.legend = TRUE, legend = 'none'
)

# Bar plot of chronic and acute pain --------------------------------------
df_chronic = read.csv('Results_wide_association/Protein_association-chronic.csv') %>% filter(p_value < 0.05/2920) %>% 
  select(exposure,outcome) %>% mutate(pain_c = 'chronic',outcome = ifelse(outcome == 'chronic_pain_status','Overall pain',outcome))
df_acute = read.csv('Results_wide_association/Protein_association-acute.csv') %>% filter(p_value < 0.05/2920) %>% 
  select(exposure,outcome) %>% mutate(pain_a = 'acute',outcome = ifelse(outcome == 'acute_pain_status','Overall pain',outcome))
df = merge(df_chronic,df_acute,by = c('exposure','outcome'),all=T)


df_plot = df %>% group_by(pain_c,pain_a,outcome) %>% summarise(n=n()) %>% 
  mutate(cat = ifelse(!is.na(pain_c)&!is.na(pain_a),'Associated with Both Acute and Chronic Pain',
                      ifelse(!is.na(pain_c),'Chronic Pain Specific','Acute Pain Specific')))
assay = readRDS('D:/Desktop/PhD/10-学术项目/34-chronic_pain/Prediction_model/DataPreparation/PRO_Data_Preprocessing/protein_assay.rds')%>% 
  mutate(exposure = paste0('P',coding))
outcomes = c(
  'Generalized_Pain','Overall pain',
  'Head_and_Face_Pain',
  'Abdominal_Pain',
  'Neck_Shoulder_and_Back_Hip_and_Knee_Pain'
  )

lable =  c(
  'Widespread pain','Overall pain',
  'Head and face pain',
  'Abdominal pain',
  'Musculoskeletal pain'
  )

df_plot = df_plot %>%filter(outcome %in% outcomes) %>% 
  mutate(outcome_label = factor(outcome, levels = rev(outcomes), labels = rev(lable)),
         cat = factor(cat,levels = rev(c('Chronic Pain Specific','Associated with Both Acute and Chronic Pain','Acute Pain Specific'))))

protein_overlap = ggplot(df_plot, aes(x = outcome_label, y = n, fill = cat)) + 
  geom_bar(stat = "identity", position = "stack",color = 'black') + 
  labs(x = "Pain", y = "Number of proteins", fill = "",title='Protein overlap') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_nejm()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'top')+
  coord_flip()


# bar_chronic -------------------------------------------------------------

df = read.csv('Results_wide_association/Protein_association-chronic.csv') %>% filter(p_value < 0.05/2923)
# assay = readRDS('H:/Prediction_model/DataPreparation/PRO_Data_Preprocessing/protein_assay.rds')%>% mutate(exposure = paste0('P',coding))
dfp = df %>% #filter(outcome == 'overall_pain') %>% 
  merge(.,assay,by='exposure') %>% group_by(Panel,outcome) %>% 
  summarise(n=n())
outcomes = c(#'overall_pain','acute_pain_status',
  'Generalized_Pain','chronic_pain_status',
  'Head_and_Face_Pain',
  'Abdominal_Pain',
  'Neck_Shoulder_and_Back_Hip_and_Knee_Pain'
  )

lable =  c(#'Overall pain','Acute pain',
  'Widespread pain','Overall pain',
  'Head and face pain',
  'Abdominal pain',
  'Musculoskeletal pain'
  )


colors <- c("Cardiometabolic" = "#374E5599",  # 色系中的第一种颜色
            "Cardiometabolic II" = "#374E554C",  # 第一种颜色的 70% 透明度
            "Inflammation" = "#DF8F4499",
            "Inflammation II" = "#DF8F444C",  # 第二种颜色的 70% 透明度
            "Neurology" = "#00A1D599",
            "Neurology II" = "#00A1D54C",  # 第三种颜色的 70% 透明度
            "Oncology" = "#B2474599",
            "Oncology II" = "#B247454C")  # 第四种颜色的 70% 透明度
total_counts <- dfp %>%filter(outcome %in% outcomes) %>% 
  mutate(outcome_label = factor(outcome, levels = outcomes, labels = lable)) %>% 
  group_by(outcome_label) %>%
  summarise(total = sum(n))

bar_chronic = dfp %>%
  filter(outcome %in% outcomes) %>% 
  mutate(outcome = factor(outcome, levels = outcomes)) %>%
  mutate(outcome_label = factor(outcome, levels = outcomes, labels = lable)) %>%
  ggplot(aes(x = outcome_label, y = n, fill = Panel)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  ylim(0,180)+
  labs(x = "", y = "Count of proteins", title = "") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5) + 
  # geom_text(data = total_counts, aes(x = outcome_label, y = total, label = paste("Total:", total)), 
  #           vjust = -1, size = 3.5, color = "black", position = position_dodge(width = 0.9)) + 
  theme(#axis.text.x = element_text(angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# bar acute ---------------------------------------------------------------

df = read.csv('Results_wide_association/Protein_association-acute.csv') %>% filter(p_value < 0.05/2923)
# assay = readRDS('H:/Prediction_model/DataPreparation/PRO_Data_Preprocessing/protein_assay.rds')%>% mutate(exposure = paste0('P',coding))
dfp = df %>% #filter(outcome == 'overall_pain') %>% 
  merge(.,assay,by='exposure') %>% group_by(Panel,outcome) %>% 
  summarise(n=n())
outcomes = c(#'overall_pain','acute_pain_status',
  'Generalized_Pain','acute_pain_status',
  'Head_and_Face_Pain',
  'Abdominal_Pain',
  'Neck_Shoulder_and_Back_Hip_and_Knee_Pain'
)

lable =  c(#'Overall pain','Acute pain',
  'Widespread pain','Overall pain',
  'Head and face pain',
  'Abdominal pain',
  'Musculoskeletal pain'
)

##### overall pain relative importance #####
colors <- c("Cardiometabolic" = "#374E5599",  # 色系中的第一种颜色
            "Cardiometabolic II" = "#374E554C",  # 第一种颜色的 70% 透明度
            "Inflammation" = "#DF8F4499",
            "Inflammation II" = "#DF8F444C",  # 第二种颜色的 70% 透明度
            "Neurology" = "#00A1D599",
            "Neurology II" = "#00A1D54C",  # 第三种颜色的 70% 透明度
            "Oncology" = "#B2474599",
            "Oncology II" = "#B247454C")  # 第四种颜色的 70% 透明度
total_counts <- dfp %>%filter(outcome %in% outcomes) %>% 
  mutate(outcome_label = factor(outcome, levels = outcomes, labels = lable)) %>% 
  group_by(outcome_label) %>%
  summarise(total = sum(n))

bar_acute = dfp %>%
  filter(outcome %in% outcomes) %>% 
  mutate(outcome = factor(outcome, levels = outcomes)) %>%
  mutate(outcome_label = factor(outcome, levels = outcomes, labels = lable)) %>%
  ggplot(aes(x = outcome_label, y = n, fill = Panel)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  ylim(0,180)+
  labs(x = "", y = "Count of proteins", title = "") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5) + 
  # geom_text(data = total_counts, aes(x = outcome_label, y = total, label = paste("Total:", total)), 
  #           vjust = -1, size = 3.5, color = "black", position = position_dodge(width = 0.9)) + 
  theme(#axis.text.x = element_text(angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#### plot
bar = ggarrange(
  protein_overlap, 
  ggarrange(
    bar_acute, bar_chronic,
    labels = c('Acute pain','Chronic pain'),
    nrow = 2, ncol = 1, common.legend = TRUE, legend = 'top'
  ), 
  # Labels removed as they are commented out in the original code
  nrow = 1, ncol = 2, common.legend = F,widths = c(0.35, 0.65)
)
ggarrange(
  bar, volcano_chronic, 
  # Labels removed as they are commented out in the original code
  nrow = 2, ncol = 1, common.legend = F,widths = c(0.5, 0.5)
)

export::graph2ppt(file = paste("results/02_Volcano_pritein.pptx", sep = ""), width = 12, height = 13.5, append = TRUE)
volcano_acute
export::graph2ppt(file = paste("results/01_Volcano_pritein.pptx", sep = ""), width = 12, height = 13.5/2, append = TRUE)

