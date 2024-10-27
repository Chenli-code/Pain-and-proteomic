load('D:/Desktop/PhD/10-学术项目/34-chronic_pain/Temp_folder2/Model.results_PROS_chronic.RData')


outcomes = c('chronic_pain_status',
             'Head_and_Face_Pain',
             'Abdominal_Pain','Neck_Shoulder_and_Back_Hip_and_Knee_Pain',
             'Generalized_Pain')
df_roc = data.frame()
df_auc = data.frame()
for (out in outcomes) {
  for (ss in c('_life.prob2','_ProRS.prob2','_pro.prob2','_simple.prob2','_simple_life.prob2')) {
    roc = rocit(results[[paste0(out,ss)]], as.numeric(results[[paste0(out,'_test')]][,out])-1)
    roc_t = data.frame(tpr = roc$TPR,fpr = roc$FPR)
    roc_t$model = ss;roc_t$outcome = out
    df_roc = rbind(df_roc,roc_t)
    auc = results[[paste0(out,gsub('prob2','auc2',ss))]]
    auc2 = results[[paste0(out,gsub('prob2','ci_auc2',ss))]]
    auc_t = data.frame(auc = auc,l = auc2[1],h=auc2[2])
    auc_t$model = ss;auc_t$outcome = out
    df_auc = rbind(df_auc,auc_t)
  }
}


outcomes = c(#'overall_pain','acute_pain_status',
  'chronic_pain_status',
  'Head_and_Face_Pain',
  'Abdominal_Pain',
  'Neck_Shoulder_and_Back_Hip_and_Knee_Pain',
  'Generalized_Pain')

lable =  c(#'Overall pain','Acute pain',
  'Chronic pain',
  'Head and face pain',
  'Abdominal pain',
  'Musculoskeletal pain',
  'Widespread pain')

#####  ROC ##### 
for (i in 1:5) {
  out = outcomes[i]
  label_out = lable[i]
  roc = df_roc %>% 
    mutate(model = factor(model,
                          levels = c('_life.prob2','_simple.prob2','_I-PS.prob2', '_simple_life.prob2','_pro.prob2'),
                          labels=c('CS','S-ProtS','I-ProtS','CS+S-ProtS','CS+I-ProtS')))
  custom_colors <- c("CS" = "#374E55FF", "I-ProtS" = "#DF8F44FF", "CS+I-ProtS" = "#00A1D5FF",
                     'S-ProtS' = '#925E9FFF','CS+S-ProtS'='#AD002AFF')
  
  auc = df_auc %>% filter(outcome == out)
  auc_text <- sprintf("CS: %.2f (%.2f-%.2f)\nS-ProtS: %.2f (%.2f-%.2f)\nI-ProtS: %.2f (%.2f-%.2f)\nCS+S-ProtS: %.2f (%.2f-%.2f)\nCS+I-ProtS: %.2f (%.2f-%.2f)",
                      auc$auc[1], auc$l[1], auc$h[1],
                      auc$auc[4], auc$l[4], auc$h[4],
                      auc$auc[2], auc$l[2], auc$h[2],
                      auc$auc[5], auc$l[5], auc$h[5],
                      auc$auc[3], auc$l[3], auc$h[3])
  #plot roc
  assign(paste0('roc_',out),
         roc %>%  filter(outcome == out) %>% 
           ggplot(aes(x = fpr, y = tpr,color = model)) +
           geom_line(linewidth = 1.1) + 
           geom_abline(linetype = "dashed") +
           labs(x = "False Positive Rate (FPR)",
                y = "True Positive Rate (TPR)",
                color = 'Model',
                title = label_out) +
           scale_color_manual(values = custom_colors) +
           theme_bw()+
           theme(legend.position = "none",
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())+
           annotate("text", x = 0.95, y = 0.05, label = auc_text, size = 3, color = "black", hjust = 1, vjust = 0)
  )
  
}

# library(ggpubr)
roc = ggarrange(roc_chronic_pain_status,roc_Head_and_Face_Pain,roc_Abdominal_Pain,
                roc_Neck_Shoulder_and_Back_Hip_and_Knee_Pain,
                # roc_Generalized_Pain,
                ncol = 1, nrow = 4)



##### calibration curve ##### 

cali_fun = function(benchp.prob_positive,ob){
  d = data.frame(
    prob = benchp.prob_positive,
    observed = as.numeric(as.character(ob)),
    bin = cut(benchp.prob_positive, 
              breaks=unique(as.numeric(quantile(benchp.prob_positive,seq(0,1,0.1)))), 
              include.lowest=TRUE, right=TRUE)) %>%
    group_by(bin) %>%
    summarise(
      mean_predicted = mean(prob, na.rm = TRUE),
      observed_freq = mean(observed, na.rm = TRUE))
  return(d)
}

for (i in 1:5) {
  out = outcomes[i]
  label_out = lable[i]
  # Calculate observed frequency for each bin
  cali.life <- cali_fun(results[[paste0(out,'_life.prob2')]], results[[paste0(out,'_test')]][,out])
  cali.IPS <- cali_fun(results[[paste0(out,'_ProRS.prob2')]], results[[paste0(out,'_test')]][,out])
  cali.pro <- cali_fun(results[[paste0(out,'_pro.prob2')]], results[[paste0(out,'_test')]][,out])
  cali.simple <- cali_fun(results[[paste0(out,'_simple.prob2')]], results[[paste0(out,'_test')]][,out])
  cali.simpletop <- cali_fun(results[[paste0(out,'_simple_life.prob2')]], results[[paste0(out,'_test')]][,out])
  
  calibration_stats = rbind(cali.life,cali.IPS,cali.pro,cali.simple,cali.simpletop)
  calibration_stats$Model = c(rep('CS',nrow(cali.life)),
                              rep('I-ProtS',nrow(cali.IPS)),
                              rep('CS+I-ProtS',nrow(cali.pro)),
                              rep('S-ProtS',nrow(cali.simple)),
                              rep('CS+S-ProtS',nrow(cali.simpletop)))
  
  library(CalibrationCurves)
  p = valProbggplot(results[[paste0(out,'_pro.prob2')]], 
                    as.numeric(as.character(results[[paste0(out,'_test')]][,out])))
  # 提取校准结果及置信区间
  calibration_intercept <- p$Calibration$Intercept[1]
  calibration_ci_lower <- p$Calibration$Intercept[2]
  calibration_ci_upper <- p$Calibration$Intercept[3]
  
  # 格式化校准结果和置信区间
  calibration_text <- sprintf("Intercept: %.2f (%.2f, %.2f)", 
                              calibration_intercept, calibration_ci_lower, calibration_ci_upper)
  
  # 提取校准结果及置信区间
  calibration_Slope <- p$Calibration$Slope[1]
  Slope_ci_lower <- p$Calibration$Slope[2]
  Slope_ci_upper <- p$Calibration$Slope[3]
  
  # 格式化校准结果和置信区间
  Slope_text <- sprintf("Slope: %.2f (%.2f, %.2f)", 
                        calibration_Slope, Slope_ci_lower, Slope_ci_upper)
  
  assign(paste0('cali_',out),
         calibration_stats %>% 
           mutate(Model = factor(Model,levels=c('CS','I-ProtS','CS+I-ProtS','S-ProtS','CS+S-ProtS'))) %>% 
           ggplot(aes(x = mean_predicted, y = observed_freq,color=Model)) +
           geom_point(size = 2)+
           geom_line(linewidth = 0.8) +
           geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
           xlim(0, NA) +ylim(0, NA) +
           labs(x = "Predicted Probability",
                y = "Observed event rate",
                color = 'Model',
                title = str_to_title(gsub('status','',gsub('_',' ',out)))) +#
           theme_bw()+
           # ggsci::scale_color_jama()+
           
           scale_color_manual(values = custom_colors) +
           theme(legend.position = "none",
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())+
           annotate("text", x = -Inf, y = Inf, 
                    label = paste("Calibration\n",calibration_text,'\n',Slope_text),
                    size = 3, color = "black", hjust = -0.1, vjust = 1.1)
  )
}


library(ggpubr)
cali = ggarrange(cali_chronic_pain_status,cali_Head_and_Face_Pain,cali_Abdominal_Pain,
                 cali_Neck_Shoulder_and_Back_Hip_and_Knee_Pain,
                 # cali_Generalized_Pain,
                 ncol = 1, nrow = 4)
# cali = cali_chronic_pain_status

pp = ggarrange(roc,cali,ncol = 2, nrow = 1)
pp
export::graph2ppt(file=paste0('results/calibration_all2.pptx'),width = 8,height =15, append = T)

pp = ggarrange(roc_Generalized_Pain,
               ggarrange(cali_Generalized_Pain,cali_Generalized_Pain,ncol = 1, nrow = 2),
               widths = c(0.6,0.4),
               ncol = 2, nrow = 1)
pp
export::graph2ppt(file=paste0('results/calibration.pptx'),width = 8,height =5, append = T)



for (out in outcomes) {
  
  # simple_pro = pro_im[which(pro_im$outcome == out),]$pro
  ##### standardized net benefit curve ##### 
  library('rmda')
  # The author sets seed at 123; I do the same throughout this document to replicate the results
  set.seed(123)
  
  test = results[[paste0(out,'_test')]]
  test[,out] = as.numeric(as.character(test[,out]))
  # the below model is run with default settings and 50 bootstraps
  set.seed(123)
  life.model <- decision_curve(as.formula(paste(out,"~", paste(c(life_style), collapse = " + "))), #fitting a logistic model
                               data = test[,c(out,life_style)],
                               # study.design = "cohort",
                               policy = "opt-in", #default
                               bootstraps = 50)
  set.seed(123)
  I_PS.model <- decision_curve(as.formula(paste(out,"~", paste(c(paste0(out,'_ProRS')), collapse = " + "))), #fitting a logistic model
                               data = test[,c(out,paste0(out,'_ProRS'))],
                               # study.design = "cohort",
                               policy = "opt-in", #default
                               bootstraps = 50)
  set.seed(123)
  pro.model <- decision_curve(as.formula(paste(out,"~", paste(c(life_style,paste0(out,'_ProRS')), collapse = " + "))), #fitting a logistic model
                              data = test[,c(out,life_style,paste0(out,'_ProRS'))],
                              # study.design = "cohort",
                              policy = "opt-in", #default
                              bootstraps = 50)
  
  set.seed(123)
  simple.model <- decision_curve(as.formula(paste(out,"~", paste(c(paste0(out, '_simpleRS')), collapse = " + "))), #fitting a logistic model
                                 data = test[,c(out,paste0(out, '_simpleRS'))],
                                 # study.design = "cohort",
                                 policy = "opt-in", #default
                                 bootstraps = 50)
  
  set.seed(123)
  simpletop.model <- decision_curve(as.formula(paste(out,"~", paste(c(life_style,paste0(out, '_simpleRS')), collapse = " + "))), #fitting a logistic model
                                    data = test[,c(out,life_style,paste0(out, '_simpleRS'))],
                                    # study.design = "cohort",
                                    policy = "opt-in", #default
                                    bootstraps = 50)
  
  
  # png("results/net benefit curve.png")
  net <- plot_decision_curve( 
    list(life.model,I_PS.model,pro.model,simple.model,simpletop.model), 
    confidence.intervals = F,cost.benefit.xlab = F,standardize = T,cost.benefit.axis = FALSE,
    curve.names = c('CS','I-PS','CS+I-PS','S-PS','CS+S-PS'),
    col = c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#925E9FFF", "#AD002AFF"),
    legend.position = 'none')
  
  # export::graph2ppt(file=paste0('results/benefit.pptx'),width = 5.8,height =6.2, append = T)
  export::graph2ppt(file=paste0('results/benefit3_',out,'.pptx'),width = 4.7,height =4, append = T)
  
}
