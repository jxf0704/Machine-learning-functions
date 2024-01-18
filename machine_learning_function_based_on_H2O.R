library(WGCNA)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(h2o)
h2o.init()
library(pROC)
library(dplyr)

external_validation_h2o<-function(motif_data_frame=motif_table,including_motif=motif_name,response,training_unit,testing_unit,time_point="P2",last_time="no"){
  if(last_time=="no"){
    training_data_frame<-motif_data_frame[motif_data_frame$unit %in% training_unit & motif_data_frame$Time_point %in% time_point,]
    testing_data_frame<-motif_data_frame[motif_data_frame$unit %in% testing_unit & motif_data_frame$Time_point %in% time_point,]
    training_data_frame<-training_data_frame[,c(including_motif,response)]
    testing_data_frame<-testing_data_frame[,c(including_motif,response)]
    #training_data_frame<-na.omit(training_data_frame)
    #testing_data_frame<-na.omit(testing_data_frame)
    training_data_frame[,response]<-as.factor(training_data_frame[,response])
    testing_data_frame[,response]<-as.factor(testing_data_frame[,response])
    target<-response
    features<-setdiff(colnames(training_data_frame), target)
    training_h2o<-as.h2o(training_data_frame)
    testing_h2o<-as.h2o(testing_data_frame)
  }
  
  if(last_time=="yes"){
    training_data_frame<-motif_data_frame[motif_data_frame$unit %in% training_unit & motif_data_frame$last_time==1,]
    testing_data_frame<-motif_data_frame[motif_data_frame$unit %in% testing_unit & motif_data_frame$last_time==1,]
    training_data_frame<-training_data_frame[,c(including_motif,response)]
    testing_data_frame<-testing_data_frame[,c(including_motif,response)]
    training_data_frame<-na.omit(training_data_frame)
    testing_data_frame<-na.omit(testing_data_frame)
    training_data_frame[,response]<-as.factor(training_data_frame[,response])
    testing_data_frame[,response]<-as.factor(testing_data_frame[,response])
    target<-response
    features<-setdiff(colnames(training_data_frame), target)
    training_h2o<-as.h2o(training_data_frame)
    testing_h2o<-as.h2o(testing_data_frame)
  }
  result<-list()
  result[["training_h2o"]]<-training_h2o
  result[["testing_h2o"]]<-testing_h2o
  result[["target"]]<-target
  result[["features"]]<-features
  result[["training_table"]]<-motif_data_frame[motif_data_frame$unit %in% training_unit & motif_data_frame$Time_point %in% time_point,]
  result[["testing_table"]]<-motif_data_frame[motif_data_frame$unit %in% testing_unit & motif_data_frame$Time_point %in% time_point,]
  
  return(result)
}

h2o_glm_construction<-function(h2o_data_list,n_seed=12345,project_name,max_models=5,split_type="Stratified",cutoff_field="accuracy",n_fold=5){
  h2o_model = h2o.glm(x = h2o_data_list$features,
                      y =h2o_data_list$target,
                      training_frame = h2o_data_list$training_h2o,
                      nfolds = n_fold,                        # 5-fold Cross-Validation
                      #max_models = max_models,                   # Max number of models
                      #stopping_metric = "logloss",       # Metric to optimize
                      #project_name = project_name,  # Specify a name so you can add more models later
                      #sort_metric = "logloss",
                      seed = n_seed,
                      fold_assignment=split_type,
                      keep_cross_validation_fold_assignment = TRUE,
                      keep_cross_validation_predictions= TRUE)
  lead_model<-h2o_model
  model_data<-lead_model@model
  used_model<-paste(model_data$model_summary$link,model_data$model_summary$regularization,model_data$model_summary$number_of_predictors_total,model_data$model_summary$number_of_active_predictors,sep="_")
  training_AUC<-model_data$training_metrics@metrics$AUC
  training_threshold<-as.data.frame(model_data$training_metrics@metrics$thresholds_and_metric_scores)
  training_threshold<-training_threshold[order(training_threshold[,cutoff_field],decreasing=T),]
  top_training_cutoff_field<-as.data.frame(training_threshold[1,])
  top_training_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_training_cutoff_field$group<-"training"
  top_training_cutoff_field$auc<-training_AUC
  
  
  CV_auc<-model_data$cross_validation_metrics@metrics$AUC
  cv_threshold<-as.data.frame(model_data$cross_validation_metrics@metrics$thresholds_and_metric_scores)
  cv_threshold<-cv_threshold[order(cv_threshold[,cutoff_field],decreasing=T),]
  
  top_cv_cutoff_field<-as.data.frame(cv_threshold[1,])
  top_cv_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_cv_cutoff_field$group<-"CV"
  top_cv_cutoff_field$auc<-CV_auc
  
  cv_prediction<-as.data.frame(h2o.getFrame(model_data[["cross_validation_holdout_predictions_frame_id"]][["name"]]))
  #cv_prediction$PCR_status<-ky087_motif_annotation_fragment_220_400_arm_T4_103$PCR_status
  #roc(cv_prediction$PCR_status,cv_prediction$p1)
  cv_prediction<-as.data.frame(cv_prediction)
  cv_prediction$Sample_ID<-h2o_data_list$training_table$Sample_ID
  cv_prediction$SampleID<-h2o_data_list$training_table$SampleID
  cv_prediction$response<-h2o_data_list$training_table[,h2o_data_list$target]
  feature_contribution<-as.data.frame(h2o.varimp(lead_model))
  
  
  perf <- h2o.performance(lead_model, h2o_data_list$testing_h2o)
  testing_auc<-h2o.auc(perf)
  pred <- h2o.predict(object = lead_model, newdata = h2o_data_list$testing_h2o)
  pred<-as.data.frame(pred)
  pred$Sample_ID<-h2o_data_list$testing_table$Sample_ID
  pred$SampleID<-h2o_data_list$testing_table$SampleID
  pred$response<-h2o_data_list$testing_table[,h2o_data_list$target]
  
  testing_threshold<-as.data.frame(perf@metrics$thresholds_and_metric_scores)
  testing_threshold<-testing_threshold[order(testing_threshold[,cutoff_field],decreasing=T),]
  
  top_testing_cutoff_field<-as.data.frame(testing_threshold[1,])
  top_testing_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_testing_cutoff_field$group<-"testing"
  top_testing_cutoff_field$auc<-testing_auc
  
  cv_prediction$detail_group<-h2o_data_list$training_table$Detail_Group
  cv_prediction$three_group<-cv_prediction$detail_group
  #cv_prediction$three_group[cv_prediction$three_group=="Healthy-Healthy"]<-"Healthy"
  #cv_prediction$three_group[cv_prediction$three_group=="Disease-lung nodules"]<-"Nodule"
  #cv_prediction$three_group[cv_prediction$three_group=="Lung"]<-"Cancer"
  cv_prediction$train_group<-"training"
  
  pred$detail_group<-h2o_data_list$testing_table$Detail_Group
  pred$three_group<-pred$detail_group
  #pred$three_group[pred$three_group=="Healthy-Healthy"]<-"Healthy"
  #pred$three_group[pred$three_group=="Disease-lung nodules"]<-"Nodule"
  #pred$three_group[pred$three_group=="Lung"]<-"Cancer"
  pred$train_group<-"testing"
  total_prediction_table<-rbind(cv_prediction,pred)
  
  result<-list()
  result_table<-rbind(top_training_cutoff_field,top_cv_cutoff_field,top_testing_cutoff_field)
  result_table$project_name<-project_name
  result_table$model<-used_model
  result[["auc_table"]]<-result_table
  result[["lead_model"]]<-lead_model
  result[["model_data"]]<-model_data
  result[["cv_prediction"]]<-cv_prediction
  result[["feature_contribution"]]<-feature_contribution
  result[["test_prediction"]]<-pred
  result[["total_prediction"]]<-total_prediction_table
  return(result)
}

h2o_gbm_construction<-function(h2o_data_list,n_seed=12345,project_name,max_models=5,split_type="Stratified",cutoff_field="accuracy",n_fold=5){
  h2o_model = h2o.gbm(x = h2o_data_list$features,
                      y =h2o_data_list$target,
                      training_frame = h2o_data_list$training_h2o,
                      nfolds = n_fold,                        # 5-fold Cross-Validation
                      #max_models = max_models,                   # Max number of models
                      #stopping_metric = "logloss",       # Metric to optimize
                      #project_name = project_name,  # Specify a name so you can add more models later
                      #sort_metric = "logloss",
                      seed = n_seed,
                      fold_assignment=split_type,
                      keep_cross_validation_fold_assignment = TRUE,
                      keep_cross_validation_predictions= TRUE)
  lead_model<-h2o_model
  model_data<-lead_model@model
  used_model<-paste(model_data$model_summary$number_of_trees,model_data$model_summary$mean_depth,model_data$model_summary$min_leaves,model_data$model_summary$max_leaves,sep="_")
  training_AUC<-model_data$training_metrics@metrics$AUC
  training_threshold<-as.data.frame(model_data$training_metrics@metrics$thresholds_and_metric_scores)
  training_threshold<-training_threshold[order(training_threshold[,cutoff_field],decreasing=T),]
  top_training_cutoff_field<-as.data.frame(training_threshold[1,])
  top_training_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_training_cutoff_field$group<-"training"
  top_training_cutoff_field$auc<-training_AUC
  
  
  CV_auc<-model_data$cross_validation_metrics@metrics$AUC
  cv_threshold<-as.data.frame(model_data$cross_validation_metrics@metrics$thresholds_and_metric_scores)
  cv_threshold<-cv_threshold[order(cv_threshold[,cutoff_field],decreasing=T),]
  
  top_cv_cutoff_field<-as.data.frame(cv_threshold[1,])
  top_cv_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_cv_cutoff_field$group<-"CV"
  top_cv_cutoff_field$auc<-CV_auc
  
  cv_prediction<-as.data.frame(h2o.getFrame(model_data[["cross_validation_holdout_predictions_frame_id"]][["name"]]))
  #cv_prediction$PCR_status<-ky087_motif_annotation_fragment_220_400_arm_T4_103$PCR_status
  #roc(cv_prediction$PCR_status,cv_prediction$p1)
  p<-cv_prediction$p1
  cv_prediction<-as.data.frame(cv_prediction)
  cv_prediction$Sample_ID<-h2o_data_list$training_table$Sample_ID
  cv_prediction$SampleID<-h2o_data_list$training_table$SampleID
  cv_prediction$response<-h2o_data_list$training_table[,h2o_data_list$target]
  feature_contribution<-as.data.frame(h2o.varimp(lead_model))
  
  
  perf <- h2o.performance(lead_model, h2o_data_list$testing_h2o)
  testing_auc<-h2o.auc(perf)
  pred <- h2o.predict(object = lead_model, newdata = h2o_data_list$testing_h2o)
  pred<-as.data.frame(pred)
  pred$Sample_ID<-h2o_data_list$testing_table$Sample_ID
  pred$SampleID<-h2o_data_list$testing_table$SampleID
  pred$response<-h2o_data_list$testing_table[,h2o_data_list$target]
  
  testing_threshold<-as.data.frame(perf@metrics$thresholds_and_metric_scores)
  testing_threshold<-testing_threshold[order(testing_threshold[,cutoff_field],decreasing=T),]
  
  top_testing_cutoff_field<-as.data.frame(testing_threshold[1,])
  top_testing_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_testing_cutoff_field$group<-"testing"
  top_testing_cutoff_field$auc<-testing_auc
  cv_prediction$detail_group<-h2o_data_list$training_table$Detail_Group
  cv_prediction$three_group<-cv_prediction$detail_group
  #cv_prediction$three_group[cv_prediction$three_group=="Healthy-Healthy"]<-"Healthy"
  #cv_prediction$three_group[cv_prediction$three_group=="Disease-lung nodules"]<-"Nodule"
  #cv_prediction$three_group[cv_prediction$three_group=="Lung"]<-"Cancer"
  cv_prediction$train_group<-"training"
  
  pred$detail_group<-h2o_data_list$testing_table$Detail_Group
  pred$three_group<-pred$detail_group
  #pred$three_group[pred$three_group=="Healthy-Healthy"]<-"Healthy"
  #pred$three_group[pred$three_group=="Disease-lung nodules"]<-"Nodule"
  #pred$three_group[pred$three_group=="Lung"]<-"Cancer"
  pred$train_group<-"testing"
  total_prediction_table<-rbind(cv_prediction,pred)
  
  result<-list()
  result_table<-rbind(top_training_cutoff_field,top_cv_cutoff_field,top_testing_cutoff_field)
  result_table$project_name<-project_name
  result_table$model<-used_model
  result[["auc_table"]]<-result_table
  result[["lead_model"]]<-lead_model
  result[["model_data"]]<-model_data
  result[["cv_prediction"]]<-cv_prediction
  result[["feature_contribution"]]<-feature_contribution
  result[["test_prediction"]]<-pred
  result[["total_prediction"]]<-total_prediction_table
  return(result)
}

h2o_rf_construction<-function(h2o_data_list,n_seed=12345,project_name,max_models=5,split_type="Stratified",cutoff_field="accuracy",n_fold=5){
  h2o_model = h2o.randomForest(x = h2o_data_list$features,
                               y =h2o_data_list$target,
                               training_frame = h2o_data_list$training_h2o,
                               nfolds = n_fold,                        # 5-fold Cross-Validation
                               #max_models = max_models,                   # Max number of models
                               #stopping_metric = "logloss",       # Metric to optimize
                               #project_name = project_name,  # Specify a name so you can add more models later
                               #sort_metric = "logloss",
                               seed = n_seed,
                               fold_assignment=split_type,
                               keep_cross_validation_fold_assignment = TRUE,
                               keep_cross_validation_predictions= TRUE)
  lead_model<-h2o_model
  model_data<-lead_model@model
  used_model<-paste(model_data$model_summary$number_of_trees,model_data$model_summary$mean_depth,model_data$model_summary$min_leaves,model_data$model_summary$max_leaves,sep="_")
  training_AUC<-model_data$training_metrics@metrics$AUC
  training_threshold<-as.data.frame(model_data$training_metrics@metrics$thresholds_and_metric_scores)
  training_threshold<-training_threshold[order(training_threshold[,cutoff_field],decreasing=T),]
  top_training_cutoff_field<-as.data.frame(training_threshold[1,])
  top_training_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_training_cutoff_field$group<-"training"
  top_training_cutoff_field$auc<-training_AUC
  
  
  CV_auc<-model_data$cross_validation_metrics@metrics$AUC
  cv_threshold<-as.data.frame(model_data$cross_validation_metrics@metrics$thresholds_and_metric_scores)
  cv_threshold<-cv_threshold[order(cv_threshold[,cutoff_field],decreasing=T),]
  
  top_cv_cutoff_field<-as.data.frame(cv_threshold[1,])
  top_cv_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_cv_cutoff_field$group<-"CV"
  top_cv_cutoff_field$auc<-CV_auc
  
  cv_prediction<-as.data.frame(h2o.getFrame(model_data[["cross_validation_holdout_predictions_frame_id"]][["name"]]))
  #cv_prediction$PCR_status<-ky087_motif_annotation_fragment_220_400_arm_T4_103$PCR_status
  #roc(cv_prediction$PCR_status,cv_prediction$p1)
  p<-cv_prediction$p1
  cv_prediction<-as.data.frame(cv_prediction)
  cv_prediction$Sample_ID<-h2o_data_list$training_table$Sample_ID
  cv_prediction$SampleID<-h2o_data_list$training_table$SampleID
  cv_prediction$response<-h2o_data_list$training_table[,h2o_data_list$target]
  feature_contribution<-as.data.frame(h2o.varimp(lead_model))
  
  
  perf <- h2o.performance(lead_model, h2o_data_list$testing_h2o)
  testing_auc<-h2o.auc(perf)
  pred <- h2o.predict(object = lead_model, newdata = h2o_data_list$testing_h2o)
  pred<-as.data.frame(pred)
  pred$Sample_ID<-h2o_data_list$testing_table$Sample_ID
  pred$SampleID<-h2o_data_list$testing_table$SampleID
  pred$response<-h2o_data_list$testing_table[,h2o_data_list$target]
  
  testing_threshold<-as.data.frame(perf@metrics$thresholds_and_metric_scores)
  testing_threshold<-testing_threshold[order(testing_threshold[,cutoff_field],decreasing=T),]
  
  top_testing_cutoff_field<-as.data.frame(testing_threshold[1,])
  top_testing_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_testing_cutoff_field$group<-"testing"
  top_testing_cutoff_field$auc<-testing_auc
  
  cv_prediction$detail_group<-h2o_data_list$training_table$Detail_Group
  cv_prediction$three_group<-cv_prediction$detail_group
  #cv_prediction$three_group[cv_prediction$three_group=="Healthy-Healthy"]<-"Healthy"
  #cv_prediction$three_group[cv_prediction$three_group=="Disease-lung nodules"]<-"Nodule"
  #cv_prediction$three_group[cv_prediction$three_group=="Lung"]<-"Cancer"
  cv_prediction$train_group<-"training"
  
  pred$detail_group<-h2o_data_list$testing_table$Detail_Group
  pred$three_group<-pred$detail_group
  #pred$three_group[pred$three_group=="Healthy-Healthy"]<-"Healthy"
  #pred$three_group[pred$three_group=="Disease-lung nodules"]<-"Nodule"
  #pred$three_group[pred$three_group=="Lung"]<-"Cancer"
  pred$train_group<-"testing"
  total_prediction_table<-rbind(cv_prediction,pred)
  
  result<-list()
  result_table<-rbind(top_training_cutoff_field,top_cv_cutoff_field,top_testing_cutoff_field)
  result_table$project_name<-project_name
  result_table$model<-used_model
  result[["auc_table"]]<-result_table
  result[["lead_model"]]<-lead_model
  result[["model_data"]]<-model_data
  result[["cv_prediction"]]<-cv_prediction
  result[["feature_contribution"]]<-feature_contribution
  result[["test_prediction"]]<-pred
  result[["total_prediction"]]<-total_prediction_table
  return(result)
}

h2o_dl_construction<-function(h2o_data_list,n_seed=12345,project_name,max_models=5,split_type="Stratified",cutoff_field="accuracy",n_fold=5,repro=TRUE){
  h2o_model = h2o.deeplearning(x = h2o_data_list$features,
                               y =h2o_data_list$target,
                               training_frame = h2o_data_list$training_h2o,
                               nfolds = n_fold,                        # 5-fold Cross-Validation
                               #max_models = max_models,                   # Max number of models
                               #stopping_metric = "logloss",       # Metric to optimize
                               #project_name = project_name,  # Specify a name so you can add more models later
                               #sort_metric = "logloss",
                               seed = n_seed,
                               fold_assignment=split_type,
                               keep_cross_validation_fold_assignment = TRUE,
                               epochs = 50,
                               reproducible=repro,
                               keep_cross_validation_predictions= TRUE)
  lead_model<-h2o_model
  model_data<-lead_model@model
  used_model<-"deep_learning"
  training_AUC<-model_data$training_metrics@metrics$AUC
  training_threshold<-as.data.frame(model_data$training_metrics@metrics$thresholds_and_metric_scores)
  training_threshold<-training_threshold[order(training_threshold[,cutoff_field],decreasing=T),]
  top_training_cutoff_field<-as.data.frame(training_threshold[1,])
  top_training_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_training_cutoff_field$group<-"training"
  top_training_cutoff_field$auc<-training_AUC
  
  
  CV_auc<-model_data$cross_validation_metrics@metrics$AUC
  cv_threshold<-as.data.frame(model_data$cross_validation_metrics@metrics$thresholds_and_metric_scores)
  cv_threshold<-cv_threshold[order(cv_threshold[,cutoff_field],decreasing=T),]
  
  top_cv_cutoff_field<-as.data.frame(cv_threshold[1,])
  top_cv_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_cv_cutoff_field$group<-"CV"
  top_cv_cutoff_field$auc<-CV_auc
  
  cv_prediction<-as.data.frame(h2o.getFrame(model_data[["cross_validation_holdout_predictions_frame_id"]][["name"]]))
  #cv_prediction$PCR_status<-ky087_motif_annotation_fragment_220_400_arm_T4_103$PCR_status
  #roc(cv_prediction$PCR_status,cv_prediction$p1)
  p<-cv_prediction$p1
  cv_prediction<-as.data.frame(cv_prediction)
  cv_prediction$Sample_ID<-h2o_data_list$training_table$Sample_ID
  cv_prediction$SampleID<-h2o_data_list$training_table$SampleID
  cv_prediction$response<-h2o_data_list$training_table[,h2o_data_list$target]
  feature_contribution<-as.data.frame(h2o.varimp(lead_model))
  
  
  perf <- h2o.performance(lead_model, h2o_data_list$testing_h2o)
  testing_auc<-h2o.auc(perf)
  pred <- h2o.predict(object = lead_model, newdata = h2o_data_list$testing_h2o)
  pred<-as.data.frame(pred)
  pred$Sample_ID<-h2o_data_list$testing_table$Sample_ID
  pred$SampleID<-h2o_data_list$testing_table$SampleID
  pred$response<-h2o_data_list$testing_table[,h2o_data_list$target]
  
  testing_threshold<-as.data.frame(perf@metrics$thresholds_and_metric_scores)
  testing_threshold<-testing_threshold[order(testing_threshold[,cutoff_field],decreasing=T),]
  
  top_testing_cutoff_field<-as.data.frame(testing_threshold[1,])
  top_testing_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_testing_cutoff_field$group<-"testing"
  top_testing_cutoff_field$auc<-testing_auc
  
  cv_prediction$detail_group<-h2o_data_list$training_table$Detail_Group
  cv_prediction$three_group<-cv_prediction$detail_group
  #cv_prediction$three_group[cv_prediction$three_group=="Healthy-Healthy"]<-"Healthy"
  #cv_prediction$three_group[cv_prediction$three_group=="Disease-lung nodules"]<-"Nodule"
  #cv_prediction$three_group[cv_prediction$three_group=="Lung"]<-"Cancer"
  cv_prediction$train_group<-"training"
  
  pred$detail_group<-h2o_data_list$testing_table$Detail_Group
  pred$three_group<-pred$detail_group
  #pred$three_group[pred$three_group=="Healthy-Healthy"]<-"Healthy"
  #pred$three_group[pred$three_group=="Disease-lung nodules"]<-"Nodule"
  #pred$three_group[pred$three_group=="Lung"]<-"Cancer"
  pred$train_group<-"testing"
  total_prediction_table<-rbind(cv_prediction,pred)
  
  result<-list()
  result_table<-rbind(top_training_cutoff_field,top_cv_cutoff_field,top_testing_cutoff_field)
  result_table$project_name<-project_name
  result_table$model<-used_model
  result[["auc_table"]]<-result_table
  result[["lead_model"]]<-lead_model
  result[["model_data"]]<-model_data
  result[["cv_prediction"]]<-cv_prediction
  result[["feature_contribution"]]<-feature_contribution
  result[["test_prediction"]]<-pred
  result[["total_prediction"]]<-total_prediction_table
  return(result)
}

h2o_bayes_construction<-function(h2o_data_list,n_seed=12345,project_name,max_models=5,split_type="Stratified",cutoff_field="accuracy",n_fold=5,repro=TRUE){
  h2o_model = h2o.naiveBayes(x = h2o_data_list$features,
                             y =h2o_data_list$target,
                             training_frame = h2o_data_list$training_h2o,
                             nfolds = n_fold,                        # 5-fold Cross-Validation
                             #max_models = max_models,                   # Max number of models
                             #stopping_metric = "logloss",       # Metric to optimize
                             #project_name = project_name,  # Specify a name so you can add more models later
                             #sort_metric = "logloss",
                             seed = n_seed,
                             fold_assignment=split_type,
                             keep_cross_validation_fold_assignment = TRUE,
                             #epochs = 50,
                             #reproducible=repro,
                             keep_cross_validation_predictions= TRUE)
  lead_model<-h2o_model
  model_data<-lead_model@model
  used_model<-"deep_learning"
  training_AUC<-model_data$training_metrics@metrics$AUC
  training_threshold<-as.data.frame(model_data$training_metrics@metrics$thresholds_and_metric_scores)
  training_threshold<-training_threshold[order(training_threshold[,cutoff_field],decreasing=T),]
  top_training_cutoff_field<-as.data.frame(training_threshold[1,])
  top_training_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_training_cutoff_field$group<-"training"
  top_training_cutoff_field$auc<-training_AUC
  
  
  CV_auc<-model_data$cross_validation_metrics@metrics$AUC
  cv_threshold<-as.data.frame(model_data$cross_validation_metrics@metrics$thresholds_and_metric_scores)
  cv_threshold<-cv_threshold[order(cv_threshold[,cutoff_field],decreasing=T),]
  
  top_cv_cutoff_field<-as.data.frame(cv_threshold[1,])
  top_cv_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_cv_cutoff_field$group<-"CV"
  top_cv_cutoff_field$auc<-CV_auc
  
  cv_prediction<-as.data.frame(h2o.getFrame(model_data[["cross_validation_holdout_predictions_frame_id"]][["name"]]))
  #cv_prediction$PCR_status<-ky087_motif_annotation_fragment_220_400_arm_T4_103$PCR_status
  #roc(cv_prediction$PCR_status,cv_prediction$p1)
  p<-cv_prediction$p1
  cv_prediction<-as.data.frame(cv_prediction)
  cv_prediction$Sample_ID<-h2o_data_list$training_table$Sample_ID
  cv_prediction$SampleID<-h2o_data_list$training_table$SampleID
  cv_prediction$response<-h2o_data_list$training_table[,h2o_data_list$target]
  feature_contribution<-as.data.frame(h2o.varimp(lead_model))
  
  
  perf <- h2o.performance(lead_model, h2o_data_list$testing_h2o)
  testing_auc<-h2o.auc(perf)
  pred <- h2o.predict(object = lead_model, newdata = h2o_data_list$testing_h2o)
  pred<-as.data.frame(pred)
  pred$Sample_ID<-h2o_data_list$testing_table$Sample_ID
  pred$SampleID<-h2o_data_list$testing_table$SampleID
  pred$response<-h2o_data_list$testing_table[,h2o_data_list$target]
  
  testing_threshold<-as.data.frame(perf@metrics$thresholds_and_metric_scores)
  testing_threshold<-testing_threshold[order(testing_threshold[,cutoff_field],decreasing=T),]
  
  top_testing_cutoff_field<-as.data.frame(testing_threshold[1,])
  top_testing_cutoff_field$cutoff_standard<-paste("max",cutoff_field,sep="_")
  top_testing_cutoff_field$group<-"testing"
  top_testing_cutoff_field$auc<-testing_auc
  
  #cv_prediction$detail_group<-h2o_data_list$training_table$Detail_Group
  #cv_prediction$three_group<-cv_prediction$detail_group
  #cv_prediction$three_group[cv_prediction$three_group=="Healthy-Healthy"]<-"Healthy"
  #cv_prediction$three_group[cv_prediction$three_group=="Disease-lung nodules"]<-"Nodule"
  #cv_prediction$three_group[cv_prediction$three_group=="Lung"]<-"Cancer"
  cv_prediction$train_group<-"training"
  
  #pred$detail_group<-h2o_data_list$testing_table$Detail_Group
  #pred$three_group<-pred$detail_group
  #pred$three_group[pred$three_group=="Healthy-Healthy"]<-"Healthy"
  #pred$three_group[pred$three_group=="Disease-lung nodules"]<-"Nodule"
  #pred$three_group[pred$three_group=="Lung"]<-"Cancer"
  pred$train_group<-"testing"
  
  total_prediction_table<-rbind(cv_prediction,pred)
  
  
  result<-list()
  result_table<-rbind(top_training_cutoff_field,top_cv_cutoff_field,top_testing_cutoff_field)
  result_table$project_name<-project_name
  result_table$model<-used_model
  result[["auc_table"]]<-result_table
  result[["lead_model"]]<-lead_model
  result[["model_data"]]<-model_data
  result[["cv_prediction"]]<-cv_prediction
  result[["feature_contribution"]]<-feature_contribution
  result[["test_prediction"]]<-pred
  result[["total_prediction"]]<-total_prediction_table
  return(result)
}

sum_is_na<-function(x){
  res<-sum(is.na(x))
  return(res)
}

PCA_feature_generate<-function(feature_data_frame=feature_data,feature_list=study_feature$frag_arm,trainingunit="training",testingunit="testing",pc_number=50){
  training_data_frame<-feature_data_frame[feature_data_frame$unit==trainingunit,feature_list]
  
  testing_data_frame<-feature_data_frame[feature_data_frame$unit==testingunit,feature_list]
  is.na_arm<-as.data.frame(apply(training_data_frame,2,sum_is_na))
  colnames(is.na_arm)[[1]]<-"NA_count"
  keep_arm<-rownames(is.na_arm)[is.na_arm$NA_count!=nrow(training_data_frame)]
  training_data_frame<-training_data_frame[,keep_arm]
  training_data_frame[is.na(training_data_frame)]<-0
  testing_data_frame<-testing_data_frame[,keep_arm]
  testing_data_frame[is.na(testing_data_frame)]<-0
  
  res.pca<- PCA(training_data_frame,ncp=pc_number, scale.unit =T, graph = FALSE)
  pca_pred<-predict(res.pca, testing_data_frame)
  
  training_data.frame<-feature_data_frame[feature_data_frame$unit==trainingunit,]
  testing_data.frame<-feature_data_frame[feature_data_frame$unit==testingunit,]
  training_data.frame<-cbind(training_data.frame,res.pca$ind$coord)
  testing_data.frame<-cbind(testing_data.frame,pca_pred$coord)
  final_data_frame<-rbind(training_data.frame,testing_data.frame)
  pca_result<-list()
  pca_result[["training_pca_object"]]<-res.pca
  pca_result[["testing_pca_object"]]<-pca_pred
  pca_result[["pca_data_frame"]]<-final_data_frame
  pca_result[["pca_feature"]]<-colnames(res.pca$ind$coord)
  return(pca_result)
}

autoencoder_feature_generate<-function(feature_data_frame=feature_data$cluster_feature_dataframe,feature_list=study_feature$frag_arm,trainingunit="training",testingunit="testing",pc_number=50,seed_num){
  training_data_frame<-feature_data_frame[feature_data_frame$unit==trainingunit,feature_list]
  
  testing_data_frame<-feature_data_frame[feature_data_frame$unit==testingunit,feature_list]
  #testing_data_frame[is.na(testing_data_frame)]<-0
  
  is.na_arm<-as.data.frame(apply(training_data_frame,2,sum_is_na))
  colnames(is.na_arm)[[1]]<-"NA_count"
  keep_arm<-rownames(is.na_arm)[is.na_arm$NA_count!=nrow(training_data_frame)]
  training_data_frame<-training_data_frame[,keep_arm]
  training_data_frame[is.na(training_data_frame)]<-0
  testing_data_frame<-testing_data_frame[,keep_arm]
  testing_data_frame[is.na(testing_data_frame)]<-0
  #testing_data_frame<-testing_data_frame[1,]
  
  training_feature <- as.h2o(training_data_frame)
  testing_feature<-as.h2o(testing_data_frame)
  
  dl_model <- h2o.deeplearning(
    x = seq_along(training_feature),
    training_frame = training_feature,
    autoencoder = TRUE,
    hidden = pc_number,
    activation = 'Tanh',
    sparse = TRUE,
    reproducible=T,
    seed=seed_num
  )
  
  training_PC <- h2o.deepfeatures(dl_model, training_feature, layer = 1)
  testing_PC<- h2o.deepfeatures(dl_model, testing_feature, layer = 1)
  
  training_data.frame<-feature_data_frame[feature_data_frame$unit==trainingunit,]
  testing_data.frame<-feature_data_frame[feature_data_frame$unit==testingunit,]
  training_data.frame<-cbind(training_data.frame,as.data.frame(training_PC))
  testing_data.frame<-cbind(testing_data.frame,as.data.frame(testing_PC))
  # testing_PC_data<-testing_data.frame[,c("Sample_ID",colnames(testing_PC))]
  
  
  final_data_frame<-rbind(training_data.frame,testing_data.frame)
  pca_result<-list()
  pca_result[["training_autoencoder_object"]]<-dl_model
  #pca_result[["testing_pca_object"]]<-pca_pred
  pca_result[["autoencoder_data_frame"]]<-final_data_frame
  pca_result[["autoencoder_feature"]]<-colnames(training_PC)
  return(pca_result)
}

wgcna_cluster_feature<-function(feature_data_frame=feature_data,feature_list=study_feature$frag_arm,trainingunit="training",testingunit="testing"){
  training_data_frame<-feature_data_frame[feature_data_frame$unit==trainingunit,feature_list]
  testing_data_frame<-feature_data_frame[feature_data_frame$unit==testingunit,feature_list]
  wgcna_modules<-blockwiseModules(as.data.frame(training_data_frame),corType = "bicor", maxPOutliers = 1,power = 12,networkType = "signed",deepSplit = 2,minModuleSize =5,TOMtype = "signed")
  cluster_feature<-data.frame(feature=colnames(training_data_frame),cluster=wgcna_modules$colors)
  feature_positive<-c()
  highest_feature<-c()
  for (me in colnames(wgcna_modules$MEs)){
    me_ori<-me
    me<-gsub("ME","",me)
    me_feature<-cluster_feature$feature[cluster_feature$cluster==me]
    me_matrix<-as.data.frame(training_data_frame)[,me_feature]
    me_cor<-cor(me_matrix,wgcna_modules$MEs[,me_ori],use = 'pairwise.complete.obs')
    me_cor_positive<-rownames(me_cor)[me_cor[,1]>0 & !is.na(me_cor[,1])]
    high_cor_positive<-names(me_cor[order(me_cor[,1],decreasing=T),])[[1]]
    feature_positive<-c(feature_positive,me_cor_positive)
    highest_feature<-c(highest_feature,high_cor_positive)
  }
  
  cluster_feature<-cluster_feature[cluster_feature$feature %in% feature_positive,]
  
  feature_means<-list()
  for (cluster in unique(cluster_feature$cluster)){
    cluster_feature_list<-cluster_feature$feature[cluster_feature$cluster==cluster]
    feature_matrix<-training_data_frame[,cluster_feature_list]
    feature_mean<-apply(feature_matrix,1,mean,na.rm=T)
    feature_means[[cluster]]<-as.numeric(feature_mean)
  }
  feature_means_total_training<-as.data.frame(do.call(cbind,feature_means))
  feature_means<-list()
  for (cluster in unique(cluster_feature$cluster)){
    cluster_feature_list<-cluster_feature$feature[cluster_feature$cluster==cluster]
    feature_matrix<-testing_data_frame[,cluster_feature_list]
    feature_mean<-apply(feature_matrix,1,mean,na.rm=T)
    feature_means[[cluster]]<-as.numeric(feature_mean)
  }
  feature_means_total_testing<-as.data.frame(do.call(cbind,feature_means))
  training_data<-feature_data_frame[feature_data_frame$unit==trainingunit,]
  testing_data<-feature_data_frame[feature_data_frame$unit==testingunit,]
  training_data<-cbind(training_data,feature_means_total_training)
  testing_data<-cbind(testing_data,feature_means_total_testing)
  cluster_feature_data<-rbind(training_data,testing_data)
  cluster_grey_rm<-as.character(unique(cluster_feature$cluster[cluster_feature$cluster!="grey"]))
  cluster_all<-as.character(unique(cluster_feature$cluster))
  cluster_grey_rm_grey_feature<-c(as.character(cluster_grey_rm),as.character(cluster_feature$feature[cluster_feature$cluster=="grey"]))
  final_result<-list()
  final_result[["cluster_feature_dataframe"]]<-cluster_feature_data
  final_result[["cluster_feature"]]<-cluster_feature
  final_result[["cluster_all"]]<-cluster_all
  final_result[["cluster_grey_rm"]]<-cluster_grey_rm
  final_result[["cluster_grey_rm_grey_feature"]]<-cluster_grey_rm_grey_feature
  final_result[["highest_positive_feature"]]<-highest_feature
  return(final_result)
}

cross_validation_h2o_repeat_glm<-function(feature_table=local_msk_logCPM,feature_list=study_feature_MT$all_MT_gene,repeat_num=10,time="P2",response_field="PCR_binary",trainingunit=c("training"),testingunit=c("testing"),n.fold=5,feature_name="all_MT_gene"){
  auc_table_total<-list()
  cv_pred_prob<-list()
  testing_pred_prob<-list()
  feature_important<-list()
  input_table<-external_validation_h2o(motif_data_frame=feature_table,including_motif=feature_list,response=response_field,training_unit=trainingunit,testing_unit=testingunit,time_point=time,last_time="no")
  
  for (i in 1:repeat_num){
    
    glm_model<-h2o_glm_construction(h2o_data_list=input_table,n_seed=seeds[[i]],project_name=paste(feature_name,"glm",sep="_"),max_models=5,n_fold=n.fold)
    glm_auc_table<-glm_model$auc_table
    glm_auc_table$feature<-feature_name
    glm_auc_table$method<-"glm"
    glm_auc_table$seed<-seeds[[i]]
    glm_auc_table$repeat_num<-i
    auc_table_total[[i]]<-glm_auc_table
    
    glm_model$cv_prediction$feature<-feature_name
    glm_model$cv_prediction$method<-"glm"
    glm_model$cv_prediction$seed<-seeds[[i]]
    glm_model$cv_prediction$repeat_num<-i
    
    glm_model$test_prediction$feature<-feature_name
    glm_model$test_prediction$method<-"glm"
    glm_model$test_prediction$seed<-seeds[[i]]
    glm_model$test_prediction$repeat_num<-i
    
    glm_model$feature_contribution$feature<-feature_name
    glm_model$feature_contribution$method<-"glm"
    glm_model$feature_contribution$seed<-seeds[[i]]
    glm_model$feature_contribution$repeat_num<-i
    
    auc_table_total[[i]]<-glm_auc_table
    cv_pred_prob[[i]]<-glm_model$cv_prediction
    testing_pred_prob[[i]]<-glm_model$test_prediction
    feature_important[[i]]<-glm_model$feature_contribution
  }
  
  total_auc<-do.call(rbind,auc_table_total)
  cv_pred_prob<-do.call(rbind,cv_pred_prob)
  test_pred_prob<-do.call(rbind,testing_pred_prob)
  feature_imp<-do.call(rbind,feature_important)
  feature_imp<-feature_imp[feature_imp$scaled_importance!=0,]
  
  cv_average<-as.data.frame(aggregate(cv_pred_prob$p1, list(cv_pred_prob$Sample_ID), FUN=mean))
  colnames(cv_average)<-c("Sample_ID","average_p1")
  cv_average_table<-cv_pred_prob[,c("SampleID","response","detail_group","three_group","train_group",       "feature","Sample_ID")]
  cv_average_table<-unique(cv_average_table)
  cv_average_table<-merge(cv_average_table,cv_average,by="Sample_ID")
  cv_average_table<-cv_average_table[order(cv_average_table$average_p1,decreasing = T),]
  cv_average_table$p1_order<-1:nrow(cv_average_table)
  average_roc_object<-roc(cv_average_table[,"response"],cv_average_table[,"average_p1"])
  
  test_pred_prob<-test_pred_prob[test_pred_prob$repeat_num==1,]
  test_pred_prob<-test_pred_prob[order(test_pred_prob$p1,decreasing=T),]
  test_pred_prob$p1_order<-1:nrow(test_pred_prob)
  test_roc_object<-roc(test_pred_prob$response,test_pred_prob$p1)
  
  auc_table_cv<-total_auc[total_auc$group=="CV",]
  
  result<-list()
  average_auc<-mean(auc_table_cv$auc,na.rm=T)
  sd_auc<-sd(auc_table_cv$auc,na.rm=T)
  result[["auc_table"]]<-total_auc
  result[["auc_table_cv"]]<-auc_table_cv
  result[["cv_average_table"]]<-cv_average_table
  result[["cv_average_roc_object"]]<-average_roc_object
  result[["test_pred_prob"]]<-test_pred_prob
  result[["test_roc_object"]]<-test_roc_object
  result[["cv_probability"]]<-cv_pred_prob
  result[["feature_rank"]]<-feature_imp
  #result[["roc_result"]]<-roc_result
  result[["average_auc"]]<-average_auc
  result[["sd_auc"]]<-sd_auc
  
  return(result)
}

cross_validation_h2o_repeat_gbm<-function(feature_table=local_msk_logCPM,feature_list=study_feature_MT$all_MT_gene,repeat_num=10,time="P2",response_field="PCR_binary",trainingunit=c("training"),testingunit=c("testing"),n.fold=5,feature_name="all_MT_gene"){
  auc_table_total<-list()
  cv_pred_prob<-list()
  testing_pred_prob<-list()
  feature_important<-list()
  input_table<-external_validation_h2o(motif_data_frame=feature_table,including_motif=feature_list,response=response_field,training_unit=trainingunit,testing_unit=testingunit,time_point=time,last_time="no")
  
  for (i in 1:repeat_num){
    
    glm_model<-h2o_gbm_construction(h2o_data_list=input_table,n_seed=seeds[[i]],project_name=paste(feature_name,"glm",sep="_"),max_models=5,n_fold=n.fold)
    glm_auc_table<-glm_model$auc_table
    glm_auc_table$feature<-feature_name
    glm_auc_table$method<-"gbm"
    glm_auc_table$seed<-seeds[[i]]
    glm_auc_table$repeat_num<-i
    auc_table_total[[i]]<-glm_auc_table
    
    glm_model$cv_prediction$feature<-feature_name
    glm_model$cv_prediction$method<-"gbm"
    glm_model$cv_prediction$seed<-seeds[[i]]
    glm_model$cv_prediction$repeat_num<-i
    
    glm_model$test_prediction$feature<-feature_name
    glm_model$test_prediction$method<-"gbm"
    glm_model$test_prediction$seed<-seeds[[i]]
    glm_model$test_prediction$repeat_num<-i
    
    glm_model$feature_contribution$feature<-feature_name
    glm_model$feature_contribution$method<-"gbm"
    glm_model$feature_contribution$seed<-seeds[[i]]
    glm_model$feature_contribution$repeat_num<-i
    
    auc_table_total[[i]]<-glm_auc_table
    cv_pred_prob[[i]]<-glm_model$cv_prediction
    testing_pred_prob[[i]]<-glm_model$test_prediction
    feature_important[[i]]<-glm_model$feature_contribution
  }
  
  total_auc<-do.call(rbind,auc_table_total)
  cv_pred_prob<-do.call(rbind,cv_pred_prob)
  test_pred_prob<-do.call(rbind,testing_pred_prob)
  feature_imp<-do.call(rbind,feature_important)
  feature_imp<-feature_imp[feature_imp$scaled_importance!=0,]
  
  cv_average<-as.data.frame(aggregate(cv_pred_prob$p1, list(cv_pred_prob$Sample_ID), FUN=mean))
  colnames(cv_average)<-c("Sample_ID","average_p1")
  cv_average_table<-cv_pred_prob[,c("SampleID","response","detail_group","three_group","train_group",       "feature","Sample_ID")]
  cv_average_table<-unique(cv_average_table)
  cv_average_table<-merge(cv_average_table,cv_average,by="Sample_ID")
  cv_average_table<-cv_average_table[order(cv_average_table$average_p1,decreasing = T),]
  cv_average_table$p1_order<-1:nrow(cv_average_table)
  average_roc_object<-roc(cv_average_table[,"response"],cv_average_table[,"average_p1"])
  
  test_pred_prob<-test_pred_prob[test_pred_prob$repeat_num==1,]
  test_pred_prob<-test_pred_prob[order(test_pred_prob$p1,decreasing=T),]
  test_pred_prob$p1_order<-1:nrow(test_pred_prob)
  test_roc_object<-roc(test_pred_prob$response,test_pred_prob$p1)
  
  auc_table_cv<-total_auc[total_auc$group=="CV",]
  
  result<-list()
  average_auc<-mean(auc_table_cv$auc,na.rm=T)
  sd_auc<-sd(auc_table_cv$auc,na.rm=T)
  result[["auc_table"]]<-total_auc
  result[["auc_table_cv"]]<-auc_table_cv
  result[["cv_average_table"]]<-cv_average_table
  result[["cv_average_roc_object"]]<-average_roc_object
  result[["test_pred_prob"]]<-test_pred_prob
  result[["test_roc_object"]]<-test_roc_object
  result[["cv_probability"]]<-cv_pred_prob
  result[["feature_rank"]]<-feature_imp
  #result[["roc_result"]]<-roc_result
  result[["average_auc"]]<-average_auc
  result[["sd_auc"]]<-sd_auc
  
  return(result)
}

cross_validation_h2o_repeat_rf<-function(feature_table=local_msk_logCPM,feature_list=study_feature_MT$all_MT_gene,repeat_num=10,time="P2",response_field="PCR_binary",trainingunit=c("training"),testingunit=c("testing"),n.fold=5,feature_name="all_MT_gene"){
  auc_table_total<-list()
  cv_pred_prob<-list()
  testing_pred_prob<-list()
  feature_important<-list()
  input_table<-external_validation_h2o(motif_data_frame=feature_table,including_motif=feature_list,response=response_field,training_unit=trainingunit,testing_unit=testingunit,time_point=time,last_time="no")
  
  for (i in 1:repeat_num){
    
    glm_model<-h2o_rf_construction(h2o_data_list=input_table,n_seed=seeds[[i]],project_name=paste(feature_name,"glm",sep="_"),max_models=5,n_fold=n.fold)
    glm_auc_table<-glm_model$auc_table
    glm_auc_table$feature<-feature_name
    glm_auc_table$method<-"rf"
    glm_auc_table$seed<-seeds[[i]]
    glm_auc_table$repeat_num<-i
    auc_table_total[[i]]<-glm_auc_table
    
    glm_model$cv_prediction$feature<-feature_name
    glm_model$cv_prediction$method<-"rf"
    glm_model$cv_prediction$seed<-seeds[[i]]
    glm_model$cv_prediction$repeat_num<-i
    
    glm_model$test_prediction$feature<-feature_name
    glm_model$test_prediction$method<-"rf"
    glm_model$test_prediction$seed<-seeds[[i]]
    glm_model$test_prediction$repeat_num<-i
    
    glm_model$feature_contribution$feature<-feature_name
    glm_model$feature_contribution$method<-"rf"
    glm_model$feature_contribution$seed<-seeds[[i]]
    glm_model$feature_contribution$repeat_num<-i
    
    auc_table_total[[i]]<-glm_auc_table
    cv_pred_prob[[i]]<-glm_model$cv_prediction
    testing_pred_prob[[i]]<-glm_model$test_prediction
    feature_important[[i]]<-glm_model$feature_contribution
  }
  
  total_auc<-do.call(rbind,auc_table_total)
  cv_pred_prob<-do.call(rbind,cv_pred_prob)
  test_pred_prob<-do.call(rbind,testing_pred_prob)
  feature_imp<-do.call(rbind,feature_important)
  feature_imp<-feature_imp[feature_imp$scaled_importance!=0,]
  
  cv_average<-as.data.frame(aggregate(cv_pred_prob$p1, list(cv_pred_prob$Sample_ID), FUN=mean))
  colnames(cv_average)<-c("Sample_ID","average_p1")
  cv_average_table<-cv_pred_prob[,c("SampleID","response","detail_group","three_group","train_group",       "feature","Sample_ID")]
  cv_average_table<-unique(cv_average_table)
  cv_average_table<-merge(cv_average_table,cv_average,by="Sample_ID")
  cv_average_table<-cv_average_table[order(cv_average_table$average_p1,decreasing = T),]
  cv_average_table$p1_order<-1:nrow(cv_average_table)
  average_roc_object<-roc(cv_average_table[,"response"],cv_average_table[,"average_p1"])
  
  test_pred_prob<-test_pred_prob[test_pred_prob$repeat_num==1,]
  test_pred_prob<-test_pred_prob[order(test_pred_prob$p1,decreasing=T),]
  test_pred_prob$p1_order<-1:nrow(test_pred_prob)
  test_roc_object<-roc(test_pred_prob$response,test_pred_prob$p1)
  
  auc_table_cv<-total_auc[total_auc$group=="CV",]
  
  result<-list()
  average_auc<-mean(auc_table_cv$auc,na.rm=T)
  sd_auc<-sd(auc_table_cv$auc,na.rm=T)
  result[["auc_table"]]<-total_auc
  result[["auc_table_cv"]]<-auc_table_cv
  result[["cv_average_table"]]<-cv_average_table
  result[["cv_average_roc_object"]]<-average_roc_object
  result[["test_pred_prob"]]<-test_pred_prob
  result[["test_roc_object"]]<-test_roc_object
  result[["cv_probability"]]<-cv_pred_prob
  result[["feature_rank"]]<-feature_imp
  #result[["roc_result"]]<-roc_result
  result[["average_auc"]]<-average_auc
  result[["sd_auc"]]<-sd_auc
  
  return(result)
}

cross_validation_h2o_repeat_dl<-function(feature_table=local_msk_logCPM,feature_list=study_feature_MT$all_MT_gene,repeat_num=10,time="P2",response_field="PCR_binary",trainingunit=c("training"),testingunit=c("testing"),n.fold=5,feature_name="all_MT_gene"){
  auc_table_total<-list()
  cv_pred_prob<-list()
  testing_pred_prob<-list()
  feature_important<-list()
  input_table<-external_validation_h2o(motif_data_frame=feature_table,including_motif=feature_list,response=response_field,training_unit=trainingunit,testing_unit=testingunit,time_point=time,last_time="no")
  
  for (i in 1:repeat_num){
    
    glm_model<-h2o_dl_construction(h2o_data_list=input_table,n_seed=seeds[[i]],project_name=paste(feature_name,"glm",sep="_"),max_models=5,n_fold=n.fold)
    glm_auc_table<-glm_model$auc_table
    glm_auc_table$feature<-feature_name
    glm_auc_table$method<-"dl"
    glm_auc_table$seed<-seeds[[i]]
    glm_auc_table$repeat_num<-i
    auc_table_total[[i]]<-glm_auc_table
    
    glm_model$cv_prediction$feature<-feature_name
    glm_model$cv_prediction$method<-"dl"
    glm_model$cv_prediction$seed<-seeds[[i]]
    glm_model$cv_prediction$repeat_num<-i
    
    glm_model$test_prediction$feature<-feature_name
    glm_model$test_prediction$method<-"dl"
    glm_model$test_prediction$seed<-seeds[[i]]
    glm_model$test_prediction$repeat_num<-i
    
    glm_model$feature_contribution$feature<-feature_name
    glm_model$feature_contribution$method<-"rf"
    glm_model$feature_contribution$seed<-seeds[[i]]
    glm_model$feature_contribution$repeat_num<-i
    
    auc_table_total[[i]]<-glm_auc_table
    cv_pred_prob[[i]]<-glm_model$cv_prediction
    testing_pred_prob[[i]]<-glm_model$test_prediction
    feature_important[[i]]<-glm_model$feature_contribution
  }
  
  total_auc<-do.call(rbind,auc_table_total)
  cv_pred_prob<-do.call(rbind,cv_pred_prob)
  test_pred_prob<-do.call(rbind,testing_pred_prob)
  feature_imp<-do.call(rbind,feature_important)
  feature_imp<-feature_imp[feature_imp$scaled_importance!=0,]
  
  cv_average<-as.data.frame(aggregate(cv_pred_prob$p1, list(cv_pred_prob$Sample_ID), FUN=mean))
  colnames(cv_average)<-c("Sample_ID","average_p1")
  cv_average_table<-cv_pred_prob[,c("SampleID","response","detail_group","three_group","train_group",       "feature","Sample_ID")]
  cv_average_table<-unique(cv_average_table)
  cv_average_table<-merge(cv_average_table,cv_average,by="Sample_ID")
  cv_average_table<-cv_average_table[order(cv_average_table$average_p1,decreasing = T),]
  cv_average_table$p1_order<-1:nrow(cv_average_table)
  average_roc_object<-roc(cv_average_table[,"response"],cv_average_table[,"average_p1"])
  
  test_pred_prob<-test_pred_prob[test_pred_prob$repeat_num==1,]
  test_pred_prob<-test_pred_prob[order(test_pred_prob$p1,decreasing=T),]
  test_pred_prob$p1_order<-1:nrow(test_pred_prob)
  test_roc_object<-roc(test_pred_prob$response,test_pred_prob$p1)
  
  auc_table_cv<-total_auc[total_auc$group=="CV",]
  
  result<-list()
  average_auc<-mean(auc_table_cv$auc,na.rm=T)
  sd_auc<-sd(auc_table_cv$auc,na.rm=T)
  result[["auc_table"]]<-total_auc
  result[["auc_table_cv"]]<-auc_table_cv
  result[["cv_average_table"]]<-cv_average_table
  result[["cv_average_roc_object"]]<-average_roc_object
  result[["test_pred_prob"]]<-test_pred_prob
  result[["test_roc_object"]]<-test_roc_object
  result[["cv_probability"]]<-cv_pred_prob
  result[["feature_rank"]]<-feature_imp
  #result[["roc_result"]]<-roc_result
  result[["average_auc"]]<-average_auc
  result[["sd_auc"]]<-sd_auc
  
  return(result)
}

cohort_feature_validation<-function(cohort=local_msk_logCPM,including_feature,pc_num=50,repeat_number=1,feature.name){
  
  cohort_wgcna_pca<-PCA_feature_generate(feature_data_frame=cohort,feature_list=including_feature,pc_number = 50)
  cohort_wgcna_pca_autoencoder<-autoencoder_feature_generate(feature_data_frame=cohort_wgcna_pca$pca_data_frame,feature_list=including_feature,seed_num=1730,pc_number=50)
  feature_glm_model<-cross_validation_h2o_repeat_glm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=including_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  feature_gbm_model<-cross_validation_h2o_repeat_gbm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=including_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  feature_dl_model<-cross_validation_h2o_repeat_dl(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=including_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  feature_rf_model<-cross_validation_h2o_repeat_rf(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=including_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  cluster_glm_model<-cross_validation_h2o_repeat_glm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna$cluster_all, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  cluster_gbm_model<-cross_validation_h2o_repeat_gbm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna$cluster_all, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  cluster_dl_model<-cross_validation_h2o_repeat_dl(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna$cluster_all, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  cluster_rf_model<-cross_validation_h2o_repeat_rf(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna$cluster_all, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  pca_glm_model<-cross_validation_h2o_repeat_glm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca$pca_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  pca_gbm_model<-cross_validation_h2o_repeat_gbm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca$pca_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  pca_dl_model<-cross_validation_h2o_repeat_dl(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca$pca_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  pca_rf_model<-cross_validation_h2o_repeat_rf(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca$pca_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  autoencoder_glm_model<-cross_validation_h2o_repeat_glm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca_autoencoder$autoencoder_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  autoencoder_gbm_model<-cross_validation_h2o_repeat_gbm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca_autoencoder$autoencoder_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  autoencoder_dl_model<-cross_validation_h2o_repeat_dl(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca_autoencoder$autoencoder_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  autoencoder_rf_model<-cross_validation_h2o_repeat_rf(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca_autoencoder$autoencoder_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  result<-list()
  result[["feature_glm_model"]]<-feature_glm_model
  result[["feature_gbm_model"]]<-feature_gbm_model
  result[["feature_rf_model"]]<-feature_rf_model
  result[["feature_dl_model"]]<-feature_dl_model
  result[["cluster_glm_model"]]<-cluster_glm_model
  result[["cluster_gbm_model"]]<-cluster_gbm_model
  result[["cluster_rf_model"]]<-cluster_rf_model
  result[["cluster_dl_model"]]<-cluster_dl_model
  result[["pca_glm_model"]]<-pca_glm_model
  result[["pca_gbm_model"]]<-pca_gbm_model
  result[["pca_rf_model"]]<-pca_rf_model
  result[["pca_dl_model"]]<-pca_dl_model
  result[["autoencoder_glm_model"]]<-autoencoder_glm_model
  result[["autoencoder_gbm_model"]]<-autoencoder_gbm_model
  result[["autoencoder_rf_model"]]<-autoencoder_rf_model
  result[["autoencoder_dl_model"]]<-autoencoder_dl_model
  return(result)
  
}


cohort_feature_validation_l2<-function(cohort=local_msk_logCPM,including_feature,pc_num=50,repeat_number=1,feature.name){
  cohort_wgcna<-wgcna_cluster_feature(feature_data_frame=cohort,feature_list=including_feature)
  cohort_wgcna_pca<-PCA_feature_generate(feature_data_frame=cohort,feature_list=including_feature,pc_number = 50)
  cohort_wgcna_pca_autoencoder<-autoencoder_feature_generate(feature_data_frame=cohort_wgcna_pca$pca_data_frame,feature_list=including_feature,seed_num=1730,pc_number=50)
  feature_glm_model<-cross_validation_h2o_repeat_glm(feature_table=cohort,feature_list=including_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  feature_gbm_model<-cross_validation_h2o_repeat_gbm(feature_table=cohort,feature_list=including_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  feature_dl_model<-cross_validation_h2o_repeat_dl(feature_table=cohort,feature_list=including_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  feature_rf_model<-cross_validation_h2o_repeat_rf(feature_table=cohort,feature_list=including_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  cluster_glm_model<-cross_validation_h2o_repeat_glm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna$cluster_all, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  cluster_gbm_model<-cross_validation_h2o_repeat_gbm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna$cluster_all, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  cluster_dl_model<-cross_validation_h2o_repeat_dl(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna$cluster_all, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  cluster_rf_model<-cross_validation_h2o_repeat_rf(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna$cluster_all, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  pca_glm_model<-cross_validation_h2o_repeat_glm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca$pca_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  pca_gbm_model<-cross_validation_h2o_repeat_gbm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca$pca_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  pca_dl_model<-cross_validation_h2o_repeat_dl(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca$pca_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  pca_rf_model<-cross_validation_h2o_repeat_rf(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca$pca_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  autoencoder_glm_model<-cross_validation_h2o_repeat_glm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca_autoencoder$autoencoder_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  autoencoder_gbm_model<-cross_validation_h2o_repeat_gbm(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca_autoencoder$autoencoder_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  autoencoder_dl_model<-cross_validation_h2o_repeat_dl(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca_autoencoder$autoencoder_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  autoencoder_rf_model<-cross_validation_h2o_repeat_rf(feature_table=cohort_wgcna_pca_autoencoder$autoencoder_data_frame,feature_list=cohort_wgcna_pca_autoencoder$autoencoder_feature, repeat_num=repeat_number,feature_name=feature.name,trainingunit = "training",testingunit="testing")
  result<-list()
  result[["feature_glm_model"]]<-feature_glm_model
  result[["feature_gbm_model"]]<-feature_gbm_model
  result[["feature_rf_model"]]<-feature_rf_model
  result[["feature_dl_model"]]<-feature_dl_model
  result[["cluster_glm_model"]]<-cluster_glm_model
  result[["cluster_gbm_model"]]<-cluster_gbm_model
  result[["cluster_rf_model"]]<-cluster_rf_model
  result[["cluster_dl_model"]]<-cluster_dl_model
  result[["pca_glm_model"]]<-pca_glm_model
  result[["pca_gbm_model"]]<-pca_gbm_model
  result[["pca_rf_model"]]<-pca_rf_model
  result[["pca_dl_model"]]<-pca_dl_model
  result[["autoencoder_glm_model"]]<-autoencoder_glm_model
  result[["autoencoder_gbm_model"]]<-autoencoder_gbm_model
  result[["autoencoder_rf_model"]]<-autoencoder_rf_model
  result[["autoencoder_dl_model"]]<-autoencoder_dl_model
  return(result)
  
}

training_test_dataframe_l2<-function(sample_data=sample_info,feature_data=feature_table,response_field="MSI",positive_response=c(1),negative_response=c(0),training_sample=kz38_sample_info$SampleID,testing_sample=kz38_sample_info$SampleID,feature_list=study_feature_kz37$griffin,featurename="griffin"){
  #real_feature<-intersect(intersect(colnames(trainingset),colnames(testingset)),feature_list)
  sample_data_frame<-sample_data[sample_data$SampleID %in% c(training_sample,testing_sample),]
  sample_data_frame<-sample_data_frame[,c("Sample_ID","SampleID","Patient_ID",response_field,"Time_point","Detail_Group")]
  sample_data_frame<-merge(sample_data_frame,feature_data,by="SampleID")
  sample_data_frame<-sample_data_frame[!is.na(sample_data_frame[[response_field]]),]
  sample_data_frame$PCR_binary<-NA
  sample_data_frame$PCR_binary[sample_data_frame[,response_field] %in% positive_response]<-1
  sample_data_frame$PCR_binary[sample_data_frame[,response_field] %in% negative_response]<-0
  sample_data_frame$PCR_binary<-as.factor(sample_data_frame$PCR_binary)
  training_data<-sample_data_frame[sample_data_frame$SampleID %in% training_sample,]
  training_data$unit<-"training"
  testing_data<-sample_data_frame[sample_data_frame$SampleID %in% testing_sample,]
  testing_data$unit<-"testing"
  sample_data_frame<-rbind(training_data,testing_data)
  
  
  result<-cohort_feature_validation_l2(cohort=sample_data_frame,including_feature=feature_list,feature.name = featurename)
  return(result)
}
