load("feature_table.RData")
library(caret)
library(splitTools)
library(e1071)
library(glmnet)
library(pROC)
library(ggplot2)

#total_local_cohort including all of 103 patients
total_local_cohort

#training_local_cohort including 68 patients for model training
training_local_cohort

#testing_local_cohort including 35 patient for model testing
testing_local_cohort

#ct_DNA_undetection_cohort including 74 patients with ctDNA undetection
ct_DNA_undetection_cohort

#study_feature is a list containing 22 kinds of feature or feature combinations
study_feature


seeds<-sample(1:5000,30, replace=F)

# used for cross validation repeating 30 times
cross_validation_repeat<-function(motif_table=training_local_cohort,include_feature=study_feature$T234_motif_6mer,repeat_number=30){
  
  model_variable<-list()
  best_parameter<-list()
  all_prediction_repeat<-list()
  metric_repeat_youden<-list()
  metric_repeat_all_threshold<-list()
  metric_repeat_accuracy<-list()
  
  for (j in 1:repeat_number){
    folds<-create_folds(motif_table$PCR_status, k = 5,type = c("stratified"),seed=seeds[[j]])
    all_prediction<-list()
    for (i in 1:5){
      #train.data<-motif_table[40:103,]
      #test.data<-motif_table[1:39,]
      #train.data<-ky087_training_data
      #test.data<-ky087_testing_data
      
      train.data  <- motif_table[folds[[i]], ]
      test.data <- motif_table[-folds[[i]], ]
      test_patient_id<-test.data$Patient.ID
      
      
      train.data<-train.data[,c(include_feature,"PCR_status")]
      test.data<-test.data[,c(include_feature,"PCR_status")]
      train.data$PCR_status<-as.factor(train.data$PCR_status)
      test.data$PCR_status<-as.factor(test.data$PCR_status)
      train.data$PCR_status<-paste("X",train.data$PCR_status,sep="")
      train.data$PCR_status<-as.factor(train.data$PCR_status)
      test.data$PCR_status<-paste("X",test.data$PCR_status,sep="")
      test.data$PCR_status<-as.factor(test.data$PCR_status)
      model <- train(
        PCR_status ~., data = train.data, method = "glmnet",
        trControl = trainControl("cv", number = 5,classProbs = TRUE,summaryFunction = twoClassSummary),family = "binomial",metric = "ROC",
        tuneLength = 10
      )
      best_alpha_lambda<-as.data.frame(model$bestTune)
      best_alpha_lambda$model<-paste(j,i,sep="_")
      coeff<- as.data.frame(as.matrix(coef(model$finalModel, model$bestTune$lambda)))
      coeff$model<-paste(j,i,sep="_")
      coeff$feature<-rownames(coeff)
      model_variable[[paste(j,i,sep="_")]]<-coeff
      
      best_parameter[[paste(j,i,sep="_")]]<-best_alpha_lambda
      
      x.test <- model.matrix(PCR_status ~., test.data)[,-1]
      predictions <- predict(model, x.test,type = "prob")
      prediction_table<-as.data.frame(predictions)
      prediction_table$PCR_status<-test.data$PCR_status
      prediction_table$Patient_ID<-test_patient_id
      prediction_table$model<-paste(j,i,sep="_")
      all_prediction[[i]]<-prediction_table
    }
    
    all_prediction_combine<-do.call(rbind,all_prediction)
    all_prediction_combine$repeat_num<-j
    all_prediction_combine$seed<-seeds[[j]]
    all_prediction_repeat[[j]]<-all_prediction_combine
    test_roc<-roc(all_prediction_combine$PCR_status,all_prediction_combine$X1) 
    rets <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv",
              "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")
    metric<-as.data.frame(coords(test_roc, "best", ret=rets, transpose = FALSE))
    metric$auc<-test_roc$auc
    metric$cutoff<-"Youden_index"
    metric$repeat_num<-j
    metric_repeat_youden[[j]]<-metric
    all_threshold<-as.data.frame(coords(test_roc, x="all","threshold", ret=rets, transpose = FALSE))
    all_threshold$repeat_num<-j
    all_threshold<-all_threshold[order(all_threshold$accuracy,decreasing=T),]
    metric_repeat_all_threshold[[j]]<-all_threshold
    metric_highest_accuracy<-as.data.frame(all_threshold[1,])
    metric_highest_accuracy$repeat_num<-j
    metric_highest_accuracy$cutoff<-"highest_accuracy"
    metric_repeat_accuracy[[j]]<-metric_highest_accuracy
  }
  model_variable_total<-do.call(rbind,model_variable)
  model_variable_total<-model_variable_total[model_variable_total$`1`!=0,]
  best_parameter<-do.call(rbind,best_parameter)
  metric_youden_total<-do.call(rbind,metric_repeat_youden)
  metric_accuracy_total<-do.call(rbind,metric_repeat_accuracy)
  threshold_total<-do.call(rbind,metric_repeat_all_threshold)
  all_prediction<-do.call(rbind,all_prediction_repeat)
  mean_auc<-mean(metric_youden_total$auc,na.rm=T)
  sd_auc<-sd(metric_youden_total$auc,na.rm=T)
  final_result<-list()
  final_result[["model_variable"]]<-model_variable_total
  final_result[["best_parameter"]]<-best_parameter
  final_result[["metric_youden"]]<-metric_youden_total
  final_result[["metric_accuracy"]]<-metric_accuracy_total
  final_result[["all_prediction"]]<-all_prediction
  final_result[["threshold_total"]]<-threshold_total
  final_result[["average_auc"]]<-mean_auc
  final_result[["sd_auc"]]<-sd_auc
  return(final_result)
}

#used for loocv validation
loocv_validation<-function(motif_table=training_local_cohort,include_feature=study_feature$T234_motif_6mer){
  
  model_variable<-list()
  best_parameter<-list()
  all_prediction_repeat<-list()
  metric_repeat_youden<-list()
  metric_repeat_all_threshold<-list()
  metric_repeat_accuracy<-list()
  
  for (j in 1:1){
    all_prediction<-list()
    for (i in 1:nrow(motif_table)){
      train.data  <- motif_table[-i,]
      test.data <- motif_table[i,]
      test_patient_id<-test.data$Patient.ID
      
      train.data<-train.data[,c(include_feature,"PCR_status")]
      test.data<-test.data[,c(include_feature,"PCR_status")]
      train.data$PCR_status<-as.factor(train.data$PCR_status)
      test.data$PCR_status<-as.factor(test.data$PCR_status)
      train.data$PCR_status<-paste("X",train.data$PCR_status,sep="")
      train.data$PCR_status<-as.factor(train.data$PCR_status)
      test.data$PCR_status<-paste("X",test.data$PCR_status,sep="")
      test.data$PCR_status<-as.factor(test.data$PCR_status)
      model <- train(
        PCR_status ~., data = train.data, method = "glmnet",
        trControl = trainControl("cv", number = 5,classProbs = TRUE,summaryFunction = twoClassSummary),family = "binomial",metric = "ROC",
        tuneLength = 10
      )
      best_alpha_lambda<-as.data.frame(model$bestTune)
      best_alpha_lambda$model<-paste(j,i,sep="_")
      coeff<- as.data.frame(as.matrix(coef(model$finalModel, model$bestTune$lambda)))
      coeff$model<-paste(j,i,sep="_")
      coeff$feature<-rownames(coeff)
      model_variable[[paste(j,i,sep="_")]]<-coeff
      
      best_parameter[[paste(j,i,sep="_")]]<-best_alpha_lambda
      
      #x.test <- model.matrix(PCR_status ~., test.data)[,-1]
      x.test<-test.data[,-ncol(test.data)]
      predictions <- predict(model, x.test,type = "prob")
      prediction_table<-as.data.frame(predictions)
      prediction_table$PCR_status<-test.data$PCR_status
      prediction_table$Patient_ID<-test_patient_id
      prediction_table$model<-paste(j,i,sep="_")
      all_prediction[[i]]<-prediction_table
    }
    
    all_prediction_combine<-do.call(rbind,all_prediction)
    all_prediction_combine$repeat_num<-j
    all_prediction_repeat[[j]]<-all_prediction_combine
    test_roc<-roc(all_prediction_combine$PCR_status,all_prediction_combine$X1) 
    rets <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv",
              "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")
    metric<-as.data.frame(coords(test_roc, "best", ret=rets, transpose = FALSE))
    metric$auc<-test_roc$auc
    metric$cutoff<-"Youden_index"
    metric$repeat_num<-j
    metric_repeat_youden[[j]]<-metric
    all_threshold<-as.data.frame(coords(test_roc, x="all","threshold", ret=rets, transpose = FALSE))
    all_threshold$repeat_num<-j
    all_threshold<-all_threshold[order(all_threshold$accuracy,decreasing=T),]
    metric_repeat_all_threshold[[j]]<-all_threshold
    metric_highest_accuracy<-as.data.frame(all_threshold[1,])
    metric_highest_accuracy$repeat_num<-j
    metric_highest_accuracy$cutoff<-"highest_accuracy"
    metric_repeat_accuracy[[j]]<-metric_highest_accuracy
  }
  model_variable_total<-do.call(rbind,model_variable)
  model_variable_total<-model_variable_total[model_variable_total$`1`!=0,]
  best_parameter<-do.call(rbind,best_parameter)
  metric_youden_total<-do.call(rbind,metric_repeat_youden)
  metric_accuracy_total<-do.call(rbind,metric_repeat_accuracy)
  threshold_total<-do.call(rbind,metric_repeat_all_threshold)
  all_prediction<-do.call(rbind,all_prediction_repeat)
  mean_auc<-mean(metric_youden_total$auc,na.rm=T)
  sd_auc<-sd(metric_youden_total$auc,na.rm=T)
  final_result<-list()
  final_result[["model_variable"]]<-model_variable_total
  final_result[["best_parameter"]]<-best_parameter
  final_result[["metric_youden"]]<-metric_youden_total
  final_result[["metric_accuracy"]]<-metric_accuracy_total
  final_result[["all_prediction"]]<-all_prediction
  final_result[["threshold_total"]]<-threshold_total
  final_result[["average_auc"]]<-mean_auc
  final_result[["sd_auc"]]<-sd_auc
  return(final_result)
}


#used for testing validation. 
testing_validation<-function(motif_table=training_local_cohort,include_feature=study_feature$T234_motif_6mer,testing_data=testing_local_cohort){
  
  model_variable<-list()
  best_parameter<-list()
  all_prediction_repeat<-list()
  metric_repeat_youden<-list()
  metric_repeat_all_threshold<-list()
  metric_repeat_accuracy<-list()
  
  for (j in 1:1){
    all_prediction<-list()
    for (i in 1:1){
      train.data  <- motif_table
      test.data <- testing_data
      test_patient_id<-test.data$Patient.ID
      
      train.data<-train.data[,c(include_feature,"PCR_status")]
      test.data<-test.data[,c(include_feature,"PCR_status")]
      train.data$PCR_status<-as.factor(train.data$PCR_status)
      test.data$PCR_status<-as.factor(test.data$PCR_status)
      train.data$PCR_status<-paste("X",train.data$PCR_status,sep="")
      train.data$PCR_status<-as.factor(train.data$PCR_status)
      test.data$PCR_status<-paste("X",test.data$PCR_status,sep="")
      test.data$PCR_status<-as.factor(test.data$PCR_status)
      model <- train(
        PCR_status ~., data = train.data, method = "glmnet",
        trControl = trainControl("cv", number = 5,classProbs = TRUE,summaryFunction = twoClassSummary),family = "binomial",metric = "ROC",
        tuneLength = 10
      )
      best_alpha_lambda<-as.data.frame(model$bestTune)
      best_alpha_lambda$model<-paste(j,i,sep="_")
      coeff<- as.data.frame(as.matrix(coef(model$finalModel, model$bestTune$lambda)))
      coeff$model<-paste(j,i,sep="_")
      coeff$feature<-rownames(coeff)
      model_variable[[paste(j,i,sep="_")]]<-coeff
      
      best_parameter[[paste(j,i,sep="_")]]<-best_alpha_lambda
      
      x.test <- model.matrix(PCR_status ~., test.data)[,-1]
      predictions <- predict(model, x.test,type = "prob")
      prediction_table<-as.data.frame(predictions)
      prediction_table$PCR_status<-test.data$PCR_status
      prediction_table$Patient_ID<-test_patient_id
      prediction_table$model<-paste(j,i,sep="_")
      all_prediction[[i]]<-prediction_table
    }
    
    all_prediction_combine<-do.call(rbind,all_prediction)
    all_prediction_combine$repeat_num<-j
    all_prediction_repeat[[j]]<-all_prediction_combine
    test_roc<-roc(all_prediction_combine$PCR_status,all_prediction_combine$X1) 
    rets <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv",
              "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")
    metric<-as.data.frame(coords(test_roc, "best", ret=rets, transpose = FALSE))
    metric$auc<-test_roc$auc
    metric$cutoff<-"Youden_index"
    metric$repeat_num<-j
    metric_repeat_youden[[j]]<-metric
    all_threshold<-as.data.frame(coords(test_roc, x="all","threshold", ret=rets, transpose = FALSE))
    all_threshold$repeat_num<-j
    all_threshold<-all_threshold[order(all_threshold$accuracy,decreasing=T),]
    metric_repeat_all_threshold[[j]]<-all_threshold
    metric_highest_accuracy<-as.data.frame(all_threshold[1,])
    metric_highest_accuracy$repeat_num<-j
    metric_highest_accuracy$cutoff<-"highest_accuracy"
    metric_repeat_accuracy[[j]]<-metric_highest_accuracy
  }
  model_variable_total<-do.call(rbind,model_variable)
  model_variable_total<-model_variable_total[model_variable_total$`1`!=0,]
  best_parameter<-do.call(rbind,best_parameter)
  metric_youden_total<-do.call(rbind,metric_repeat_youden)
  metric_accuracy_total<-do.call(rbind,metric_repeat_accuracy)
  threshold_total<-do.call(rbind,metric_repeat_all_threshold)
  all_prediction<-do.call(rbind,all_prediction_repeat)
  mean_auc<-mean(metric_youden_total$auc,na.rm=T)
  sd_auc<-sd(metric_youden_total$auc,na.rm=T)
  final_result<-list()
  final_result[["model_variable"]]<-model_variable_total
  final_result[["best_parameter"]]<-best_parameter
  final_result[["metric_youden"]]<-metric_youden_total
  final_result[["metric_accuracy"]]<-metric_accuracy_total
  final_result[["all_prediction"]]<-all_prediction
  final_result[["threshold_total"]]<-threshold_total
  final_result[["average_auc"]]<-mean_auc
  final_result[["sd_auc"]]<-sd_auc
  return(final_result)
}


