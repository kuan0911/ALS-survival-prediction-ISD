source('Data/ALSdata-preprocessing-v2.R')
source('Data/clinical_variables.R')
library(reticulate)
use_condaenv('base')
source('analysisMaster.R')
source('Evaluations/aveMetrics.R')



go = function() {
  results_dir = 'results_temp/'
  if(!dir.exists(results_dir)){dir.create(results_dir)}
  
  survivalDataset = ALS_preprocessing_v2(PSCID = F,reduce_feature = F,image_feature=T,only_image_feature=F)
  clinical_vars = clinical_vars_func()
  
  survivalDataset = validateAndClean(survivalDataset, imputeZero=F)
  
  sum(survivalDataset$delta)
  
  variable = colnames(survivalDataset)
  variable_mean = rep(NA,length(variable))
  variable_sd = rep(NA,length(variable))
  for(k in 1:length(variable)) {
    if(is.numeric(survivalDataset[,k])) {
      variable_mean[k] = round(mean(survivalDataset[,k],na.rm=T),2)
      variable_sd[k] = round(sd(survivalDataset[,k],na.rm=T),2)
    }
  }
  variable_stats = cbind(variable,variable_mean,variable_sd)
  write.csv(variable_stats,paste0(results_dir,"variable_stats_all.csv"),row.names=F)
  
  
  
  numberOfFolds = 5
  modelname_vec = c('AFT','CoxEN','MTLR','RSF')
  FS_vec = c('none','image-FS','clinical-FS','seperate','together','image-only','clinical-only','image-FS-only','clinical-FS-only')
  loss_mat = matrix(-1,numberOfFolds,length(modelname_vec)*length(FS_vec))
  loss_mat_filename = paste0(results_dir,'loss_mat.csv')
  hyperparameter_setting_filename = paste0(results_dir,'hyperparameter_setting.csv')
  loss_mat_colnames = rep('',length(modelname_vec)*length(FS_vec))
  ISD_all = data.frame()
  
  if(!file.exists(hyperparameter_setting_filename)) {
    setting_cv = data.frame(modelname=character(),FS=character(),valid=logical(),stringsAsFactors=F)
    write.csv(setting_cv,hyperparameter_setting_filename,row.names=F)
  }
  if(!file.exists(loss_mat_filename)) {
    write.csv(loss_mat,loss_mat_filename,row.names=F)
  }

  
  foldsAndNormalizedData = createFoldsAndNormalize(survivalDataset, numberOfFolds=numberOfFolds)
  originalIndexing = foldsAndNormalizedData[[1]]
  normalizedData = foldsAndNormalizedData[[2]]
  for(i in 1:numberOfFolds) {
    index=1
    training = normalizedData[[1]][[i]]
    testing = normalizedData[[2]][[i]]
    
    setting_cv = read.csv(hyperparameter_setting_filename)
    if(!is.na(setting_cv[i,'valid'])) {
      best_modelname = setting_cv[i,'modelname']
      best_FS = setting_cv[i,'FS']
      print(paste0('fold ',i,' use previous cv setting'))
    }else {
      for(modelname in modelname_vec) {
        for(FS in FS_vec) {
          loss_mat = read.csv(loss_mat_filename)
          if(loss_mat[i,index]==-1) {
            ISD = analysisMaster(training, modelname=modelname,additionalFS=FS, FS = T, numberOfFolds = 5,useAllData=F,clinical_vars=clinical_vars,foldNum=i,results_dir=results_dir)
            res = aveMetrics(ISD,std=T)
            loss_mat[i,index] = round(res$L1Loss,3)
            write.csv(loss_mat,loss_mat_filename,row.names=F)
          }
          index = index +1
        }
      }
      loss_vec = loss_mat[i,]
      best_index = which.min(loss_vec)-1
      best_modelname = modelname_vec[best_index/length(FS_vec)+1]
      best_FS = FS_vec[best_index%%length(FS_vec)+1]
      setting_cv[i,] = c(best_modelname,best_FS,T)
      write.csv(setting_cv,hyperparameter_setting_filename,row.names=F)
    }
    ISD_all_temp = analysisMaster(use_splited_data=T,training_input = training,testing_input = testing, modelname=setting_cv[i,'modelname'],additionalFS=setting_cv[i,'FS'], FS = T, numberOfFolds = 5,useAllData=F,clinical_vars=clinical_vars,foldNum=i,results_dir=results_dir)
    
    ISD_all = rbind(ISD_all,ISD_all_temp)
  }
  
  loss_mat_final = matrix(-1,1,length(modelname_vec)*length(FS_vec))
  loss_mat_final_filename = paste0(results_dir,'loss_mat_final.csv')
  if(!file.exists(loss_mat_final_filename)) {
    write.csv(loss_mat_final,loss_mat_final_filename,row.names=F)
  }
  index=1
  for(modelname in modelname_vec) {
    for(FS in FS_vec) {
      loss_mat_final = read.csv(loss_mat_final_filename)
      if(loss_mat_final[1,index]==-1) {
        ISD_final_cv = analysisMaster(survivalDataset, modelname=modelname,additionalFS=FS, FS = T, numberOfFolds = 5,useAllData=F,clinical_vars=clinical_vars,foldNum='all_cv',results_dir=results_dir)
        res_final_cv = aveMetrics(ISD_final_cv,std=T)
        loss_mat_final[1,index] = round(res_final_cv$L1Loss,3)
        write.csv(loss_mat_final,loss_mat_final_filename,row.names=F)
      }
      index = index +1
    }
  }
  best_index_final = which.min(loss_mat_final[1,])-1
  best_modelname = modelname_vec[best_index_final/length(FS_vec)+1]
  best_FS = FS_vec[best_index_final%%length(FS_vec)+1]
  print(paste0('Final CV. Best model: ',best_modelname,'. Best FS: ',best_FS))
  print(ISD_all)
  print(paste0('L1 loss final: ',mean(ISD_all$L1Loss),' C-index: ',mean(ISD_all$Concordance)))
  
  return(ISD_all)
}

plot_figures = function() {
  results_dir = 'results_temp/'
  survivalDataset = ALS_preprocessing_v2(PSCID = F,reduce_feature = F,image_feature=T,only_image_feature=F)
  clinical_vars = clinical_vars_func()
  
  survivalDataset = validateAndClean(survivalDataset, imputeZero=F)
  
  hist(survivalDataset[survivalDataset$delta==1,]$time/30,ylim=c(0,40),xlim = c(0,max(survivalDataset$time/30)),breaks=10,main = "Survival Time Distribution",xlab="Time (month)",ylab="Population")
  hist(survivalDataset[survivalDataset$delta==0,]$time/30,ylim=c(0,40),xlim = c(0,max(survivalDataset$time/30)),breaks=10,main = "Censor Time Distribution",xlab="Time (month)",ylab="Population")
  hist(survivalDataset$Age,main = "Age",ylim=c(0,40),xlab="Years",ylab="Population")
  
  
  loss_mat_final_filename = paste0(results_dir,'loss_mat_final.csv')
  loss_mat_final = read.csv(loss_mat_final_filename)
  modelname_vec = c('AFT','CoxEN','MTLR','RSF')
  FS_vec = c('none','image-FS','clinical-FS','seperate','together','image-only','clinical-only','image-FS-only','clinical-FS-only')
  best_index_final = which.min(loss_mat_final[1,])-1
  best_modelname = modelname_vec[best_index_final/length(FS_vec)+1]
  best_FS = FS_vec[best_index_final%%length(FS_vec)+1]
  print(paste0('Final CV. Best model: ',best_modelname,'. Best FS: ',best_FS))
  
  ISD_plot = analysisMaster(survivalDataset, modelname=best_modelname,additionalFS=best_FS, FS = T, numberOfFolds = 5,useAllData=T,clinical_vars=clinical_vars,foldNum='plot',results_dir=results_dir)
  ISD_plot$survivalCurves[[best_modelname]]$time = ISD_plot$survivalCurves[[best_modelname]]$time/30
  
  
  plot_survivalCurves = ISD_plot$survivalCurves[[best_modelname]]
  plot_time_index = floor(length(plot_survivalCurves$time)/2)
  plot_median_time = plot_survivalCurves$time[plot_time_index]
  plot_survival_prob_vec = unlist(plot_survivalCurves[plot_time_index,2:ncol(plot_survivalCurves)])
  plot_index_min = which.min(plot_survival_prob_vec+(1-survivalDataset$delta)*1e5)
  plot_index_max = which.max(plot_survival_prob_vec-(1-survivalDataset$delta)*1e5)
  plot_censor_median = sort(survivalDataset[survivalDataset$delta==0,]$time)[floor(nrow(survivalDataset[survivalDataset$delta==0,])*9/16)]
  plot_index_censor = match(plot_censor_median,survivalDataset$time)
  
  
  plot_point_time_min = round(survivalDataset$time[plot_index_min]/30,1)
  plot_point_prob_min = round(predictProbabilityFromCurve(plot_survivalCurves[,plot_index_min+1],plot_survivalCurves$time,plot_point_time_min),2)
  plot_point_time_max = round(survivalDataset$time[plot_index_max]/30,1)
  plot_point_prob_max = round(predictProbabilityFromCurve(plot_survivalCurves[,plot_index_max+1],plot_survivalCurves$time,plot_point_time_max),2)
  plot_point_time_censor = round(survivalDataset$time[plot_index_censor]/30,1)
  plot_point_prob_censor = round(predictProbabilityFromCurve(plot_survivalCurves[,plot_index_censor+1],plot_survivalCurves$time,plot_point_time_censor),2)
  
  plot_point_median_time_min = round(predictMedianSurvivalTimeSpline(plot_survivalCurves[,plot_index_min+1],plot_survivalCurves$time),1)
  plot_point_median_time_max = round(predictMedianSurvivalTimeSpline(plot_survivalCurves[,plot_index_max+1],plot_survivalCurves$time),1)
  plot_point_median_time_censor = round(predictMedianSurvivalTimeSpline(plot_survivalCurves[,plot_index_censor+1],plot_survivalCurves$time),1)
  
  plot_spline = splinefun(plot_survivalCurves[,plot_index_max+1],plot_survivalCurves$time, method = "hyman")
  plot_point_time_ref = round(plot_spline(0.9),1)
  
  plot_points = data.frame(x = c(plot_point_time_ref,plot_point_median_time_min,plot_point_median_time_max,plot_point_median_time_censor,plot_point_time_min,plot_point_time_max,plot_point_time_censor),y=c(0.9,0.5,0.5,0.5,plot_point_prob_min,plot_point_prob_max,plot_point_prob_censor),filled=c(F,F,F,F,T,T,F),colour=c(rep('darkgreen',4),rep('black',3)))
  
  
  kmMod = prodlim(Surv(time/30,delta)~1, data = survivalDataset)
  KM = predict(kmMod,plot_survivalCurves$time)
  plot_survivalCurves  = cbind(plot_survivalCurves,KM)
  plot_km_index = ncol(plot_survivalCurves)-1
  KM_median_time = round(predictMedianSurvivalTimeSpline(KM,plot_survivalCurves$time),1)
  curve_color = c('grey','green','orange','blue')
  
  
  plotSurvivalCurves(plot_survivalCurves, c(plot_km_index,plot_index_min,plot_index_max,plot_index_censor), title='Individual Survival Distribution',plot_points=plot_points,color = curve_color)
  
  print(paste0('Plot: min: (',plot_point_time_min,', ',plot_point_prob_min,'). max: (',plot_point_time_max,', ',plot_point_prob_max,'). censor: (',plot_point_time_censor,', ',plot_point_prob_censor,')'))
}

misc = function() {
  results_dir = 'results_temp/'
  survivalDataset = ALS_preprocessing_v2(PSCID = F,reduce_feature = F,image_feature=T,only_image_feature=F)
  clinical_vars = clinical_vars_func()
  
  survivalDataset = validateAndClean(survivalDataset, imputeZero=F)
  
  loss_mat_final_filename = paste0(results_dir,'loss_mat_final.csv')
  loss_mat_final = read.csv(loss_mat_final_filename)
  modelname_vec = c('AFT','CoxEN','MTLR','RSF')
  FS_vec = c('none','image-FS','clinical-FS','seperate','together','image-only','clinical-only','image-FS-only','clinical-FS-only')
  best_index_final = which.min(loss_mat_final[1,])-1
  best_modelname = modelname_vec[best_index_final/length(FS_vec)+1]
  best_FS = FS_vec[best_index_final%%length(FS_vec)+1]
  print(paste0('Final CV. Best model: ',best_modelname,'. Best FS: ',best_FS))
  
  ISD_plot = analysisMaster(survivalDataset, modelname=best_modelname,additionalFS=best_FS, FS = T, numberOfFolds = 5,useAllData=T,clinical_vars=clinical_vars,foldNum='plot',results_dir=results_dir)
  ISD_plot$survivalCurves[[best_modelname]]$time = ISD_plot$survivalCurves[[best_modelname]]$time/30
  
  
  plot_censor_median = sort(survivalDataset[survivalDataset$delta==0,]$time)[floor(nrow(survivalDataset[survivalDataset$delta==0,])*9/16)]
  plot_index_censor = match(plot_censor_median,survivalDataset$time)
  
  plot_survivalCurves = ISD_plot$survivalCurves[[best_modelname]]
  plot_point_time_censor = round(survivalDataset$time[plot_index_censor]/30,1)
  plot_point_prob_censor = round(predictProbabilityFromCurve(plot_survivalCurves[,plot_index_censor+1],plot_survivalCurves$time,plot_point_time_censor),2)
  
  kmMod = prodlim(Surv(time/30,delta)~1, data = survivalDataset)
  KM = predict(kmMod,plot_survivalCurves$time)
  
  # calculate MAE-margin death time estimation for censored data
  MAE_margin_prob = KM/plot_point_prob_censor
  KM_median_time = round(predictMedianSurvivalTimeSpline(MAE_margin_prob,plot_survivalCurves$time),1)
  
  # baseline kaplan-meier estimator result
  ISD_KM = analysisMaster(survivalDataset, modelname='KM',additionalFS='none', FS = T, numberOfFolds = 5,useAllData=F,clinical_vars=clinical_vars,foldNum='all_cv',results_dir=results_dir)
  res_KM = aveMetrics(ISD_KM,std=T)
}

compare_img = function() {
  results_dir = 'results_temp/'
  if(!dir.exists(results_dir)){dir.create(results_dir)}
  
  survivalDataset = ALS_preprocessing_v2(PSCID = F,reduce_feature = F,image_feature=T,only_image_feature=F)
  clinical_vars = clinical_vars_func()
  survivalDataset = validateAndClean(survivalDataset, imputeZero=F)
  
  numberOfFolds = 5
  modelname_vec = c('AFT','CoxEN','MTLR','RSF')
  FS_vec_clinical = c('clinical-only','clinical-FS-only')
  FS_vec_image = c('image-only','image-FS-only')
  loss_mat_clinical = matrix(-1,numberOfFolds,length(modelname_vec)*length(FS_vec_clinical))
  loss_mat_clinical_filename = paste0(results_dir,'loss_mat_clinical.csv')
  loss_mat_image = matrix(-1,numberOfFolds,length(modelname_vec)*length(FS_vec_image))
  loss_mat_image_filename = paste0(results_dir,'loss_mat_image.csv')
  hyperparameter_setting_clinical_filename = paste0(results_dir,'hyperparameter_setting_clinical.csv')
  hyperparameter_setting_image_filename = paste0(results_dir,'hyperparameter_setting_image.csv')
  
  ISD_all_clinical = data.frame()
  ISD_all_image = data.frame()
  
  if(!file.exists(hyperparameter_setting_clinical_filename)) {
    setting_cv = data.frame(modelname=character(),FS=character(),valid=logical(),stringsAsFactors=F)
    write.csv(setting_cv,hyperparameter_setting_clinical_filename,row.names=F)
  }
  if(!file.exists(loss_mat_clinical_filename)) {
    write.csv(loss_mat_clinical,loss_mat_clinical_filename,row.names=F)
  }
  if(!file.exists(hyperparameter_setting_image_filename)) {
    setting_cv = data.frame(modelname=character(),FS=character(),valid=logical(),stringsAsFactors=F)
    write.csv(setting_cv,hyperparameter_setting_image_filename,row.names=F)
  }
  if(!file.exists(loss_mat_image_filename)) {
    write.csv(loss_mat_image,loss_mat_image_filename,row.names=F)
  }
  
  
  foldsAndNormalizedData = createFoldsAndNormalize(survivalDataset, numberOfFolds=numberOfFolds)
  originalIndexing = foldsAndNormalizedData[[1]]
  normalizedData = foldsAndNormalizedData[[2]]
  for(i in 1:numberOfFolds) {
    training = normalizedData[[1]][[i]]
    testing = normalizedData[[2]][[i]]

    for(var_type in c('clinical','image')) {
      if(var_type=='clinical') {
        FS_vec = FS_vec_clinical
        loss_mat_filename = loss_mat_clinical_filename
        hyperparameter_setting_filename = hyperparameter_setting_clinical_filename
      }else if(var_type=='image') {
        FS_vec = FS_vec_image
        loss_mat_filename = loss_mat_image_filename
        hyperparameter_setting_filename = hyperparameter_setting_image_filename
      }
      
      setting_cv = read.csv(hyperparameter_setting_filename)
      if(!is.na(setting_cv[i,'valid'])) {
        best_modelname = setting_cv[i,'modelname']
        best_FS = setting_cv[i,'FS']
        print(paste0('fold ',i,' use previous cv setting'))
      }else {
        index = 1
        for(modelname in modelname_vec) {
          for(FS in FS_vec) {
            loss_mat = read.csv(loss_mat_filename)
            if(loss_mat[i,index]==-1) {
              ISD = analysisMaster(training, modelname=modelname,additionalFS=FS, FS = T, numberOfFolds = 5,useAllData=F,clinical_vars=clinical_vars,foldNum=i,results_dir=results_dir)
              res = aveMetrics(ISD,std=T,show=T)
              loss_mat[i,index] = round(res$L1Loss,3)
              write.csv(loss_mat,loss_mat_filename,row.names=F)
            }
            index = index +1
          }
        }
        loss_vec = loss_mat[i,]
        best_index = which.min(loss_vec)-1
        best_modelname = modelname_vec[best_index/length(FS_vec)+1]
        best_FS = FS_vec[best_index%%length(FS_vec)+1]
        setting_cv[i,] = c(best_modelname,best_FS,T)
        write.csv(setting_cv,hyperparameter_setting_filename,row.names=F)
      }
      ISD_all_temp = analysisMaster(use_splited_data=T,training_input = training,testing_input = testing, modelname=setting_cv[i,'modelname'],additionalFS=setting_cv[i,'FS'], FS = T, numberOfFolds = 5,useAllData=F,clinical_vars=clinical_vars,foldNum=i,results_dir=results_dir)
      
      if(var_type == 'clinical') {ISD_all_clinical = rbind(ISD_all_clinical,ISD_all_temp)}
      else if(var_type == 'image') {ISD_all_image = rbind(ISD_all_image,ISD_all_temp)}
    }
  }
  
  print('clinical')
  print(ISD_all_clinical)
  print('image')
  print(ISD_all_image)
  
  print(paste0('L1 loss: clinical ',mean(ISD_all_clinical$L1Loss),' image ',mean(ISD_all_image$L1Loss)))
  print(paste0('C-index: clinical ',mean(ISD_all_clinical$Concordance),' image ',mean(ISD_all_image$Concordance)))
  return(list(ISD_all_clinical = ISD_all_clinical, ISD_all_image = ISD_all_image))
}


variable_stats = function() {
  results_dir = 'results_temp/'
  if(!dir.exists(results_dir)){dir.create(results_dir)}
  
  survivalDataset = ALS_preprocessing_v2(PSCID = F,reduce_feature = F,image_feature=T,only_image_feature=F)
  clinical_vars = clinical_vars_func()
  
  survivalDataset = validateAndClean(survivalDataset, imputeZero=F)
  
  sum(survivalDataset$delta)
  
  variable = colnames(survivalDataset)
  variable_mean = rep(NA,length(variable))
  variable_sd = rep(NA,length(variable))
  for(k in 1:length(variable)) {
    if(is.numeric(survivalDataset[,k])) {
      variable_mean[k] = round(mean(survivalDataset[,k],na.rm=T),2)
      variable_sd[k] = round(sd(survivalDataset[,k],na.rm=T),2)
    }
  }
  variable_stats = cbind(variable,variable_mean,variable_sd)
  write.csv(variable_stats,paste0(results_dir,"variable_stats_all.csv"),row.names=F)
}


