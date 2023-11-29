#------------------------------------------------------------
#Cumulative model with seperate network at each timepoint
#include all data in structure learning. Na if dead or censored
#
#
#------------------------------------------------------------

#library(e1071)
library(prodlim)
library(glmnet)
source("Models/fractionallyInclude.R")
source("Models/DHLRglm.R")
source("ValidateCleanCV/createFoldsAndNormalize.R")
source("Models/kerasHelper.R")

DHLR = function(training, testing, timePoints = 0,groundTruthTesting=NULL,debug = FALSE,ALSFRS=NULL){
  originalTraining = training
  originalTesting = testing
  
  kmMod = prodlim(Surv(time,delta)~1, data = originalTraining)

  print(timePoints)
  
  numTimepoint = length(timePoints)

  res = prepareDataKeras(training,timePoints)
  x=res$x
  y=res$y
  w=res$w
  oldx = x
  oldy = y
  oldw = w
  copy_y = y
  copy_w = matrix(0,nrow(x),ncol(y))

  if(!is.null(ALSFRS)) {
    for(i in 2:length(timePoints)) {
      for(k in 1:nrow(y)) {
        if(!is.na(ALSFRS[k,'V2'])&ALSFRS[k,'V2_days']<=timePoints[i]&ALSFRS[k,'V2_days']>timePoints[i-1]) {
          y[k,i] = 1
          w[k,i] = 1-ALSFRS[k,'V2']/48
          copy_y[k,i] = 0
          copy_w[k,i] = ALSFRS[k,'V2']/48
        }
        if(!is.na(ALSFRS[k,'V3'])&ALSFRS[k,'V3_days']<=timePoints[i]&ALSFRS[k,'V3_days']>timePoints[i-1]) {
          y[k,i] = 1
          w[k,i] = 1-ALSFRS[k,'V3']/48
          copy_y[k,i] = 0
          copy_w[k,i] = ALSFRS[k,'V3']/48
        }
        if(!is.na(ALSFRS[k,'V4'])&ALSFRS[k,'V4_days']<=timePoints[i]&ALSFRS[k,'V4_days']>timePoints[i-1]) {
          y[k,i] = 1
          w[k,i] = 1-ALSFRS[k,'V4']/48
          copy_y[k,i] = 0
          copy_w[k,i] = ALSFRS[k,'V4']/48
        }
      }
    }
    copy_training = training
    x = rbind(x,x)
    y = rbind(y,copy_y)
    w = rbind(w,copy_w)
    w_spare = w
    x = x[rowSums(w_spare)>0,]
    y = y[rowSums(w_spare)>0,]
    w = w[rowSums(w_spare)>0,]
    copy_training = copy_training[rowSums(w_spare[(nrow(training)+1):nrow(w_spare),])>0,]
    
    cvFoldIndex = createFoldsOfData(training, numberOfFolds=3)[[1]]
    copy_cvFoldIndex = createFoldsOfData(copy_training, numberOfFolds=3)[[1]]
    for(cviter in 1:length(cvFoldIndex)) {cvFoldIndex[[cviter]] = c(cvFoldIndex[[cviter]],copy_cvFoldIndex[[cviter]]+nrow(training))}
  }else {
    cvFoldIndex = createFoldsOfData(training, numberOfFolds=3)[[1]]
  }
  
  
  
  
  for(iter in 1:1) {
    
    fitList = DHLRglm(x,y,w,timePoints,cvFoldIndex)
    # res = EMdata(oldx,oldy,oldw,training,fitList,timePoints)
    # x=res$x
    # y=res$y
    # w=res$w
  }
  
  if(F) {
    print('start with KM imputation')
    resKM = EMdataKM(oldx,oldy,oldw,training,kmMod,timePoints)
    x=resKM$x
    y=resKM$y
    w=resKM$w
  }
  
  
  print('start predict')
  #prediction
  survivalFunctionTesting = predictFunctionLRglm(fitList,originalTesting,timePoints)
  survivalFunctionTraining = predictFunctionLRglm(fitList,training,timePoints)
  
  testCurvesToReturn = survivalFunctionTesting
  timesAndCensTest = cbind.data.frame(time = originalTesting$time, delta = originalTesting$delta)
  timesAndCensTrain = cbind.data.frame(time = originalTraining$time, delta = originalTraining$delta)
  trainingCurvesToReturn = survivalFunctionTraining
  
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn,timePoints=timePoints))  
}


  

  
  
