#### File Information #####################################################################################################################
#File Name: MTLR.R
#Date Created: October 23, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#The algorithim here is based Multi-task Logistic Regression (MTLR) for which the original paper can be found here:
#https://papers.nips.cc/paper/4210-learning-patient-specific-cancer-survival-distributions-as-a-sequence-of-dependent-regressors

#This file is used to run MTLR to generate individual survival curves. Here we are using R's base `optim` function for gradient descent
#(L-BFGS-B) specifically. In order to optimize the code we have written the gradient and objective function calculation in Rcpp files
#(see sourceCpp command below).

### Functions #############################################################################################################################

## Function 1: MTLR(training, testing, C1 = NULL, numFolds = 5)

#Inputs:
#   training: The training dataset (after normalization and imputation).
#   testing:  The testing dataset (after normalization and imputation).
#   C1:       Regularization parameter. Default is NULL. If NULL the parameter will be selected using internal cross validation.
#   numFolds: Number of folds for internal cross validation.

# Output: A list of 4 items:(1) TestCurves - The survival curves for the testing set.
#                           (2) TestData - The censor/death indicator and event time for the testing set. 
#                           (3) TrainData - The censor/death indicator and event time for the training set. 
#                           (4) TrainCurves - The survival curves for the training set.

# Usage: Train and evaluate the MTLR model.


## Function 2: internalCV_MTLR(training, numFolds)

# Inputs: See MTLR().

# Output: The optimal regularization parameter (C1).

# Usage: Perform internal cross validation for selecting the regularization parameter. The regularization parameter is chosen by finding
#        the value which optimizes the average log likelihood loss. (Helper function for MTLR).


## Function 3: avgLogLikLoss(params, dat, timePoints)

# Inputs:
#   params: The feature weights.
#   dat: The testing dataset.
#   timePoints: The time points for evaluating MTLR.

# Output: The average log-likelihood loss.

# Usage: Calculate the average log-likelhood loss on the testing data using the learned feature weights. 
#        (Helper function for internalCV_MTLR).

### Code ##################################################################################################################################
#Library Dependencies:
#We use Rcpp to source the Rcpp files used for MTLR.
library(Rcpp)

#File Dependencies
#Cpp files for MTLR:
library(MTLR)
source("ValidateCleanCV/createFoldsAndNormalize.R")
#For computing the average log likelihood loss we will need predictProbabilityFromCurve()
source("Evaluations/EvaluationHelperFunctions.R")

MTLR = function(training, testing, C1 = NULL, numFolds = 5, ALSFRS_training = NULL){
  if(is.null(C1)){
    #nintervals <- ceiling(sqrt(((5-1)/5)*nrow(training)))
    m = floor(sqrt(nrow(training)*4/5)+1)
    #m = 10
    quantileVals = seq(0,1,length.out = m+2)[-c(1,m+2)]
    timePoints = unname(quantile(training$time, quantileVals))
    
    #extend curve plot
    #timePoints = c(timePoints,max(c(training$time,testing$time)))
    
    # kmMod = prodlim(Surv(time,delta)~1, data = training)
    # step = max(training$time)/500
    # timePoints = kmTimesplitV2(m,kmMod,training,step=step)
    # print(timePoints)
    
    time_points = timePoints[!duplicated(timePoints)]
    m = length(time_points)
    C1 = mtlr_cv(Surv(time,delta)~.,data=training, loss= "conc", C1_vec = c(0.001,0.01,0.1,1,10,100,1000,10000)
, train_biases = F,train_uncensored = F,time_points = time_points)$best_C1
  }
  
  
  if(!is.null(ALSFRS_training)) {
    print('using ALSFRS')
    riskScore_mat_training = matrix(NA,nrow=nrow(training),ncol=ncol(ALSFRS_training$time))
    riskScore_mat_testing = matrix(NA,nrow=nrow(testing),ncol=ncol(ALSFRS_training$time))
    riskScore_data = training
    for(j in 1:ncol(ALSFRS_training$time)) {
      riskScore_data$time = ALSFRS_training$time[,j]
      riskScore_data$delta = ALSFRS_training$delta[,j]
      riskScore_C1 = mtlr_cv(Surv(time,delta)~.,data=riskScore_data, loss= "conc", C1_vec = c(0.001,0.01,0.1,1,10,100,1000,10000)
                   , train_biases = F,train_uncensored = F,time_points = time_points)$best_C1
      riskMod = mtlr(Surv(time,delta)~., data = riskScore_data, C1=riskScore_C1, time_points = time_points, train_biases = F, train_uncensored = F)
      riskScoreCurvesTrain = predict(riskMod)[-1]
      riskScorePredictedTimes = predict(riskMod)$time
      average_training = unlist(lapply(1:nrow(training),function(index) predictMedianSurvivalTimeSpline(riskScoreCurvesTrain[,index],riskScorePredictedTimes)))
      riskScore_mat_training[,j] = average_training

      riskScoreCurvesTest = predict(riskMod,testing)[-1]
      average_testing = unlist(lapply(1:nrow(testing),function(index) predictMedianSurvivalTimeSpline(riskScoreCurvesTest[,index],riskScorePredictedTimes)))
      riskScore_mat_testing[,j] = average_testing
    }
    training = cbind(training,scale(riskScore_mat_training))
    testing = cbind(testing,scale(riskScore_mat_testing))
    
    # split_index = nrow(training)
    # validateData = FeatureSelection(rbind(training,testing), type = "UniCox")
    # training = validateData[1:split_index,]
    # testing = validateData[split_index:nrow(validateData),]
    #print(training)
    #testing = testing[,colnames(training)]
  }
  
  C1 = mtlr_cv(Surv(time,delta)~.,data=training, loss= "conc", C1_vec = c(0.001,0.01,0.1,1,10,100,1000,10000)
               , train_biases = F,train_uncensored = F,time_points = time_points)$best_C1
  #print(paste0('MTLR C1: ', C1))
  mod = mtlr(Surv(time,delta)~., data = training, C1=C1, time_points = time_points, train_biases = F, train_uncensored = F)
  testCurvesToReturn = predict(mod,testing)
  #testCurvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilitiesTest) 
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  trainingCurvesToReturn = predict(mod)
  #trainingCurvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilitiesTrain) 
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

internalCV_MTLR = function(training, numFolds){
  foldedData = createFoldsOfData(training, numFolds)[[2]] 
  resultsMatrix = matrix(rep(0,numFolds*7), ncol = 7,nrow =5) #7 vals of C1
  
  #Much of this code is duplicated from MTLR(). See comments in that function.
  for(i in 1:numFolds){
    trainingFold = foldedData[[1]][[i]]
    testingFold = foldedData[[2]][[i]]
    
    trainingFold = trainingFold[order(trainingFold$delta),]
    testingFold = testingFold[order(testingFold$delta),]
    
    m = floor(sqrt(nrow(training))+1)
    quantileVals = seq(0,1,length.out = m+2)[-c(1,m+2)]
    timePoints = unname(quantile(trainingFold$time, quantileVals))
    timePoints = timePoints[!duplicated(timePoints)]
    d = as.matrix(trainingFold[,-c(1,2)])
    dAsZero = matrix(0,ncol = ncol(d), nrow = nrow(d))
    yval = matrix(1 -Reduce(c,Map(function(x) trainingFold$time > timePoints[x],1:length(timePoints))), ncol = nrow(trainingFold), byrow = T)
    
    #Train biases first and then all parameters.
    resultVec = c()
    for(C1 in c(0.001,0.01,0.1, 1, 10,100,1000)){
      biasPar = optim(par = rep(0,length(timePoints)*(ncol(d) +1)),fn = mtlr_objVal,gr = mtlr_grad, yval = yval, featureVal = dAsZero,C1=C1, delta = sort(trainingFold$delta), 
                      method = "L-BFGS-B", lower = -20,upper=20,control=c(maxit = 5000, factr = .45036e11))
      
      allParamsUnc = optim(par = biasPar$par,fn = mtlr_objVal,gr = mtlr_grad, yval = yval, featureVal = d,C1=C1,delta = rep(1,nrow(trainingFold)),
                           method = "L-BFGS-B", lower = -20,upper=20,control=c(maxit = 5000, factr = .45036e11))
      
      allParams = optim(par = allParamsUnc$par,fn = mtlr_objVal,gr = mtlr_grad, yval = yval, featureVal = d,C1=C1,delta = sort(trainingFold$delta),
                        method = "L-BFGS-B", lower = -20,upper=20,control=c(maxit = 5000, factr = .45036e11))
      
      resultVec = c(resultVec,avgLogLikLoss(allParams$par,testingFold,timePoints))
    }
    resultsMatrix[i,] = resultVec
  }
  meanResults = apply(resultsMatrix, 2, mean)
  bestC1 = c(0.001,0.01,0.1, 1, 10, 100,1000)[which.min(meanResults)]
  return(bestC1)
}

avgLogLikLoss = function(params, dat, timePoints){
  #For the log-likelihood loss we need to compute losses differently for censored and uncensored patients.
  #For censored patients the loss will correspond the the (log) survival probability assigned by the model at the time of censoring.
  #For uncensored patients, we will consider the log of the probability assigned to the time interval when the patient died. 
  #Then we take the negative of this loss and thus would like to minimize the loss.
  d = as.matrix(dat[,-c(1,2)])
  survivalCurves = mtlr_predict(params,d)
  NCens = sum(1-  dat$delta)
  
  logloss = 0
  #Censored patients
  censorTimes = dat$time[1:NCens]
  probAtCensorTime = sapply(seq_along(censorTimes),
                            function(index) predictProbabilityFromCurve(survivalCurves[,index],
                                                                        timePoints,
                                                                        censorTimes[index]))
  logloss = logloss - sum(log(probAtCensorTime + 1e-5))
  
  #Uncensored patients
  deathTimes = dat$time[(NCens+1):nrow(dat)]
  uncenSurvival = survivalCurves[,(NCens+1):nrow(dat),drop=FALSE]
  uncenSurvival = rbind(1,uncenSurvival,0)
  pmfProbs = -diff(uncenSurvival)
  indexRow = sapply(deathTimes, function(x) findInterval(x, timePoints)) + 1
  indexCol = 1:(nrow(dat) - NCens)
  indexMat = matrix(c(indexRow,indexCol),ncol = 2)
  probs = pmfProbs[indexMat]
  logloss = logloss  - sum(log(probs))
  
  return(logloss/nrow(dat))
}