source('Evaluations/EvaluationHelperFunctions.R')
source("ValidateCleanCV/createFoldsAndNormalize.R")
#source("Models/kerasLossFunction.R")

EMdata = function(x,y,w,data,model,timePoints) {
  survivalFunction <- data.frame(matrix(nrow = nrow(x), ncol = ncol(y)))
  hazardFunction <- data.frame(matrix(nrow = nrow(x), ncol = ncol(y)))
  previousTimepointProb = rep(1,nrow(x))
  #probInt = 1 - model %>% predict(x)
  #probInt = probInt[,1:ncol(y)]
  for(i in 1:ncol(y)) {
    if(!is.null(model[[i]])) {
      prob = 1 - predict(model[[i]],newx=x,type='response',s="lambda.min")
      hazardFunction[,i] = prob
    }else {
      prob = rep(1,nrow(x))
    }
    survivalFunction[,i] = prob*previousTimepointProb
    previousTimepointProb = survivalFunction[,i]
  }

  oldy = y
  copy_y = y
  weight = w
  copyWeight = matrix(0,nrow(w),ncol(w))
  for(k in 1:nrow(y)) {
    for(i in 1:ncol(y)) {
      if(w[k,i]<1&data[k,'delta']==0) {
        oldy[k,i] = 0
        copy_y[k,i] = 1
        survivalCurve = c(1,survivalFunction[k,])
        survivalCurveTime = c(0,timePoints)
        if(i>1) {
          a = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,timePoints[i-1])
        }else {
          a = 1
        }
        b = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,timePoints[i])
        c = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,data[k,'time'])
        if(i>1) {
          if(w[k,i-1]<1) {
            survprob = survivalFunction[k,i-1]/c
            #survprob = survivalFunction[k,i-1]+c
            #survprob = a
          }else {
            survprob = 1
          }
        }else {
          survprob = 1
        }
        prob = hazardFunction[k,i]
        #prob = b/a
        # if(i>1) {
        #   prob = 1- ((a/c)-(b/c))
        # }else {
        #   prob = b/c
        # }
        weight[k,i] = survprob*(prob)
        copyWeight[k,i] = survprob*(1-prob)
        if(is.na(weight[k,i])) {
          weight[k,i] = 0
          copyWeight[k,i] = 0
        }
      }
    }
  }
  newx = rbind(x,x)
  newy = rbind(oldy,copy_y)
  neww = rbind(weight,copyWeight)
  copyWeight = matrix(0,nrow(w),ncol(w))
  # newx = rbind(x,x)
  # newy = rbind(y,copy_y)
  # neww = rbind(w,copyWeight)
  
  # newx = newx[rowSums(neww)>0,]
  # newy = newy[rowSums(neww)>0,]
  # neww = neww[rowSums(neww)>0,]
  
  return(list(x=newx,y=newy,w=neww))
}

EMdataKM = function(x,y,w,data,kmMod,timePoints) {
  
  copy_y = y
  weight = w
  copyWeight = matrix(0,nrow(w),ncol(w))
  for(k in 1:nrow(y)) {
    for(i in 1:ncol(y)) {
      if(w[k,i]<1&data[k,'delta']==0) {
        y[k,i] = 0
        copy_y[k,i] = 1
        if(i>1) {
          a = predict(kmMod,timePoints[i-1])
        }else {
          a = 1
        }
        b = predict(kmMod,timePoints[i])
        c = predict(kmMod,data[k,'time'])
        if(i>1) {
          if(w[k,i-1]<1) {
            survprob = predict(kmMod,timePoints[i-1])/c
            #survprob = a/c
          }else {
            survprob = 1
          }
        }else {
          survprob = 1
        }
        #prob = b/a
        if(i>1) {
          prob = 1- ((a/c)-(b/c))
        }else {
          prob = b/c
        }
        weight[k,i] = 1-survprob*(1-prob)
        copyWeight[k,i] = survprob*(1-prob)
      }
    }
  }
  newx = rbind(x,x)
  newy = rbind(y,copy_y)
  neww = rbind(weight,copyWeight)
  
  return(list(x=newx,y=newy,w=neww))
}
predictFunctionLR <- function(fitList,testing,timePoints) {
  numTimepoint = length(timePoints)
  numReturnNA = 0
  numNotDecreasing = 0
  testing[,c('time','delta')] = NULL
  
  previousTimepointProb = rep(1,nrow(testing))
  
  survivalFunction <- data.frame(matrix(ncol = nrow(testing), nrow = numTimepoint))
  for(i in 1:numTimepoint) {
    if(i>length(fitList)){
      survivalFunction[i,] = previousTimepointProb
      break
    }
    fitted = fitList[[i]]
    if(is.null(fitted)) {
      prob = 1
    }else {
      testing = as.matrix(testing)
      prob = 1 - fitted %>% predict(testing)
      #prob = 1 - predict(fitted,newx=testing,type='response',s="lambda.min")
    }
    survivalFunction[i,] = prob*previousTimepointProb
    
    previousTimepointProb = survivalFunction[i,]
  }
  if(numReturnNA>0) {cat('return NA: ',numReturnNA)}
  if(numNotDecreasing>0) {cat('Not decreasing: ',numNotDecreasing)}
  
  colnames(survivalFunction) = 1:nrow(testing)
  return(survivalFunction)
}

predictFunctionLRInt <- function(model,testing,timePoints) {
  numTimepoint = length(timePoints)
  testing[,c('time','delta')] = NULL
  
  survivalFunction <- data.frame(matrix(ncol = nrow(testing), nrow = numTimepoint))
  previousTimepointProb = rep(1,nrow(testing))
  testing = as.matrix(testing)
  probInt = 1 - model %>% predict(testing)
  if(anyNA(probInt)) {print('NaN in keras prediction')}
  #probInt[is.na(probInt)] = 1
  probInt = probInt[,1:length(timePoints)]
  for(i in 1:numTimepoint) {
    prob = probInt[,i]
    survivalFunction[i,] = prob*previousTimepointProb
    previousTimepointProb = survivalFunction[i,]
  }
  colnames(survivalFunction) = 1:nrow(testing)
  return(survivalFunction)
}



prepareDataKeras = function(data,timePoints) {
  
  y = matrix(0,nrow(data),length(timePoints))
  
  for(i in 1:length(timePoints)) {
    singleTime = rep(0,nrow(data))
    singleTime[data$time <= timePoints[i] & data$delta == 1] <- 1
    singleTime[data$time <= timePoints[i] & data$delta == 0] <- NA
    y[,i] =  singleTime
  }
  x = data
  x[,c('time','delta','id')] = NULL
  x = as.matrix(x)
  
  w = matrix(1,nrow(y),ncol(y))
  for(i in 2:length(timePoints)) {
    singleweight = rep(1,nrow(y))
    singleweight[y[,i]==1 & y[,i-1]==0 ] <- 1
    singleweight[y[,i]==1 & y[,i-1]==1 ] <- 0
    singleweight[is.na(y[,i]) & y[,i-1]==0] <- 0.5
    #singleweight[is.na(y[,i]) & y[,i-1]==0] <- 0
    singleweight[is.na(y[,i]) & is.na(y[,i-1])] <- 0
    w[,i] =  singleweight
  }
  y[is.na(y)] = 0
  
  return(list(x=x,y=y,w=w))
}

hazard2survival = function(hazard,timePoints) {
  invHazard = 1-hazard
  survivalFunction <- matrix(ncol = ncol(invHazard), nrow = nrow(invHazard))
  previousTimepointProb = rep(1,nrow(invHazard))
  for(i in 1:ncol(invHazard)) {
    survivalFunction[,i] = invHazard[,i]*previousTimepointProb
    previousTimepointProb = survivalFunction[,i]
  }
  prob = cbind(1,survivalFunction)
  time = c(0,timePoints)
  return(list(prob = prob, time = time))
}

makeMod = function(TestCurves,TrainCurves, timePoints, training, testing) {
  
  #timePoints = fixtime(timePoints)
  survivalFunctionTesting = TestCurves
  survivalFunctionTesting= rbind(rep(1,nrow(testing)),survivalFunctionTesting)
  testCurvesToReturn = cbind(time = c(0,timePoints), survivalFunctionTesting)
  
  #testCurvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilitiesTest) 
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  
  survivalFunctionTraining = TrainCurves
  survivalFunctionTraining= rbind(rep(1,nrow(training)),survivalFunctionTraining)
  trainingCurvesToReturn = cbind(time = c(0,timePoints), survivalFunctionTraining)
  #trainingCurvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilitiesTrain) 
  
  curveCheck(testCurvesToReturn)
  curveCheck(trainingCurvesToReturn)
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
  
}

Evaluation = function(EMMod) {
  survivalPredictionMethod = 'Median'
  bayesConc = Concordance(EMMod, 'None',survivalPredictionMethod)
  bayesBrierInt = BrierScore(EMMod, type = "Integrated", numPoints =  1000, integratedBrierTimes = NULL)
  bayesL1 = L1(EMMod, 'Margin', F,survivalPredictionMethod)
  bayesDcal = DCalibrationCumulative(list(EMMod),10)
  return(list(Conc=bayesConc,BrierInt=bayesBrierInt,L1=bayesL1,Dcal=bayesDcal))
}

logliklyhood = function(xlim,survivalCurves,data,m=100,deltaFraction=0.1) {
  loglik = 0
  delta = data$time*deltaFraction/2
  tstep = xlim/m
  for(t in seq(0,xlim,tstep)) {
    probInt = rep(0,nrow(data))
    for(n in 1:nrow(data)) {
      if(data[n,'time']<=t & data[n,'time']>t-tstep) {
        if(data[n,'delta']==1) {
          pmf = -(predictProbabilityFromCurve(survivalCurves$prob[n,],survivalCurves$time,t+delta)-predictProbabilityFromCurve(survivalCurves$prob[n,],survivalCurves$time,t-tstep-delta))/(delta*2+tstep)
          if(pmf<0) {print('logliklyhood error: negative pmf')}
          loglik = loglik + log(pmf+0.00001)
        }else if(data[n,'delta']==0){
          survFunS = predictProbabilityFromCurve(survivalCurves$prob[n,],survivalCurves$time,t)
          loglik = loglik + survFunS
        }
      }
    }
  }
  return(loglik)
}

cvGeneral = function(data,timePoints,numberOfFolds=3) {
  #lambdaList = c(0.001,0.01,0.1,0.5,1,10)
  lambdaList = c(50,40,30,20,10,5)
  lossList = rep(0,length(lambdaList))
  cvFoldIndex = createFoldsOfData(data, numberOfFolds=numberOfFolds)[[1]]
  cat('internal cross validation: ')
  for(lambdaIter in 1:length(lambdaList)) {
    cat(lambdaList[lambdaIter]);cat(' ');
    for(cvIter in 1:length(cvFoldIndex)) {
      quantileVals = seq(0,1,length.out = lambdaList[lambdaIter]+2)[-c(1,lambdaList[lambdaIter]+2)]
      timePoints = unname(quantile(data[cvFoldIndex[[cvIter]],]$time, quantileVals))
      timePoints = timePoints[!duplicated(timePoints)]
      print(lambdaList[lambdaIter])
      res = prepareDataKeras(data,timePoints)
      x=res$x
      y=res$y
      w=res$w
      cv_x = x[-cvFoldIndex[[cvIter]],]
      cv_y = y[-cvFoldIndex[[cvIter]],]
      cv_w = w[-cvFoldIndex[[cvIter]],]
      cv_x_v = x[cvFoldIndex[[cvIter]],]
      cv_y_v = y[cvFoldIndex[[cvIter]],]
      cv_w_v = w[cvFoldIndex[[cvIter]],]
      
      allHazard = colSums(cv_y*cv_w)/colSums(cv_w)
      my_regularizer_wrapper_cv <- custom_metric("reg", function(x){my_regularizer(x, lambda1=0.1,lambda2=0.1,weights=w)})
      my_bias_regularizer_wrapper <- custom_metric("regb", function(x){my_bias_regularizer(x, allHazard=log(allHazard))})
      
      model_cv <- keras_model_sequential()
      model_cv %>%  layer_dense(units=2*ncol(cv_y), activation='sigmoid', input_shape=c(ncol(cv_x)),use_bias=TRUE,
                                #bias_initializer = initializer_random_normal(-2.19,0.5),
                                #kernel_regularizer = regularizer_l1(l = 0.5),
                                kernel_regularizer = my_regularizer_wrapper_cv,
                                bias_regularizer = my_bias_regularizer_wrapper
      )
      #summary(model)
      optimizerSGD = optimizer_sgd(
        lr = 1.0,
        momentum = 0.1,
        decay = 0.0,
        nesterov = FALSE,
        clipnorm = 1,
        clipvalue = 1
      )
      model_cv %>% compile(
        loss = customBCE,
        optimizer = optimizerSGD,
        metrics = NULL
        #sample_weight_mode='temporal',
        #weighted_metrics = 'binary_accuracy'
      )
      early_stopping = callback_early_stopping(monitor='loss', patience=300, verbose=2)
      reduceLearningRate = callback_reduce_lr_on_plateau(monitor='loss', patience=10, factor=0.9,min_lr=0.01, verbose=0)
      
      fitted <- model_cv %>% fit(
        cv_x, cbind(cv_y,cv_w),
        epochs = 6000, batch_size = 16,
        validation_split = 0,
        #validation_data = list(cv_x_v,cbind(cv_y_v,cv_w_v)),
        shuffle=T,
        verbose = 0,
        callbacks = list(early_stopping,reduceLearningRate)
      )
      
      probInt = model_cv %>% predict(cv_x_v)
      k_clear_session()
      
      probInt = probInt[,1:ncol(cv_y_v)]
      survivalCurves = hazard2survival(probInt,timePoints)
      
      loss = logliklyhood(xlim=quantile(data$time,0.8),survivalCurves,data[cvFoldIndex[[cvIter]],],m=100)
      lossList[lambdaIter] = lossList[lambdaIter]+loss
      #print(loss)
    }
    print(lossList)
  }
  bestLambda = lambdaList[which.max(lossList)]
  
  return(bestLambda)
}