#-------------------
#cross validation
#-------------------
source("Models/DHLR.R")

DHupper = function(training,testing,ALSFRS=NULL){
  print('start internal cross validation')
  #bestC1 = internalCV_BayesianNet(training,5,10)
  bestC1 = 1
  print('bestC1')
  print(bestC1)
  
  #m = floor(sqrt(nrow(training[training$delta==1,]))+1)
  m = 3
  quantileVals = seq(0,1,length.out = m+2)[-c(1,m+2)]
  #quantileVals = seq(0,1,length.out = m+2)[-1]
  timePoints = unname(quantile(training[training$delta==1,]$time, quantileVals))
  
  # kmMod = prodlim(Surv(time,delta)~1, data = training)
  # step = max(training$time)/500
  # timePoints = kmTimesplitV2(m,kmMod,training,step=step)
  #timePoints = timePoints[timePoints<max(training[training$delta==1,]$time)+1]
  #fillup = seq(max(timePoints),max(training$time),(tail(timePoints, n=2)[2]-tail(timePoints, n=2)[1])*2)
  #timePoints = c(timePoints,fillup)
  #timePoints = c(timePoints,max(training$time))
  timePoints = timePoints[!duplicated(timePoints)]
  m = length(timePoints)

  mod = DHLR(training, testing,timePoints,debug=T,ALSFRS=ALSFRS)
  
  timePoints = mod$timePoints
  survivalFunctionTesting = mod$TestCurves
  survivalFunctionTesting= rbind(rep(1,nrow(testing)),survivalFunctionTesting)
  testCurvesToReturn = cbind(time = c(0,timePoints), survivalFunctionTesting)
  
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)

  survivalFunctionTraining = mod$TrainCurves
  survivalFunctionTraining= rbind(rep(1,nrow(training)),survivalFunctionTraining)
  trainingCurvesToReturn = cbind(time = c(0,timePoints), survivalFunctionTraining)
  
  curveCheck(testCurvesToReturn)
  curveCheck(trainingCurvesToReturn)
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

curveCheck = function(survivalFunction) {
  if(anyNA(survivalFunction)) {
    print('survival function has missing value')
    return()
  }
  for(k in 2:ncol(survivalFunction)) {
    for(i in 2:nrow(survivalFunction)) {
      if(survivalFunction[(i-1),k]<survivalFunction[i,k]) {
        print('survival function check fail')
      }
    }
  }
}

