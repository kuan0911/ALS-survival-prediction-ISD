#### File Information #####################################################################################################################
#File Name: AcceleratedFailureTime.R
#Date Created: May 28, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################

#This file is used to run the Accelerated Failure Time (AFT) model for some specified distribution to create individual survival curves.
#This implimentaiton only supports the following distributions: weibull, exponential, lognormal, gaussian, loglogistic, and logistic.

### Functions #############################################################################################################################

## Function 1: AFT(training, testing, AFTDistribution)

#Inputs:
#   training: The training dataset (after normalization and imputation).
#   testing: The testing dataset (after normalization and imputation).
#   AFTDistribution: The distribution to use for the AFT model.

# Output: A list of 4 items:(1) TestCurves - The survival curves for the testing set.
#                           (2) TestData - The censor/death indicator and event time for the testing set. 
#                           (3) TrainData - The censor/death indicator and event time for the training set. 
#                           (4) TrainCurves - The survival curves for the training set.

# Usage: Train and evaluate the AFT model.


## Function 2:  survfunc(AFTMod, newdata, t)

# Inputs:
#   AFTMod: The AFT model returned by survreg().
#   newdata: The data on which to evaluate the AFT model.
#   t: The times at which to evaluate the AFT model.

# Output: The survival curves for each patient.

# Usage: Evaluate the AFT model to get survival probabilites. (Helper function for AFT).


### Code #############################################################################################################################
#Library Dependencies:
#survival is needed to get survreg for the AFT model.
library(survival)
library(iregnet)

AFT = function(training, testing, AFTDistribution){
  #Sometimes the AFT model will fail to converge (though rare) we want to catch this.
  tryCatch({
    if(ncol(training)<=3) {
      training$dummy = runif(nrow(training));
      testing$dummy = runif(nrow(testing));
    }
    X <- as.matrix(training[,!colnames(training) %in% c('time','delta')])
    y <- matrix(0, nrow(training), 2)
    y[,1] = training$time
    y[,2] = training$time
    y[training$delta==0,2] = Inf
    AFTMod <- cv.iregnet(X, y, family = AFTDistribution, alpha=0)
    
    #newdata = newdata[,!colnames(newdata) %in% c('time','delta')]
    #lp <- predict(AFTMod, newx = as.matrix(newdata),type = 'response',lambda=0)
    
    #allvars <- as.matrix(training)
    #allvars = allvars[,!colnames(allvars) %in% c('time','delta')]
    #ridge.formula <- as.formula(paste("Surv(time,delta) ~. + ridge(allvars,theta = 10)"))
    
    #AFTMod = survreg(ridge.formula, data = training, dist = AFTDistribution,control = list(maxiter=5000))
    #names(AFTMod$coefficients) = c("(Intercept)",colnames(allvars))
    #AFTMod = survreg(Surv(time,delta)~Age+Sexmale, data = training, dist = AFTDistribution,control = list(maxiter=1000))
    
    trainingTimes = sort(unique(training$time))
    if(length(trainingTimes)>20) {
      trainingTimes = trainingTimes[seq(1,length(trainingTimes),floor(length(trainingTimes)/20))]
    }
    if(0 %in% trainingTimes){
      timesToPredict = trainingTimes
    } else {
      timesToPredict = c(0,trainingTimes)
    }

    survivalCurves = survfunc(AFTMod, newdata = testing, t = timesToPredict)
    survivalCurvesTrain = survfunc(AFTMod, newdata = training, t = timesToPredict)
  },
  error = function(e) {
    message(e)
    warning("AFT failed to converge.")
  })
  if(!exists("AFTMod") | !exists("survivalCurves")){
    return(NA)
  }

  probabilities = survivalCurves$sur
  probabilitiesTrain = survivalCurvesTrain$sur
  #Since survfunc returns survival probabilities with the first time point for every individual (ordered by how the testing individuals)
  #were passed in, we can simply fill a matrix by row to have each individual curve be a column. This can be verified by checking
  #the survival probabilties (sur) for any ID_SurvivalCurves against any column, e.g. the survival probabilities for ID_SurvivalCurves == 2,
  #correspond to the probabilites found in the second column of the matrix below.
  probabilityMatrix = matrix(probabilities, ncol = nrow(testing),byrow = T)
  probabilityTrainMatrix = matrix(probabilitiesTrain, ncol = nrow(training),byrow = T)
  
  curvesToReturn = cbind.data.frame(time = timesToPredict, probabilityMatrix)
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  trainingCurvesToReturn = cbind.data.frame(time = timesToPredict, probabilityTrainMatrix)
  return(list(TestCurves = curvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

#The following was taken and altered from http://rstudio-pubs-static.s3.amazonaws.com/161203_6ee743eb28df4cd68089a519aa705123.html.
#This code is used to pass in a time point and get the predicted probability. The predict function for survreg objects only return 
#times from probabilities so we need to reverse enginner this using the code below.
survfunc = function (object, t, newdata, name = "t") {
  #Altered from origina: I am going to add an ID to every row so we can retrieve the individuals easily from the output.
  #I gave a weird ID variable name so if the original data came in with a variable ("ID") it won't break our system.
  newdata = newdata[,!colnames(newdata) %in% c('time','delta')]
  if(any(class(object)=='survreg')) {
    newdata$ID_SurvivalCurves = 1:nrow(newdata)
  }
  newdata <- do.call(rbind, rep(list(newdata), length(t)))
  t <- rep(t, each = nrow(newdata)/length(t))
  if (all(class(object) != "survreg") && all(class(object) != "iregnet")) 
    stop("not a survreg object or iregnet object")
  if(any(class(object)=='iregnet')) {
    lp <- predict(object, newx = as.matrix(newdata),type = 'min')[,1]
  }else if(any(class(object)=='survreg')) {
    lp <- predict(object, newdata = newdata, type = "lp")
  }
  
  if(any(class(object)=='cv.iregnet')) {
    dist = object$family
    scale = object$scale[which.min(object$loglik)]
  }else if(any(class(object)=='iregnet')) {
    dist = object$call$family
    scale = object$scale[length(object$scale)]
  }else if(any(class(object)=='survreg')) {
    dist = object$dist
    scale = object$scale
  }
  if (dist %in% c("weibull", "exponential")) {
    newdata$pdf <- dweibull(t, 1/scale, exp(lp))
    newdata$cdf <- ifelse(t == 0,0,
                          ifelse(is.nan(pweibull(t, 1/scale, exp(lp))),1,pweibull(t, 1/scale, exp(lp))))
    newdata$haz <- exp(dweibull(t, 1/scale, exp(lp), 
                                log = TRUE) - pweibull(t, 1/scale, exp(lp), 
                                                       lower.tail = FALSE, log.p = TRUE))
  }
  else if (dist == "lognormal") {
    newdata$pdf <- dlnorm(t, lp, scale)
    newdata$cdf <- plnorm(t, lp, scale)
    newdata$haz <- exp(dlnorm(t, lp, scale, log = TRUE) - 
                         plnorm(t, lp, scale, lower.tail = FALSE, log.p = TRUE))
  }
  else if (dist == "gaussian") {
    newdata$pdf <- dnorm(t, lp, scale)
    newdata$cdf <- pnorm(t, lp, scale)
    newdata$haz <- exp(dnorm(t, lp, scale, log = TRUE) - 
                         pnorm(t, lp, scale, lower.tail = FALSE, log.p = TRUE))
  }
  else if (dist == "loglogistic") {
    newdata$pdf <- dlogis(log(t), lp, scale)/t
    newdata$cdf <- plogis(log(t), lp, scale)
    newdata$haz <- exp(dlogis(log(t), lp, scale, log = TRUE) - 
                         log(t) - plogis(log(t), lp, scale, lower.tail = FALSE, 
                                         log.p = TRUE))
  }
  else if (dist == "logistic") {
    newdata$pdf <- dlogis(t, lp, scale)
    newdata$cdf <- plogis(t, lp, scale)
    newdata$haz <- exp(dlogis(t, lp, scale, log = TRUE) - 
                         dlogis(t, lp, scale, lower.tail = FALSE, log.p = TRUE))
  }
  else {
    stop("unknown distribution")
  }
  newdata$sur <- 1 - newdata$cdf
  newdata[name] <- t
  return(newdata)
}

