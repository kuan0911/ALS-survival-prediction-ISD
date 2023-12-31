#### File Information #####################################################################################################################

#File Name: EvaluationHelperFunctions.R
#Date Created: May 28th, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#This file contains a number of methods useful for model predictions/model evaluations.

### Functions #############################################################################################################################

## Function 1: predictProbabilityFromCurve(survivalCurve,predictedTimes, timeToPredict)

# Inputs:
#   survivalCurve: A vector of survival probabilities, e.g. (1, .9, .8, 0)
#   predictedTimes: The times corresponding to each survival probability, e.g. (0, 10, 20, 40)
#   timeToPredict: The time at which to evaluate the survival curve

# Output: The probability of a survival curve given a specified time.

# Usage: Given a survival curve we would like to extract a survival probability for a certain time. To do this we need to fit the model
#        to a spline an then extrac the survival probability.


## Function 2: predictMeanSurvivalTimeSpline(survivalCurve,predictedTimes)

# Inputs: See predictProbabilityFromCurve()

# Output: The mean survival time given a survival curve.

# Usage: Given a survival curve we extend the curve to zero using a linear line and then integrate to get the mean survival time.


## Function 3: predictMeanSurvivalTimeSpline(survivalCurve,predictedTimes)

# Inputs: See predictProbabilityFromCurve()

# Output: The median survival time given a survival curve.

# Usage: Given a survival curve we extend the curve to zero using a linear line and then get the median time. (Note the linear extension
#        does not need to go all the way to zero.)


### Code ##################################################################################################################################
#Library Dependencies: None.
#We need some type of predict function for survival curves - here we build a spline to fit the survival model curve. This spline is 
#the montotone spline using the hyman filtering of the cubic Hermite spline method,
#see https://en.wikipedia.org/wiki/Monotone_cubic_interpolation. Also see help(splinefun).
#Note that we make an alteration to the method because if the last two time points
#have the same probability (y value) then the spline is constant outside of the training data. We need this to be a decreasing function
#outside the training data so instead we take the linear fit of (0,1) and the last time point we have (p,t*) and then apply this linear
#function to all points outside of our fit.
predictProbabilityFromCurve = function(survivalCurve,predictedTimes, timeToPredict){
  spline = splinefun(predictedTimes, survivalCurve, method = "hyman")
  maxTime = max(predictedTimes)
  slope = (1-spline(maxTime))/(0 - max(predictedTimes))
  predictedProbabilities = rep(0, length(timeToPredict))
  linearChange = which(timeToPredict > maxTime)
  if(length(linearChange) > 0){
    predictedProbabilities[linearChange] = pmax(1 + timeToPredict[linearChange]*slope,0)
    predictedProbabilities[-linearChange] = spline(timeToPredict[-linearChange])
  }
  else{
    predictedProbabilities = spline(timeToPredict)
  }
  return(predictedProbabilities)
}

#We calculate the mean and median survival times assuming a monotone spline fit of the survival curve points.
predictMeanSurvivalTimeSpline = function(survivalCurve, predictedTimes){
  #If all the predicted probabilities are 1 the integral will be infinite. For this reason we slightly decrease the 
  #last value.
  if(all(survivalCurve==1)){
    return(Inf)
  }
  spline = splinefun(predictedTimes, survivalCurve, method = "hyman")
  maxTime = max(predictedTimes)
  slope = (1-spline(maxTime))/(0 - max(predictedTimes))
  zeroProbabilitiyTime = min( predictedTimes[which(survivalCurve ==0)], maxTime + (0-spline(maxTime))/slope)
  splineWithLinear = function(time) ifelse(time < maxTime, spline(time),1 + time*slope)
  area = integrate(splineWithLinear,0, zeroProbabilitiyTime,subdivisions = 1000,rel.tol = .0001)[[1]]
  return(area)
}

predictMedianSurvivalTimeSpline = function(survivalCurve, predictedTimes){
  #If all the predicted probabilities are 1 the integral will be infinite.
  if(all(survivalCurve==1)){
    return(Inf)
  }
  if(all(survivalCurve<0.001)){
    return(0)
  }
  spline = splinefun(predictedTimes, survivalCurve, method = "hyman")
  minProb = min(spline(predictedTimes))
  
  if(is.na(minProb)){return(0)}
  
  #cut = max(survivalCurve)/2
  cut=0.5
  
  if(minProb < cut){
    maximumSmallerThanMedian = predictedTimes[min(which(survivalCurve <cut))]
    minimumGreaterThanMedian = predictedTimes[max(which(survivalCurve >cut))]
    splineInv = splinefun(spline(seq(minimumGreaterThanMedian, maximumSmallerThanMedian, length.out = 1000)),
                          seq(minimumGreaterThanMedian, maximumSmallerThanMedian, length.out = 1000))
    medianProbabilityTime = splineInv(cut)
  }
  else{
    maxTime = max(predictedTimes)
    slope = (max(survivalCurve)-spline(maxTime))/(min(predictedTimes) - max(predictedTimes))
    medianProbabilityTime = maxTime + (cut-spline(maxTime))/slope
  }
  return(medianProbabilityTime)
}


