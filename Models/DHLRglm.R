DHLRglm = function(x_all,y_all,w_all,timePoints,cvFoldIndex) {
  fitList <- vector("list", length(timePoints))
  for(i in 1:length(timePoints)) {
    cat(i);cat(' ');
    
    x = x_all
    y = y_all[,i]
    w = w_all[,i]
    
    if(sum(y == 0 & w>0.5)<5 | sum(y == 1 & w>0.5)<5) {
      fitList[[i]] = NULL
      cat('skip ')
    }else {
      foldid = rep(NA,nrow(x))
      for(cvIter in 1:length(cvFoldIndex)) {foldid[cvFoldIndex[[cvIter]]]=cvIter}
      fitted = cv.glmnet(x,y,family='binomial', weights=w, foldid=foldid)
      #fitted <- glmnet(x,y,family='binomial',weights=w,lambda=0)
      fitList[[i]] <- fitted
    }
  }
  return(fitList)
}

predictFunctionLRglm <- function(fitList,testing,timePoints) {
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
      prob = 1 - predict(fitted,newx=testing,type='response',s="lambda.min")
    }
    survivalFunction[i,] = prob*previousTimepointProb
    
    previousTimepointProb = survivalFunction[i,]
  }
  if(numReturnNA>0) {cat('return NA: ',numReturnNA)}
  if(numNotDecreasing>0) {cat('Not decreasing: ',numNotDecreasing)}
  
  colnames(survivalFunction) = 1:nrow(testing)
  return(survivalFunction)
}