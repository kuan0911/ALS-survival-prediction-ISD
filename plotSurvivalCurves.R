#### File Information #####################################################################################################################
#File Name: plotSurvivalCurves.R
#Date Created: June 20, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#This file is used to plot the survival curves of each survival model.

### Functions #############################################################################################################################

## Function 1: plotSurvivalCurves(survivalCurves, indexToPlot = 1, color = c(), xlim = c())

#Inputs: 
#   survivalCurves: A matrix of survival curves. The first column must be the survival times and all other columns should represent the 
#                   survival probabilites at each of those times.
#   indexToPlot: The indexs to plot i.e. which patients curve to plot.
#   color: The survival curve color. This must match the length of then number of indexs to plot.
#   xlim: The interval to plot on the x-axis.

# Output: The plot of the survival curves.

# Usage: This function should be used to plot survival curves.

### Code ##################################################################################################################################
#Library Dependencies:
#For the plots:
library(ggplot2)
#For shaping the data into long form:
library(reshape2)
library(ggrepel)
plotSurvivalCurves = function(survivalCurves, indexToPlot = 1, color = c(), xlim = c(), title='',plot_points = data.frame(),x_title='time') {
  colorOK = T
  if(length(color) == 0)
    colorOK = F
  else if(length(color) != length(indexToPlot)){
    warning("If you would like to select custom colors please make sure the number of colors
            matches the number of curves.")
    colorOK =F 
  }
  time = survivalCurves$time
  curves = survivalCurves[,indexToPlot +1,drop=F]
  maxPlotTimes = ifelse(length(xlim)==2,xlim[2],max(time))
  plotTimes = seq(min(time),maxPlotTimes, length.out = length(time)*100)
  plotProbs = as.data.frame(sapply(curves,
                                   function(curve){
                                     curve = ifelse(curve < 1e-20,0,curve)
                                     survivialSpline = splinefun(time, curve, method = "hyman")
                                     return(pmax(predictProbabilityFromCurve(curve,time,plotTimes),0))
                                   }
                                   #predictProbabilityFromCurve(curve,time,plotTimes)
  ))
  data = cbind.data.frame(plotTimes,plotProbs)
  longFormData = melt(data,measure.vars = names(data)[-1], variable.name = "Index")
  plot = ggplot(data = longFormData, aes(x = plotTimes,y = value, colour = Index))+
    geom_line(linewidth = 1.0)
  if(colorOK)
    plot = plot + scale_color_manual(values = color) 
  if(length(xlim)==2){
    plot = plot+ xlim(c(xlim[1],xlim[2]))
  }
  plot = plot +scale_y_continuous( limits= c(0,1),breaks = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))+
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(text = element_text(size=15),
          axis.title = element_text(size = 12),
          axis.title.x = element_text(margin = margin(t = 15)),
          axis.title.y = element_text(margin = margin(r = 15))) + 
    labs(y = "Survival Probaility",x = x_title,title=title )
  #If we have too many survival curves then the legend takes up the whole plot. We have an if to catch this and remove the legend.
  if(length(indexToPlot) > 15){
    plot = plot + theme(legend.position = "None")
  }
  if(nrow(plot_points) > 0) {
    shape_vec = rep(19,nrow(plot_points))
    if(!is.null(plot_points$filled)) {shape_vec[!plot_points$filled] = 1}
    colour_vec = rep('black',nrow(plot_points))
    if(!is.null(plot_points$colour)) {colour_vec = plot_points$colour}
    plot_points$label = paste('(',plot_points$x,', ',plot_points$y,')',sep='')
    plot = plot + geom_point(data=plot_points,aes(x=x, y=y),colour = colour_vec, size = 2,shape=shape_vec,stroke=1.5) + 
      geom_text_repel(data=plot_points,aes(x=x, y=y,label = label),nudge_x=6,colour=colour_vec)
  }
  return(plot)
}




