input = read.csv("pretrained/new_input.csv")
# normalize input variables
normalized_input = fast_scale(input, scales=scales,way="scale",verbose=F)
# load pre-trained mtlr model
mtlrMode = readRDS("pretrained/mtlrMod.RDS")
# predict
result_curves = predict(mtlrMode,normalized_input)
# convert time from days to months
result_curves$time = result_curves$time/30.4
# plot curves
plotSurvivalCurves(result_curves, c(1:10), title='Individual Survival Distribution',x_title="Month")
# save curves as csv
write.csv(result_curves,"pretrained/result_curves.csv")
