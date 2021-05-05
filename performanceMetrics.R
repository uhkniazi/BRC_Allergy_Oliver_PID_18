# Name: performanceMetrics.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 4/5/2021
# Desc: report performance metrics for imported prediction data

library(ROCR)

dfData = read.csv('results/paper/predictions.csv', header=T)
ivPredict = dfData$Predicted.PA.Probability
ivTruth = dfData$Actual.Group == 'PA'

p = prediction(ivPredict, ivTruth)
s = performance(p, 'sens')
sp = performance(p, 'spec')
ppv = performance(p, 'ppv')
npv = performance(p, 'npv')

dfMetrics = round(data.frame(Cutoff=s@x.values[[1]], Sensitivity=s@y.values[[1]],
                       Specificity=sp@y.values[[1]],
                       PPV=ppv@y.values[[1]], NPV=npv@y.values[[1]]), 3)
write.csv(dfMetrics, file='results/paper/model_2_performance.csv')
