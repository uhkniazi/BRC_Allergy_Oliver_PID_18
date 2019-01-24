# Name: 05_serologyVariableSelection.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 22/01/2019
# Desc: data modelling and variable selection for serology data

############ data loading and cleaning
source('header.R')

## load the data
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 37)')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
df = read.csv(n, header=T)

# close connection after getting data
dbDisconnect(db)

dim(df)
df$Patient = factor(gsub(' ', '', as.character(df$Patient)))
df$Allergic.Status = factor(gsub(' ', '', as.character(df$Allergic.Status)), levels = c('PS', 'PA'))
df = droplevels.data.frame(df)
str(df)
data.frame(colnames(df))
df = df[,-c(18:21)]
df = na.omit(df)
df = droplevels.data.frame(df)

## make count matrix
data.frame(colnames(df))
mData = as.matrix(df[,c(4:17)])
colnames(mData)
dfSample = df
rownames(mData) = as.character(dfSample$Patient)
str(dfSample)

lData.train = list(data=mData, covariates=dfSample)
rm(df)
############ end data loading

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

########################## perform a random forest step
dfData = data.frame(log(lData.train$data+1e-4))
fGroups = lData.train$covariates$Allergic.Status

set.seed(123)
oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100)
save(oVar.r, file='temp/oVar.r_serology.rds')

plot.var.selection(oVar.r)

######################## Stan section for binomial regression approach
#dfData = data.frame(lData.train$data)
dfData = data.frame(log(lData.train$data+1e-4))
dim(dfData)
dfData$fGroups = fGroups
str(dfData)
lData = list(resp=ifelse(dfData$fGroups == 'PA', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


# fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('tau', 'betas2'), init=initf, cores=4,
#                     control=list(adapt_delta=0.99, max_treedepth = 13))
# 
# save(fit.stan, file='temp/fit.stan.binom.serology.rds')

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
# ## get the intercept at population level
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]
colnames(mCoef) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]

## function to calculate statistics for a coefficient
getDifference = function(ivData){
  # get the difference vector
  d = ivData
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(p)
}

ivPval = apply(mCoef, 2, getDifference)
hist(ivPval)
plot(colMeans(mCoef), ivPval, pch=19)
m = colMeans(mCoef)
names(m) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
text(colMeans(mCoef), ivPval, names(m), pos=1)
m = abs(m)
m = sort(m, decreasing = T)

length(m)
#p.old = par(mar=c(6,3,4,2)+0.1)
l2 = barplot(m[1:14], 
             las=2, xaxt='n', col='grey', main='Top Variables', ylab='Absolute Log Odds')
axis(1, at = l2, labels = names(m)[1:14], tick = F, las=2, cex.axis=0.7 )

## format for line plots
m = colMeans(mCoef)
names(m) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
m = sort(m, decreasing = T)
mTreatment = mCoef[,names(m)]

df = apply(mTreatment, 2, getms)
x = 1:ncol(mTreatment)

par(p.old)
plot(x, df['m',], ylim=c(min(df), max(df)), pch=20, xlab='', main='Effect on Log Odds of PA',
     ylab='Slopes', xaxt='n')
axis(1, at = x, labels = colnames(mTreatment), las=2, cex.axis=0.7)
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(h = 0, col='grey')

### there appear to be a lot of correlations in the data
## find correlated variables
mCor = cor(log(mData+1e-4), use="na.or.complete")
library(caret)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.7, names=T)
data.frame(n)
sapply(n, function(x) {
  (abs(mCor[,x]) >= 0.7)
})

### drop some of the correlated variables
dim(mTreatment)
i = sample(1:2000, 1000)
mTreatment.sub = mTreatment[i,]
pairs(mTreatment.sub[,n], pch=20, col='grey')

# > data.frame(n)
# n
# 1               Ara_h_2_sIgE
# 2                Peanut_sIgE
# 3               Ara_h_1_sIgE
# 4    Peanut_specific_acivity
# 5               Ara_h_3_sIgE
# 6   Ara_h_1_specific_acivity
# 7 rAra_h_8_specific_activity

### use variable combinations
str(dfData)
set.seed(123)
oVar.sub = CVariableSelection.ReduceModel(dfData[,-15], dfData$fGroups, boot.num = 500)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)
# print variable combinations
for (i in 1:4){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

################# use the one variable at a time to model the data
#### write the functions to fit models
# 
# fit.0 = glm(fGroups ~ ., data=dfData, family='binomial')
# summary(fit.0)

library(LearnBayes)
## lets write a custom glm using a bayesian approach
## write the log posterior function
mylogpost = function(theta, data){
  ## parameters to track/estimate
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = plogis(iFitted)
  # write the priors and likelihood 
  lp = dnorm(betas[1], 0, 10, log=T) + sum(dnorm(betas[-1], 0, 10, log=T))
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lik + lp
  return(val)
}

## binomial prediction
mypred = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  mModMatrix = data$mModMatrix
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  #iFitted = plogis(iFitted)
  return(iFitted)
}

## set the variables
dfData.all = data.frame(log(lData.train$data+1e-4))
dim(dfData.all)

lFits = lapply(colnames(dfData.all), function(cPredictor) {
  dfData  = data.frame(dfData.all[,cPredictor], fGroups)
  colnames(dfData)[1] = cPredictor
  
  lData = list(resp=ifelse(dfData$fGroups == 'PS', 0, 1), 
               mModMatrix=model.matrix(fGroups ~ ., data=dfData))
  start = c(rep(0, times=ncol(lData$mModMatrix)))
  
  fit.1 = laplace(mylogpost, start, lData)
  ### lets take a sample from this 
  ## parameters for the multivariate t density
  tpar = list(m=fit.1$mode, var=fit.1$var*2, df=4)
  ## get a sample directly and using sir (sampling importance resampling with a t proposal density)
  s = sir(mylogpost, tpar, 5000, lData)
  colnames(s) = colnames(lData$mModMatrix)
  apply(s, 2, mean)
  apply(s, 2, sd)
  pairs(s, pch=20)
  fit.1$sir = s
  return(fit.1)
})

names(lFits) = colnames(dfData.all)


## plot of coefficients
temp = lapply(lFits, function(fit.1){
  s = fit.1$sir
  df = apply(s, 2, getms)
  x = 1:ncol(s)
  
  #par(p.old)
  plot(x, df['m',], ylim=c(min(df), max(df)), pch=20, xlab='', main='Effect on Log Odds of PA',
       ylab='Coefficients', xaxt='n')
  axis(1, at = x, labels = colnames(s), las=2, cex.axis=0.7)
  for(l in 1:ncol(df)){
    lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
  }
  abline(h = 0, col='grey')
})

## create the plots for regression with each predictor (not input) fixed at its average
## see Data Analysis ... Regression & Multilevel M [Gelman] for jitter.binary function
jitter.binary = function(a, jitt=.05){
  ifelse (a==0, runif (length(a), 0, jitt), runif (length(a), 1-jitt, 1))
}

temp = lapply(names(lFits), function(cPredictor){
  allergic.jitt = jitter.binary(ifelse(lData.train$covariates$Allergic.Status == 'PS', 0, 1))
  
  plot(dfData.all[,cPredictor], allergic.jitt, pch=20, xlab='Covariate', ylab='Probability of PA class',
       main=paste0('Prediction of PA class vs ', cPredictor))
  x = seq(min(dfData.all[,cPredictor]), max(dfData.all[,cPredictor]), length.out = 100)
  m = cbind(1, x)
  c = colMeans(lFits[[cPredictor]]$sir)
  lines(x, plogis(m %*% c), col='black')
})

### once we have results from the classifier we can make some plots to see
### the performance
library(lattice)
library(car)
library(ROCR)

lPredicted = lapply(names(lFits), function(cPredictor){
  ## create model matrix
  X = as.matrix(cbind(rep(1, times=nrow(dfData.all)), dfData.all[,cPredictor]))
  colnames(X) = colnames(lFits[[cPredictor]]$sir)
  ivPredict.raw = mypred(colMeans(lFits[[cPredictor]]$sir), list(mModMatrix=X))[,1]
  ivPredict = plogis(ivPredict.raw)
  xyplot(ivPredict ~ lData.train$covariates$Allergic.Status, xlab='Actual Group', ylab='Predicted Probability of Being PA (1)',
         main=paste0('Predicted scores vs Actual groups', cPredictor))
  
  ############# ROC curve 
  ## draw a ROC curve first for calibration performance test
  ivTruth = (lData.train$covariates$Allergic.Status == 'PA')
  p = prediction(ivPredict, ivTruth)
  perf.alive = performance(p, 'tpr', 'fpr')
  dfPerf.alive = data.frame(c=perf.alive@alpha.values, t=perf.alive@y.values[[1]], f=perf.alive@x.values[[1]], 
                            r=perf.alive@y.values[[1]]/perf.alive@x.values[[1]])
  colnames(dfPerf.alive) = c('c', 't', 'f', 'r')
  plot(perf.alive, main=paste0('Classifier Performance to predict PA ', cPredictor))
  return(ivPredict.raw)
})

names(lPredicted) = names(lFits)

