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


fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('tau', 'betas2'), init=initf, cores=4,
                    control=list(adapt_delta=0.99, max_treedepth = 13))

save(fit.stan, file='temp/fit.stan.binom.serology.rds')

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

# ## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)])
# pairs(mCoef, pch=20)


### once we have results from the classifier we can make some plots to see
### the performance
library(lattice)
library(car)
## get the predicted values
dfData.new = dfData
## create model matrix
X = as.matrix(cbind(rep(1, times=nrow(dfData.new)), dfData.new[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict.raw = mypred(colMeans(mCoef), list(mModMatrix=X))[,1]
ivPredict = plogis(ivPredict.raw)
xyplot(ivPredict ~ fGroups, xlab='Actual Group', ylab='Predicted Probability of Being PA (1)')
xyplot(ivPredict ~ lData.train$covariates$Allergic.Status, xlab='Actual Group', ylab='Predicted Probability of Being PA (1)',
       main='Predicted scores vs Actual groups')
densityplot(~ ivPredict, data=dfData, type='n')
densityplot(~ ivPredict | fGroups, data=dfData, type='n', xlab='Predicted Score', main='Actual Scale')
densityplot(~ ivPredict, groups=fGroups, data=dfData, type='n', 
            xlab='Predicted Score', main='Actual Scale', auto.key = list(columns=2))

## lets check on a different scale of the score
densityplot(~ ivPredict.raw, data=dfData)
xyplot(ivPredict.raw ~ lData.train$covariates$Allergic.Status, xlab='Actual Group', ylab='Predicted Probability of Being PS (1)')
densityplot(~ ivPredict.raw, groups=fGroups, data=dfData, type='n', 
            xlab='Predicted Score', main='Logit Scale', auto.key = list(columns=2))


############# ROC curve 
## draw a ROC curve first for calibration performance test
library(ROCR)
ivTruth = fGroups == 'PA'
p = prediction(ivPredict, ivTruth)
perf.alive = performance(p, 'tpr', 'fpr')
dfPerf.alive = data.frame(c=perf.alive@alpha.values, t=perf.alive@y.values[[1]], f=perf.alive@x.values[[1]], 
                          r=perf.alive@y.values[[1]]/perf.alive@x.values[[1]])
colnames(dfPerf.alive) = c('c', 't', 'f', 'r')
plot(perf.alive, main='Classifier Performance to predict PA')
a = performance(p, 'auc')
legend('bottomright', paste0('auc ', as.numeric(a@y.values)))

############################################################################
########### univariate models
############################################################################

# ### there appear to be a lot of correlations in the data
# ## find correlated variables
# mCor = cor(log(mData+1e-4), use="na.or.complete")
# library(caret)
# ### find the columns that are correlated and should be removed
# n = findCorrelation((mCor), cutoff = 0.7, names=T)
# data.frame(n)
# sapply(n, function(x) {
#   (abs(mCor[,x]) >= 0.7)
# })
# 
# ### drop some of the correlated variables
# dim(mTreatment)
# i = sample(1:2000, 1000)
# mTreatment.sub = mTreatment[i,]
# pairs(mTreatment.sub[,n], pch=20, col='grey')
# 
# # > data.frame(n)
# # n
# # 1               Ara_h_2_sIgE
# # 2                Peanut_sIgE
# # 3               Ara_h_1_sIgE
# # 4    Peanut_specific_acivity
# # 5               Ara_h_3_sIgE
# # 6   Ara_h_1_specific_acivity
# # 7 rAra_h_8_specific_activity

### plot the data
df = dfData[,-15]
df = stack(df)
df$fGroups = fGroups
bwplot(values ~ fGroups | ind, data=df)
rm(df)

### use variable combinations
str(dfData)
set.seed(123)
oVar.sub = CVariableSelection.ReduceModel(dfData[,-15], dfData$fGroups, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)
# print variable combinations
for (i in 1:4){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

################# use the one variable at a time to model the data
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
  #pairs(s, pch=20)
  fit.1$sir = s
  return(fit.1)
})

names(lFits) = colnames(dfData.all)

par(mfrow=c(2,2))
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
       main=paste0('Variable PA ', cPredictor))
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
par(mfrow=c(2,2))
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
  plot(perf.alive, main=paste0('Performance to predict PA ', cPredictor))
  a = performance(p, 'auc')
  legend('bottomright', paste0('auc ', as.numeric(a@y.values)))
  return(ivPredict.raw)
})

names(lPredicted) = names(lFits)

## further analysis for top variable of choice
ivPredict = lPredicted$Ara_h_2_sp_ac
densityplot(~ ivPredict, groups=fGroups)
xyplot(plogis(ivPredict) ~ lData.train$covariates$Allergic.Status, xlab='Actual Group', ylab='Predicted Probability of Being PA (1)',
       main='Predicted scores vs Actual groups')
################################ section for mixture model
######## this mixture model will help us decide an appropriate cutoff for the decision rule
######## see Gelman 2013 around P18 for an example of record linking score calibration
stanDso = rstan::stan_model(file='normResponseFiniteMixture_2.stan')

## take a subset of the data
lStanData = list(Ntotal=length(ivPredict), y=ivPredict, iMixtures=2)

## give initial values if you want, look at the density plot 
initf = function(chain_id = 1) {
  list(mu = c(-5, 5), sigma = c(1, 1), iMixWeights=c(0.5, 0.5))
} 

## give initial values function to stan
# l = lapply(1, initf)
fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, cores=2, init=initf)
print(fit.stan, digi=3)
traceplot(fit.stan)
save(fit.stan, file='temp/fit.stan.mixture.rds')
## check if labelling degeneracy has occured
## see here: http://mc-stan.org/users/documentation/case-studies/identifying_mixture_models.html
params1 = as.data.frame(extract(fit.stan, permuted=FALSE)[,1,])
params2 = as.data.frame(extract(fit.stan, permuted=FALSE)[,2,])
params3 = as.data.frame(extract(fit.stan, permuted=FALSE)[,3,])
params4 = as.data.frame(extract(fit.stan, permuted=FALSE)[,4,])

## check if the means from different chains overlap
## Labeling Degeneracy by Enforcing an Ordering
par(mfrow=c(2,2))
plot(params1$`mu[1]`, params1$`mu[2]`, pch=20, col=2)
plot(params2$`mu[1]`, params2$`mu[2]`, pch=20, col=3)
plot(params3$`mu[1]`, params3$`mu[2]`, pch=20, col=4)
plot(params4$`mu[1]`, params4$`mu[2]`, pch=20, col=5)

par(mfrow=c(1,1))
plot(params1$`mu[1]`, params1$`mu[2]`, pch=20, col=2)
points(params2$`mu[1]`, params2$`mu[2]`, pch=20, col=3)
points(params3$`mu[1]`, params3$`mu[2]`, pch=20, col=4)
points(params4$`mu[1]`, params4$`mu[2]`, pch=20, col=5)


# model checks
############# extract the mcmc sample values from stan
mStan = do.call(cbind, extract(fit.stan))
mStan = mStan[,-(ncol(mStan))]
colnames(mStan) = c('mu1', 'mu2', 'sigma1', 'sigma2', 'mix1', 'mix2')
dim(mStan)
## get a sample for this distribution
########## simulate 200 test quantities
mDraws = matrix(NA, nrow = length(ivPredict), ncol=200)

for (i in 1:200){
  p = sample(1:nrow(mStan), size = 1)
  mix = mean(mStan[,'mix1'])
  ## this will take a sample from a normal mixture distribution
  sam = function() {
    ind = rbinom(1, 1, prob = mix)
    return(ind * rnorm(1, mStan[p, 'mu1'], mStan[p, 'sigma1']) + 
             (1-ind) * rnorm(1, mStan[p, 'mu2'], mStan[p, 'sigma2']))
  }
  mDraws[,i] = replicate(length(ivPredict), sam())
}

mDraws.normMix = mDraws

yresp = density(ivPredict)
yresp$y = yresp$y/max(yresp$y)
plot(yresp, xlab='', main='Fitted distribution', ylab='scaled density', lwd=2)
temp = apply(mDraws, 2, function(x) {x = density(x)
x$y = x$y/max(x$y)
lines(x, col='darkgrey', lwd=0.6)
})
lines(yresp, lwd=2)


print(fit.stan)

range(ivPredict)
## reconvert back to inverse logit scale i.e. 0 to 1 range
ivPredict = plogis(lPredicted$Ara_h_2_sp_ac)

## draw a ROC curve first for calibration performance test
ivTruth = fGroups == 'PA'
p = prediction(ivPredict, ivTruth)
perf.alive = performance(p, 'tpr', 'fpr')
dfPerf.alive = data.frame(c=perf.alive@alpha.values, t=perf.alive@y.values[[1]], f=perf.alive@x.values[[1]], 
                          r=perf.alive@y.values[[1]]/perf.alive@x.values[[1]])
colnames(dfPerf.alive) = c('c', 't', 'f', 'r')
plot(perf.alive, main='Classifier Performance to predict PA')
a = performance(p, 'auc')
legend('bottomright', paste0('auc ', as.numeric(a@y.values)))

## draw the simulation lines
## these are p-values from the mixture components
## create posterior smatter lines
grid = seq(-8, 7, length.out = 100)
f_getSmatterLines = function(m, s, g){
  return(pnorm(g, m, s, lower.tail = F))
}
y = f_getSmatterLines(4.15, 1.61, grid)
x = f_getSmatterLines(-2.07, 2.21, grid)
lines(x, y, col=2, lwd=2)

## holders for the simulated p-values
mTP = matrix(NA, nrow = length(grid), ncol = 2000)
mFP = matrix(NA, nrow = length(grid), ncol = 2000)

for (i in 1:2000){
  p = sample(1:nrow(mStan), size = 1)
  x = pnorm(grid, mStan[p, 'mu1'], mStan[p, 'sigma1'], lower.tail = F) 
  y = pnorm(grid, mStan[p, 'mu2'], mStan[p, 'sigma2'], lower.tail=F)
  lines(x, y, col='darkgrey', lwd=0.5)
  mFP[,i] = x
  mTP[,i] = y
}

plot(perf.alive, add=T, col='blue', lwd=2)

c = cbind(tp=rowMeans(mTP), fp=rowMeans(mFP))
matplot(c, type = 'l', xaxt='n', xlab='Decision Boundary', ylab='Average Rate',
        main='Simulated True Positive & False Positive Rates')
legend('topright', c('TP', 'FP'), fill=c('black', 'red'))
axis(1, 1:nrow(c), labels = round(plogis(grid), 3), cex.axis=0.6, las=2)

## calculate average scores via simulation at desired cutoff
p = sample(1:nrow(mStan), size = 2000)
x = pnorm(logit(0.57), mStan[p, 'mu1'], mStan[p, 'sigma1'], lower.tail = F)
y = pnorm(logit(0.57), mStan[p, 'mu2'], mStan[p, 'sigma2'], lower.tail=F)

hist(x, main='False Positive Rate at 0.57', xlab='')
hist(y, main='True Positive Rate at 0.57', xlab='')

fPredict = rep('PS', times=length(ivPredict))
fPredict[ivPredict >= 0.57] = 'PA'
table(fPredict, fGroups)
