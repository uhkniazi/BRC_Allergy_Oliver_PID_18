# Name: 04_responseVsDiversity.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 3/12/2018
# Desc: diversity i.e. number of allergens recognised related to class prediction

############ data loading and cleaning
source('header.R')
library(lattice)
## load the data
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 36)')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
df = read.csv(n, header=T)

# close connection after getting data
dbDisconnect(db)

dim(df)
df = na.omit(df)
df$Patient = factor(gsub(' ', '', as.character(df$Patient)))
df$Allergic.Status = factor(gsub(' ', '', as.character(df$Allergic.Status)), levels = c('PS', 'PA'))
df = droplevels.data.frame(df)
## make count matrix
mData = as.matrix(df[,-c(1:2)])
dfSample = df[,1:2]
rownames(mData) = as.character(dfSample$Patient)
str(dfSample)
## remove NAs and 0s and convert other values to 1 by adding a jitter
f = mData < 0.3
table(f)
mData[f] = 0
mData[!f] = 1
dim(na.omit(mData))
dim(mData)

lData.train = list(data=mData, covariates=dfSample)
rm(df)

dim(mData)
# calculate diversity of allergens for each class
iClass.ps = colSums(mData[dfSample$Allergic.Status == 'PS', ])
iClass.pa = colSums(mData[dfSample$Allergic.Status == 'PA', ])

i = order(iClass.pa, decreasing = F)
m = cbind(pa=iClass.pa[i], ps=iClass.ps[i])
head(m)
df = stack(data.frame(m))
head(df)
df$names = factor(names(iClass.pa[i]), levels = names(iClass.pa[i]))
barchart(names ~ values | ind, data=df, xlab='Abundance', main='Ranked (on PA) Allergen Diversity', 
         scales=list(y=list(cex=0.6)))

## calculating diversity 
## diversity currently is defined as the number of allergens found in each sample
iData = rowSums(mData)
fGroups = lData.train$covariates$Allergic.Status

## load the stan data from binomial model from previous analysis
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load(file='temp/fit.stan.binom.binary.rds')

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

dfData = data.frame(lData.train$data)
dim(dfData)
dfData$fGroups = fGroups

lData = list(resp=ifelse(dfData$fGroups == 'PA', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
# ## get the intercept at population level
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]
colnames(mCoef) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]

m = colMeans(mCoef)
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

######## choose only the positive coefficients
m = colMeans(mCoef)
i = which(m >= 1)
m[i]

mData.pos = mData[,names(i)]
dim(mData.pos)
iData.pos = rowSums(mData.pos)

## choose the negative coefficients
m = colMeans(mCoef)
i = which(m <= -1)
m[i]

mData.neg = mData[,names(i)]
dim(mData.neg)
iData.neg = rowSums(mData.neg)

### data model
dfData = data.frame(allergic=fGroups, pos=iData.pos, neg=iData.neg)
str(dfData)

fit.0 = glm(allergic ~ pos + neg, data=dfData, family='binomial')
summary(fit.0)

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

lData = list(resp=ifelse(dfData$allergic == 'PS', 0, 1), 
             mModMatrix=model.matrix(allergic ~ pos + neg, data=dfData))
start = c(rep(0, times=ncol(lData$mModMatrix)))

mylogpost(start, lData)

fit.1 = laplace(mylogpost, start, lData)
fit.1
data.frame(coef(fit.0), fit.1$mode)
se = sqrt(diag(fit.1$var))

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

## plot of coefficients
df = apply(s, 2, getms)
x = 1:ncol(s)

par(p.old)
plot(x, df['m',], ylim=c(min(df), max(df)), pch=20, xlab='', main='Effect on Log Odds of PA',
     ylab='Coefficients', xaxt='n')
axis(1, at = x, labels = colnames(s), las=2, cex.axis=1)
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(h = 0, col='grey')

## create the plots for regression with each predictor (not input) fixed at its average
## see Data Analysis ... Regression & Multilevel M [Gelman] for jitter.binary function
jitter.binary = function(a, jitt=.05){
  ifelse (a==0, runif (length(a), 0, jitt), runif (length(a), 1-jitt, 1))
}

allergic.jitt = jitter.binary(lData$resp)

plot(dfData$pos, allergic.jitt, pch=20, xlab='Positive Effect', ylab='Probability of PA class',
     main='Prediction of PA class vs Abundance of Positive Allergens')
x = seq(min(dfData$pos), max(dfData$pos), length.out = 100)
m = cbind(1, x, mean(dfData$neg))
c = colMeans(fit.1$sir)
lines(x, plogis(m %*% c), col='black')
m = cbind(1, x, min(dfData$neg))
lines(x, plogis(m %*% c), col='red')
m = cbind(1, x, max(dfData$neg))
lines(x, plogis(m %*% c), col='green')
legend('bottomright', legend = c('Min Neg', 'Average Neg', 'Max Neg'), fill=c('red', 'black', 'green'))

## about 22% class PS
table(dfData$allergic[dfData$neg < 1])

## second predictor fixed i.e. pos
plot(dfData$neg, allergic.jitt, pch=20, xlab='Negative Effect', ylab='Probability of PA class',
     main='Prediction of PA class vs Abundance of Negative Allergens')
x = seq(min(dfData$neg), max(dfData$neg), length.out = 100)
m = cbind(1, mean(dfData$pos), x)
c = colMeans(fit.1$sir)
lines(x, plogis(m %*% c), col='black')
m = cbind(1, min(dfData$pos), x)
lines(x, plogis(m %*% c), col='red')
m = cbind(1, max(dfData$neg), x)
lines(x, plogis(m %*% c), col='green')
legend('right', legend = c('Min Pos', 'Average Pos', 'Max Pos'), fill=c('red', 'black', 'green'))

### once we have results from the classifier we can make some plots to see
### the performance
library(lattice)
library(car)
## get the predicted values
dfData.new = dfData
str(dfData.new)
## create model matrix
X = as.matrix(cbind(rep(1, times=nrow(dfData.new)), dfData.new[,colnames(fit.1$sir)[-1]]))
colnames(X) = colnames(fit.1$sir)
head(X)
ivPredict.raw = mypred(colMeans(fit.1$sir), list(mModMatrix=X))[,1]
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

# convert to logit scale for model fitting
ivPredict = ivPredict.raw

################################ section for mixture model
######## this mixture model will help us decide an appropriate cutoff for the decision rule
######## see Gelman 2013 around P18 for an example of record linking score calibration
stanDso = rstan::stan_model(file='normResponseFiniteMixture_2.stan')

## take a subset of the data
lStanData = list(Ntotal=length(ivPredict), y=ivPredict, iMixtures=2)

## give initial values if you want, look at the density plot 
initf = function(chain_id = 1) {
  list(mu = c(-5, 7), sigma = c(1, 1), iMixWeights=c(0.5, 0.5))
} 

## give initial values function to stan
# l = lapply(1, initf)
fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, cores=4, init=initf)
print(fit.stan, digi=3)
traceplot(fit.stan)
save(fit.stan, file='temp/fit.stan.mixture_diversity.rds')
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
ivPredict = plogis(ivPredict.raw)

## draw a ROC curve first for calibration performance test
ivTruth = fGroups == 'PA'
p = prediction(ivPredict, ivTruth)
perf.alive = performance(p, 'tpr', 'fpr')
dfPerf.alive = data.frame(c=perf.alive@alpha.values, t=perf.alive@y.values[[1]], f=perf.alive@x.values[[1]], 
                          r=perf.alive@y.values[[1]]/perf.alive@x.values[[1]])
colnames(dfPerf.alive) = c('c', 't', 'f', 'r')
plot(perf.alive, main='Classifier Performance to predict PA')

## draw the simulation lines
## these are p-values from the mixture components
## create posterior smatter lines
grid = seq(-16, 14, length.out = 100)
f_getSmatterLines = function(m, s, g){
  return(pnorm(g, m, s, lower.tail = F))
}
y = f_getSmatterLines(5.58, 3.98, grid)
x = f_getSmatterLines(-4.43, 4.21, grid)
lines(x, y, col=2, lwd=2)

## holders for the simulated p-values
mTP = matrix(NA, nrow = length(grid), ncol = 200)
mFP = matrix(NA, nrow = length(grid), ncol = 200)

for (i in 1:200){
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
x = pnorm(logit(0.5), mStan[p, 'mu1'], mStan[p, 'sigma1'], lower.tail = F)
y = pnorm(logit(0.5), mStan[p, 'mu2'], mStan[p, 'sigma2'], lower.tail=F)

hist(x, main='False Positive Rate at 0.5', xlab='')
hist(y, main='True Positive Rate at 0.5', xlab='')

fPredict = rep('PS', times=length(ivPredict))
fPredict[ivPredict >= 0.5] = 'PA'
table(fPredict, fGroups)
# TP = 46/(46+4)
# FP = 3/(3+37)
