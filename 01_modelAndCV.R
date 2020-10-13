# Name: 01_modelAndCV.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 13/10/2020
# Desc: model fitting with cross validation and roc curves

############ data loading and cleaning
source('header.R')

## load the training set
df = read.csv(file.choose(), header=T)
dim(df)

# remove NA patients
f = is.na(df$Allergic.Status)
table(f)
df = df[!f,]
dim(df)
# remove any additi
df = na.omit(df)
dim(df)
## remove white space
colnames(df)
df$Sample = factor(gsub(' ', '', as.character(df$Sample)))
df$Allergic.Status = factor(gsub(' ', '', as.character(df$Allergic.Status)), levels = c('PS', 'PA'))
df = droplevels.data.frame(df)

## make count matrix
mData = as.matrix(df[,-c(1:2)])
dfSample = df[,1:2]
rownames(mData) = as.character(dfSample$Sample)
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

#### load the test/validation data
df = read.csv(file.choose(), header=T)
dim(df)

# remove NA patients
f = is.na(df$Allergic.Status)
table(f)
df = df[!f,]
dim(df)
# remove any additi
df = na.omit(df)
dim(df)
## remove white space
colnames(df)
df$Sample = factor(gsub(' ', '', as.character(df$Sample)))
df$Allergic.Status = factor(gsub(' ', '', as.character(df$Allergic.Status)), levels = c('PS', 'PA'))
df = droplevels.data.frame(df)

## make count matrix
mData = as.matrix(df[,-c(1:2)])
dfSample = df[,1:2]
rownames(mData) = as.character(dfSample$Sample)
str(dfSample)

## remove NAs and 0s and convert other values to 1 by adding a jitter
f = mData < 0.3
table(f)
mData[f] = 0
mData[!f] = 1
dim(na.omit(mData))
dim(mData)

lData.test = list(data=mData, covariates=dfSample)
rm(df)

## sanity check
table(colnames(lData.train$data) %in% colnames(lData.test$data))
############ end data loading

######################## Stan section for binomial regression approach
### fit 2 models one standard binomial and one with robust version
dfData = data.frame(lData.train$data)
dim(dfData)
dfData$fGroups = lData.train$covariates$Allergic.Status
lData = list(resp=ifelse(dfData$fGroups == 'PA', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

library(rethinking)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

# ## give initial values
# initf = function(chain_id = 1) {
#   list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
# }


fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=4, pars=c('tau', 'betas2', 'log_lik'), cores=4,# init=initf,
                    control=list(adapt_delta=0.99, max_treedepth = 13))

#save(fit.stan, file='temp/fit.stan.binom.binary.rds')

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

## line plot of coefficients
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

## export the coefficients to results file
m = c(mean(iIntercept), colMeans(mTreatment))
s = c(sd(iIntercept), apply(mTreatment, 2, sd))
r = signif(cbind(m, s), 3)
colnames(r) = c('Coefficient', 'SE')
rownames(r)[1] = 'Intercept'

write.csv(r, file = 'results/model_1_robust_train.csv')

################# predictions and comparison with robust model
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

mCoef = extract(fit.stan)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)])
library(lattice)
## get the predicted values
## create model matrix
X = as.matrix(cbind(rep(1, times=nrow(dfData)), dfData[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
xyplot(ivPredict ~ fGroups, xlab='Actual Group', 
       ylab='Predicted Probability of Being PA (1) - Train',
       data=dfData)
# outlier samples
i = which(ivPredict < 0.5 & dfData$fGroups == 'PA')
cvOutliers = names(i)

## repeat on test data
## create model matrix
dfData = data.frame(lData.test$data)
dim(dfData)
dfData$fGroups = lData.test$covariates$Allergic.Status 
X = as.matrix(cbind(rep(1, times=nrow(dfData)), dfData[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
xyplot(ivPredict ~ fGroups, xlab='Actual Group', 
       ylab='Predicted Probability of Being PA (1) - Test',
       data=dfData)
# outlier samples
i = which(ivPredict < 0.5 & dfData$fGroups == 'PA')
cvOutliers = names(i)

#### compare the 2 models
m.s = fit.stan
m.r = fit.stan
plot(compare(m.s, m.r))
# # remove appropriate covariate or samples
# m = lData.train$data
# colnames(m)
# m = m[,-6]
# lData.train$data = m

## remove outliers
# i = which(lData.train$covariates$Patient %in% cvOutliers)
# lData.train$data = lData.train$data[-i,]
# lData.train$covariates = lData.train$covariates[-i,]
## plots of coeftab, from various models
## to show that the guessing parameter controls for outliers
ct = coeftab(m.s, m.r)
rn = rownames(ct@coefs)
i = grep('betas', rn)
i = i[-c(1, length(i))]
rownames(ct@coefs)[i] = colnames(mData)
rownames(ct@se)[i] = colnames(mData)
plot(ct, pars=colnames(mData))

#################

dim(dfData)
# create the cross validation object
url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/bernoulli.stan'
download(url, 'bernoulli.stan')

oCV.s = CCrossValidation.StanBern(dfData[,-113], dfData[, -113], fGroups, fGroups, level.predict = 'PA',
                                  boot.num = 10, k.fold = 10, ncores = 2, nchains = 2) 

save(oCV.s, file='temp/oCV.s.rds')

plot.cv.performance(oCV.s)
unlink('bernoulli.stan')
################################################# binomial regression with mixture model section
################ fit a binomial model on the chosen model size based on previous results
## this can be another classifier as well e.g. LDA. Using this model check how is the performance 
## and using this make some calibration curves to select decision boundary


library(LearnBayes)
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


lData = list(resp=ifelse(dfData$fGroups == 'PA', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))
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

## identify possible outliers/misclassified observations
df = data.frame(fGroups, ivPredict)
i = which(df$fGroups == 'PS' & df$ivPredict > 0.4)
rownames(df)[i]
i = which(df$fGroups == 'PA' & df$ivPredict < 0.5)
rownames(df)[i]
## lets check on a different scale of the score
densityplot(~ ivPredict.raw, data=dfData)
xyplot(ivPredict.raw ~ lData.train$covariates$Allergic.Status, xlab='Actual Group', ylab='Predicted Probability of Being PS (1)')
densityplot(~ ivPredict.raw, groups=fGroups, data=dfData, type='n', 
            xlab='Predicted Score', main='Logit Scale', auto.key = list(columns=2))


############# ROC curve 
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

# convert to logit scale for model fitting
ivPredict = ivPredict.raw
################################ section for mixture model
######## this mixture model will help us decide an appropriate cutoff for the decision rule
######## see Gelman 2013 around P18 for an example of record linking score calibration
fit.stan.bin = fit.stan

stanDso = rstan::stan_model(file='normResponseFiniteMixture_2.stan')

## take a subset of the data
lStanData = list(Ntotal=length(ivPredict), y=ivPredict, iMixtures=2)

## give initial values if you want, look at the density plot 
initf = function(chain_id = 1) {
  list(mu = c(-5, 5), sigma = c(1, 1), iMixWeights=c(0.5, 0.5))
} 

## give initial values function to stan
# l = lapply(1, initf)
fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, cores=4, init=initf)
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
grid = seq(-7, 7, length.out = 100)
f_getSmatterLines = function(m, s, g){
  return(pnorm(g, m, s, lower.tail = F))
}
y = f_getSmatterLines(4.66, 3, grid)
x = f_getSmatterLines(-4.01, 1.98, grid)
lines(x, y, col=2, lwd=2)

## holders for the simulated p-values
mTP = matrix(NA, nrow = length(grid), ncol = 2000)
mFP = matrix(NA, nrow = length(grid), ncol = 2000)

for (i in 1:2000){
  p = i #sample(1:nrow(mStan), size = 1)
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
