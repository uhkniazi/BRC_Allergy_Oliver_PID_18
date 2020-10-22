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
## subset of inputs to use
cvSubset = c('Ara.h.1', 'Ara.h.2', 'Ara.h.3', 'Ara.h.6')
df$Sample = factor(gsub(' ', '', as.character(df$Sample)))
df$Allergic.Status = factor(gsub(' ', '', as.character(df$Allergic.Status)), levels = c('PS', 'PA'))
df = droplevels.data.frame(df)

## make count matrix
mData = as.matrix(df[,-c(1:2)])
dfSample = df[,1:2]
rownames(mData) = as.character(dfSample$Sample)
str(dfSample)

mData = mData[,cvSubset]
## remove NAs and 0s and convert other values to 1 by adding a jitter
# f = mData < 0.3
# table(f)
# mData[f] = 0
# mData[!f] = 1
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

mData = mData[,colnames(lData.train$data)]
## remove NAs and 0s and convert other values to 1 by adding a jitter
# f = mData < 0.3
# table(f)
# mData[f] = 0
# mData[!f] = 1
dim(na.omit(mData))
dim(mData)

lData.test = list(data=mData, covariates=dfSample)
rm(df)
rm(mData)
## sanity check
table(colnames(lData.train$data) %in% colnames(lData.test$data))
identical(colnames(lData.train$data), colnames(lData.test$data))
############ end data loading

######################## Stan section for binomial regression approach
### fit 2 models one standard binomial and one with robust version
dfData = rbind(data.frame(lData.train$data), data.frame(lData.test$data))
dim(dfData)
dfData$fGroups = factor(c(as.character(lData.train$covariates$Allergic.Status), 
                          as.character(lData.test$covariates$Allergic.Status)), levels = c('PS', 'PA'))
lData = list(resp=ifelse(dfData$fGroups == 'PA', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

library(rethinking)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialGuessMixtureRegressionSharedCoeffVariance.stan')

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

#### compare the 2 models
m.s = fit.stan
m.r = fit.stan
plot(compare(m.s, m.r))
compare(m.s, m.r)

# choose the appropriate model to work with
fit.stan = m.r

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

write.csv(r, file = 'results/model_5_subset_n100.csv')

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
       ylab='Predicted Probability of Being PA (1) - All',
       data=dfData)
# outlier samples
i = which(ivPredict < 0.5 & dfData$fGroups == 'PA')
cvOutliers = names(i)
cvOutliers

# ## repeat on test data
# ## create model matrix
# dfData = data.frame(lData.test$data)
# dim(dfData)
# dfData$fGroups = lData.test$covariates$Allergic.Status 
# X = as.matrix(cbind(rep(1, times=nrow(dfData)), dfData[,colnames(mCoef)[-1]]))
# colnames(X) = colnames(mCoef)
# head(X)
# ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
# xyplot(ivPredict ~ fGroups, xlab='Actual Group', 
#        ylab='Predicted Probability of Being PA (1) - Test',
#        data=dfData)
# # outlier samples
# i = which(ivPredict < 0.5 & dfData$fGroups == 'PA')
# cvOutliers = names(i)
# 
# #### compare the 2 models
# m.s = fit.stan
# m.r = fit.stan
# plot(compare(m.s, m.r))

rm(dfData)
rm(stanDso)

#### repeat the analysis on the new dataset only separately
dfData = data.frame(lData.test$data)
dim(dfData)
dfData$fGroups = lData.test$covariates$Allergic.Status
lData = list(resp=ifelse(dfData$fGroups == 'PA', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))
## scroll up to the model section to create figures and output file

#################
####################### cross validation under LDA and binomial models
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

## setup input data
dfData.train = data.frame(lData.train$data)
fGroups.train = lData.train$covariates$Allergic.Status

dfData.test = data.frame(lData.test$data)
fGroups.test = lData.test$covariates$Allergic.Status

dim(dfData.train)
dim(dfData.test)

# # create the cross validation object
# url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/bernoulli.stan'
# download(url, 'bernoulli.stan')

oCV.s = CCrossValidation.StanBern(train.dat = dfData.train, 
                                  test.dat = dfData.test, 
                                  test.groups = fGroups.test, 
                                  train.groups = fGroups.train,
                                  level.predict = 'PA',
                                  boot.num = 10, k.fold = 10, 
                                  ncores = 2, nchains = 2) 

save(oCV.s, file='temp/oCV.s_m5_subset.rds')

plot.cv.performance(oCV.s)
unlink('bernoulli.stan')
