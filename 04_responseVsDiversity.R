# Name: 04_responseVsDiversity.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 3/12/2018
# Desc: diversity i.e. number of allergens recognised related to class prediction

############ data loading and cleaning
source('header.R')

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

## create the plots for regression with each predictor (not input) fixed at its average
## see Data Analysis ... Regression & Multilevel M [Gelman] for jitter.binary function
jitter.binary = function(a, jitt=.05){
  ifelse (a==0, runif (length(a), 0, jitt), runif (length(a), 1-jitt, 1))
}

allergic.jitt = jitter.binary(lData$resp)

plot(dfData$pos, allergic.jitt, pch=20)
x = seq(min(dfData$pos), max(dfData$pos), length.out = 100)
m = cbind(1, x, mean(dfData$neg))
c = colMeans(fit.1$sir)
lines(x, plogis(m %*% c))
m = cbind(1, x, min(dfData$neg))
lines(x, plogis(m %*% c))
m = cbind(1, x, max(dfData$neg))
lines(x, plogis(m %*% c))
