# File: 01_EDA.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: initial data and covariate checks
# Date: 12/06/2018

source('header.R')

setwd('dataID_40/')

## load the data
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 40)')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
df = read.csv(n[1], header=T)

# close connection after getting data
dbDisconnect(db)

str(df)
df$Patient.ID = factor(gsub(' ', '', as.character(df$Patient.ID)))
df$Allergic.Status = gsub(' ', '', as.character(df$Allergic.Status))
df$Allergic.Status[is.na(df$Allergic.Status)] = 'Healthy'
df$Allergic.Status = factor(df$Allergic.Status, levels = c('Healthy', 'PS', 'PA'))
df = droplevels.data.frame(df)

dfData = df
## check where is imputation required
apply(df, 2, function(x) sum(is.na(x)))
summary(dfData)

library(lattice)
df = dfData[,-c(1,2)]
summary(df)
df = stack(df)
df$treatment = dfData$Allergic.Status

densityplot(~ values | ind, groups=treatment, data=df[df$treatment != 'Healthy',], scales=list(relation='free'), type='n',
            auto.key=list(columns=3))


densityplot(~ values | treatment, data=df[df$ind == 'CD63.Act',], scales=list(relation='free'), type='n',
            auto.key=list(columns=3))

df = dfData[dfData$Allergic.Status == 'PA',-c(1,2,16,19)]
pairs(log(df), pch=20)

dfData = dfData[dfData$Allergic.Status != 'Healthy', ]
dfData = droplevels.data.frame(dfData)
dfData = dfData[,-c(16)]
summary(dfData)

## impute the average values for missing data 
f = is.na(dfData$Peanut.SPT)
table(f)
tapply(dfData$Peanut.SPT, dfData$Allergic.Status, mean, na.rm=T)
f = which(is.na(dfData$Peanut.SPT))
dfData[f,c('Allergic.Status', 'Peanut.SPT')]
dfData[c(25, 88), 'Peanut.SPT'] = 3.55
dfData[c(90), 'Peanut.SPT'] = 9.69

f = is.na(dfData$Age)
table(f)
tapply(dfData$Age, dfData$Allergic.Status, mean, na.rm=T)
f = which(is.na(dfData$Age))
dfData[f,c('Allergic.Status', 'Age')]
dfData[c(88), 'Age'] = 5.88
dfData[c(90), 'Age'] = 8.87

summary(dfData)

## save this data for further analysis
n = make.names(paste('Imputed data for mechanistic analysis data ID 40 rds'))
n2 = paste0('~/Data/MetaData/', n)
#save(dfData, file=n2)

## note: comment out as this entry has been made in db
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=40, name=n, type='rds', location='~/Data/MetaData/', comment='Olivers allergen data after imputation for PS and PA groups')
#dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)

################# preliminary model fitting
dfData.bk = dfData
dfData = dfData[dfData$Allergic.Status == 'PA',]
dfData = dfData[,-c(1,2)]

df = stack(dfData[,-c(16)])
df$y = dfData$CD63.Act

xyplot(y ~ log(values+0.5) | ind, data=df, type=c('g', 'p', 'r'), pch=19, cex=0.6,
       index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'))

xyplot(y ~ values | ind, data=df, type=c('g', 'p', 'smooth'), pch=19, cex=0.6,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'))

xyplot(y ~ log(values+0.5) | ind, data=df, type=c('g', 'p', 'smooth'), pch=19, cex=0.6,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'))


xyplot(y ~ log(values+0.5) | ind, data=df, type=c('smooth'), pch=19, cex=0.6,
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'))

xyplot(y ~ values | ind, data=df, type=c('smooth'), pch=19, cex=0.6,
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'))

##### transformation of the covariates
tapply(df$values, df$ind, function(x) quantile(x, 0:10/10))
data.frame(colnames(dfData))
dfData = dfData[,-c(9, 10)]
# transform certain values, before log transformation
dfData$Peanut.Sp.Act[dfData$Peanut.Sp.Act < 1] = 1e-4 
dfData$Peanut.Sp.Act = log(dfData$Peanut.Sp.Act)
dfData$Peanut.Sp.Act[dfData$Peanut.Sp.Act < 0] = 0

dfData$Ara.h.2.Sp.Act[dfData$Ara.h.2.Sp.Act < 1] = 1e-4 
dfData$Ara.h.2.Sp.Act = log(dfData$Ara.h.2.Sp.Act)
dfData$Ara.h.2.Sp.Act[dfData$Ara.h.2.Sp.Act < 0] = 0 

dfData$Ara.h.6.Sp.Act[dfData$Ara.h.6.Sp.Act < 1] = 1e-4 
dfData$Ara.h.6.Sp.Act = log(dfData$Ara.h.6.Sp.Act)
dfData$Ara.h.6.Sp.Act[dfData$Ara.h.6.Sp.Act < 0] = 0

dfData$f423.Ara.h.2[dfData$f423.Ara.h.2 < 1] = 1e-4
dfData$f423.Ara.h.2 = log(dfData$f423.Ara.h.2)
dfData$f423.Ara.h.2[dfData$f423.Ara.h.2 < 0] = 0

dfData$f424.rAra.h.3[dfData$f424.rAra.h.3 < 1] = 1e-4
dfData$f424.rAra.h.3 = log(dfData$f424.rAra.h.3)
dfData$f424.rAra.h.3[dfData$f424.rAra.h.3 < 0] = 0

dfData$f423.nAra.H.6[dfData$f423.nAra.H.6 < 1] = 1e-4
dfData$f423.nAra.H.6 = log(dfData$f423.nAra.H.6)
dfData$f423.nAra.H.6[dfData$f423.nAra.H.6 < 0] = 0

dfData$f13.Peanut[dfData$f13.Peanut < 1] = 1e-4
dfData$f13.Peanut = log(dfData$f13.Peanut)
dfData$f13.Peanut[dfData$f13.Peanut < 0] = 0

dfData$f422.rAra.h.1[dfData$f422.rAra.h.1 < 1] = 1e-4
dfData$f422.rAra.h.1 = log(dfData$f422.rAra.h.1)
dfData$f422.rAra.h.1[dfData$f422.rAra.h.1 < 0] = 0

dfData$total.IgE = log(dfData$total.IgE)

df = stack(dfData[,-c(14)])
df$y = dfData$CD63.Act

xyplot(y ~ values | ind, data=df, type=c('g', 'p', 'smooth'), pch=19, cex=0.6,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'))

xyplot(y ~ values | ind, data=df, type=c('smooth'), pch=19, cex=0.6,
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'))

# fit.1 = lm(CD63.Act ~ Peanut.Sp.Act, data=dfData)
# fit.2 = lm(CD63.Act ~ log(Peanut.Sp.Act), data=dfData)
# summary(fit.1)
# summary(fit.2)
# 
# s1 = simulate(fit.1, 20)
# s2 = simulate(fit.2, 20)
# par(mfrow=c(2,3))
# plot(density(dfData$CD63.Act))
# 
# plot(density(simulate(fit.1, 1)[,1]))
# apply(s1, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
# lines(density(dfData$CD63.Act), col=2)
# 
# plot(density(simulate(fit.2, 1)[,1]))
# apply(s2, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
# lines(density(dfData$CD63.Act), col=2)
# 
# dfData.new = dfData
# dfData.new$Peanut.Sp.Act = log(dfData$Peanut.Sp.Act)
# 
# ## repeat for each covariate
# fit.1 = lm(CD63.Act ~ Peanut.SPT, data=dfData)
# fit.2 = lm(CD63.Act ~ log(Peanut.SPT), data=dfData)
# summary(fit.1)
# summary(fit.2)
# 
# s1 = simulate(fit.1, 20)
# s2 = simulate(fit.2, 20)
# par(mfrow=c(2,3))
# plot(density(dfData$CD63.Act))
# 
# plot(density(simulate(fit.1, 1)[,1]))
# apply(s1, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
# lines(density(dfData$CD63.Act), col=2)
# 
# plot(density(simulate(fit.2, 1)[,1]))
# apply(s2, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
# lines(density(dfData$CD63.Act), col=2)
# 
# ###
# fit.1 = lm(CD63.Act ~ total.IgE, data=dfData)
# fit.2 = lm(CD63.Act ~ log(total.IgE), data=dfData)
# summary(fit.1)
# summary(fit.2)
# 
# s1 = simulate(fit.1, 20)
# s2 = simulate(fit.2, 20)
# par(mfrow=c(2,3))
# plot(density(dfData$CD63.Act))
# 
# plot(density(simulate(fit.1, 1)[,1]))
# apply(s1, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
# lines(density(dfData$CD63.Act), col=2)
# 
# plot(density(simulate(fit.2, 1)[,1]))
# apply(s2, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
# lines(density(dfData$CD63.Act), col=2)
# 
# dfData.new$total.IgE = log(dfData$total.IgE)
# 

###
fit.1 = lm(CD63.Act ~ ., data=dfData)
summary(fit.1)
s1 = simulate(fit.1, 20)

par(mfrow=c(2,2))
plot(density(dfData$CD63.Act))

plot(density(simulate(fit.1, 1)[,1]))
apply(s1, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
lines(density(dfData$CD63.Act), col=2)

#################### use stan to generate MCMC sample
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='tResponseRegression.stan')

m = model.matrix(CD63.Act ~ ., data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 y=dfData$CD63.Act)

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'mu', 'sigmaPop', 'nu'),
                    cores=2)
print(fit.stan, c('betas', 'sigmaPop', 'nu'), digits=3)

# some diagnostics for stan
traceplot(fit.stan, c('sigmaPop', 'nu'), ncol=1, inc_warmup=F)
pairs(fit.stan, pars = c("sigmaPop", "nu", "lp__"))
pairs(fit.stan, pars = c("betas", "lp__"))

m = extract(fit.stan, 'betas')
betas = colMeans(m$betas)
names(betas) = colnames(lStanData$X)
# compare with lm 
data.frame(coef(fit.1), betas)

s = cbind(extract(fit.stan)$betas, extract(fit.stan)$sigmaPop)
colnames(s) = c(colnames(lStanData$X), 'sigmaPop')
pairs(s, pch=20)

### partial pooling of coefficients
stanDso.2 = rstan::stan_model(file='tResponseRegression_partialPooling.stan')

m = model.matrix(CD63.Act ~ ., data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 y=dfData$CD63.Act)

fit.stan.2 = sampling(stanDso.2, data=lStanData, iter=1000, chains=2, pars=c('betas', 'mu', 'sigmaPop', 'nu', 'sigmaRan'),
                    cores=2)
print(fit.stan.2, c('betas', 'sigmaPop', 'sigmaRan', 'nu'), digits=3)

# some diagnostics for stan
pairs(fit.stan.2, pars = c("sigmaPop", "sigmaRan", 'betas[1]', 'nu', "lp__"))
pairs(fit.stan.2, pars = c("betas", "lp__"))

m = extract(fit.stan.2, 'betas')
betas.2 = colMeans(m$betas)
names(betas.2) = colnames(lStanData$X)
# compare with lm 
data.frame(coef(fit.1), betas, betas.2)

s2 = cbind(extract(fit.stan.2)$betas, extract(fit.stan.2)$sigmaPop, extract(fit.stan.2)$sigmaRan)
colnames(s2) = c(colnames(lStanData$X), 'sigmaPop', 'sigmaRan')
pairs(s2, pch=20)

### partial pooling of batches of coefficients
stanDso.2 = rstan::stan_model(file='tResponseRegression_partialPoolingBatches.stan')

m = model.matrix(CD63.Act ~ ., data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                 NBatchMap = c(1, 2, rep(3, times=9), 4, 4), NscaleBatches=4,
                 y=dfData$CD63.Act)

fit.stan.2 = sampling(stanDso.2, data=lStanData, iter=5000, chains=4, pars=c('betas', 'mu', 'sigmaPop', 'nu', 'sigmaRan'),
                      cores=4)
print(fit.stan.2, c('betas', 'sigmaPop', 'sigmaRan', 'nu'), digits=3)

traceplot(fit.stan.2, c('sigmaRan'))
traceplot(fit.stan.2, c('betas'))
# some diagnostics for stan
pairs(fit.stan.2, pars = c("sigmaPop", "sigmaRan", 'betas[1]', 'nu', "lp__"))
#pairs(fit.stan.2, pars = c("betas", "lp__"))

m = extract(fit.stan.2, 'betas')
betas.2 = colMeans(m$betas)
names(betas.2) = colnames(lStanData$X)
# compare with lm 
data.frame(coef(fit.1), betas, betas.2)

s2 = cbind(extract(fit.stan.2)$betas)
colnames(s2) = c(colnames(lStanData$X))
pairs(s2, pch=20)

###########################################################
###### 2 component finite mixture model
###########################################################
library(flexmix)
fit.flex = flexmix(CD63.Act ~ 1, data=dfData, k=2)
summary(fit.flex)
## fitted coefficients
parameters(fit.flex)
plot(density(dfData$CD63.Act))
## give initial values if you want, look at the density plot 
initf = function(chain_id = 1) {
  list(mu = c(0, 36), sigmaPop = c(1, 1), iMixWeights=c(0.5, 0.5))
} 

stanDso.3 = rstan::stan_model(file='normResponseFiniteMixture_partialPooling.stan')

m = model.matrix(CD63.Act ~ . -1, data=dfData)


lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m, 
                 iMixtures=2,
                 #iIntercepts=c(0.7, 40),
                 y=dfData$CD63.Act)

fit.stan.3 = sampling(stanDso.3, data=lStanData, iter=10000, chains=4, pars=c('betas', 'mu', 'muFitted', 'sigmaPop',
                                                                            'iMixWeights'),
                      cores=2, control=list(adapt_delta=0.99, max_treedepth = 12), init=initf)
print(fit.stan.3, c('betas', 'sigmaPop', 'iMixWeights', 'mu'), digits=3)
traceplot(fit.stan.3, c('sigmaPop', 'mu', 'iMixWeights'))
pairs(fit.stan.3, pars = c("sigmaPop", 'mu', "lp__"))

############# extract the mcmc sample values from stan
l = extract(fit.stan.3)
names(l)
mStan = cbind(l$mu, l$sigmaPop, l$iMixWeights)
dim(mStan)
colnames(mStan) = c('mu1', 'mu2', 'sigma1', 'sigma2', 'mix1', 'mix2')
dim(mStan)
colMeans(mStan)
mf = l$muFitted
dim(mf)
## get a sample for this distribution
########## simulate 200 test quantities
mDraws = matrix(NA, nrow = length(dfData$CD63.Act), ncol=200)

for (i in 1:200){
  p = sample(1:nrow(mStan), size = 1)
  mix = mean(mStan[,'mix1'])
  ## this will take a sample from a normal mixture distribution
  sam = function() {
    ind = rbinom(ncol(mf), 1, prob = mix)
    return(ind * rnorm(ncol(mf), mStan[p, 'mu1']+mf[p,], mStan[p, 'sigma1']) + 
             (1-ind) * rnorm(ncol(mf), mStan[p, 'mu2']+mf[p,], mStan[p, 'sigma2']))
  }
  mDraws[,i] = sam()
}

mDraws.normMix = mDraws

betas.3 = colMeans(extract(fit.stan.3)$betas)
data.frame(coef(fit.1), betas, betas.2, c(NA, betas.3))

###### end finite mixture model




###########################################################
##### Plot Coefficients
###########################################################
## get the coefficient of interest
mCoef = extract(fit.stan)$betas
dim(mCoef)
## get the intercept 
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]
colnames(mCoef) = colnames(lStanData$X)[2:ncol(lStanData$X)]

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
#names(m) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
text(colMeans(mCoef), ivPval, names(m), pos=1)
m = abs(m)
m = sort(m, decreasing = T)

par(mar=c(6,3,4,2)+0.1)
l2 = barplot(m[1:length(m)], 
             las=2, xaxt='n', col='grey', main='Top Variables', ylab='Absolute Coefficient')
axis(1, at = l2, labels = names(m)[1:length(m)], tick = F, las=2, cex.axis=0.7 )

## format for line plots
m = colMeans(mCoef)
#names(m) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
m = sort(m, decreasing = T)
mTreatment = mCoef[,names(m)]

df = apply(mTreatment, 2, getms)
x = 1:ncol(mTreatment)

par(p.old)
plot(x, df['m',], ylim=c(min(df), max(df)), pch=20, xlab='', main='Average Effect on Activation',
     ylab='Slopes', xaxt='n')
axis(1, at = x, labels = colnames(mTreatment), las=2, cex.axis=0.7)
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(h = 0, col='grey')

# ## export the coefficients to results file
# m = c(mean(iIntercept), colMeans(mTreatment))
# s = c(sd(iIntercept), apply(mTreatment, 2, sd))
# r = signif(cbind(m, s), 3)
# colnames(r) = c('Coefficient', 'SE')
# rownames(r)[1] = 'Intercept'
# 
# write.csv(r, file = 'results/ISACCoef.csv')
######## end plot coefficients
###########################################################

###########################################################
######## model checks for residuals
###########################################################
mFitted = extract(fit.stan.2)$mu
fitted = colMeans(mFitted)
# get residuals that is response minus fitted values
iResid = (dfData$CD63.Act - fitted)
plot(fitted, iResid, pch=20, cex=0.5)
lines(lowess(fitted, iResid), col=2, lwd=2)

## calculate standardized residuals
## these are useful to detect non-normality
## see equation 14.7 in Gelman 2013
## for t-distribution it is sqrt((scale^2)*(nu/(nu-2)))
s = mean(extract(fit.stan.2)$sigmaPop)
nu = mean(extract(fit.stan.2)$nu)
s = sqrt((s^2)*(nu/(nu-2)))
plot(fitted, iResid/s, pch=20, cex=0.5, main='standardized residuals')
lines(lowess(fitted, iResid/s), col=2, lwd=2)

## checking for non-linearity
n = colnames(dfData)
n = n[-(length(n))]
par(mfrow=c(2,2))
sapply(n, function(x){
  plot(dfData[,x], iResid, main=paste(x))
  lines(lowess(dfData[,x], iResid), col=2)
})

## unequal variances
sapply(n, function(x){
  plot(dfData[,x], abs(iResid), main=paste(x))
  lines(lowess(dfData[,x], abs(iResid)), col=2)
})

### generate some posterior predictive data
## follow the algorithm in section 14.3 page 363 in Gelman 2013
simulateOne = function(betas, sigma, mModMatrix, groupMap){
  f = mModMatrix %*% betas
  yrep = rnorm(length(f), f, sigma[groupMap])
  return(yrep)
}

runRegression = function(yrep, lStanData){
  lStanData$y = yrep
  f.s = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'mu', 'sigmaPop'),
                 cores=2)
  return(f.s)
}

# returns residuals
getResiduals = function(fit.object, mModMatrix, yrep){
  b = colMeans(extract(fit.object)$betas)
  f = mModMatrix %*% b
  r = (yrep - f)
  return(r)
}


## sample n values, 1000 times
mDraws.sim = matrix(NA, nrow = nrow(dfData), ncol=1000)
mDraws.fitted = matrix(NA, nrow = nrow(dfData), ncol=1000)
mDraws.res = matrix(NA, nrow = nrow(dfData), ncol=1000)
dim(mStan)
for (i in 1:1000){
  p = sample(1:nrow(mStan), 1)
  sigma = mStan[p,5:6]
  betas = mStan[p, 1:4]
  mDraws.sim[,i] = simulateOne(betas, sigma, lStanData$X, as.numeric(dfData$treatment))
  f.s = runRegression(mDraws.sim[,i], lStanData)
  mDraws.fitted[,i] = apply(extract(f.s)$mu, 2, mean)
  mDraws.res[,i] = getResiduals(f.s, lStanData$X,
                                mDraws.sim[,i])
}

### visual checks
ivResp = dfData$response
yresp = density(ivResp)
plot(yresp, xlab='', main='Fitted distribution', ylab='density', lwd=2)#, ylim=c(0, 1))
hist(ivResp, prob=T, main='Original Data with simulated data', xlab='Response Variable')
temp = apply(mDraws.sim, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='lightgrey', lwd=0.6)
})
lines(yresp)

### plot the residuals
plot(dfData$response, iResid, pch=c(2,20)[as.numeric(dfData$treatment)], cex=0.5, main='MCMC')
lines(lowess(dfData$response, iResid))

plot(density(iResid))
g = apply(mDraws.res, 2, function(x) lines(density(x), lwd=0.5, col=2))
lines(density(iResid))



####### end model checks residuals
###########################################################

###########################################################
######## model checks for data distribution
###########################################################
## calculate bayesian p-value for a test statistic
getPValue = function(Trep, Tobs){
  left = sum(Trep <= Tobs)/length(Trep)
  right = sum(Trep >= Tobs)/length(Trep)
  return(min(left, right))
}
## define some test quantities to measure the lack of fit
## define a test quantity T(y, theta)
## variance
T1_var = function(Y) return(var(Y))

## min quantity
T1_min = function(Y){
  return(min(Y))
} 

## max quantity
T1_max = function(Y){
  return(max(Y))
} 

## mean quantity
T1_mean = function(Y){
  return(mean(Y))
} 

## median quantity
T1_median = function(Y){
  return(median(Y))
} 

## mChecks
mChecks = matrix(NA, nrow=5, ncol=2)
rownames(mChecks) = c('Variance', 'Median', 'Max', 'Min', 'Mean')
colnames(mChecks) = c('student', 'mixture')

## generate random samples from alternative t-distribution parameterization
## see https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu

ivResp = dfData$CD63.Act
########## simulate 200 test replications of test quantities
mDraws = matrix(NA, nrow = length(ivResp), ncol=200)
mThetas = matrix(NA, nrow=200, ncol=4)
colnames(mThetas) = c('mu', 'median', 'scale', 'nu')

# sample of values from the simulation
l = extract(fit.stan.2)
names(l)

sigSample = l$sigmaPop
muSample = l$mu
nuSample = l$nu
dim(muSample)
dim(mDraws)

for (i in 1:200){
  p = sample(1:nrow(muSample), size = 1)
  s = sigSample[p]
  m = muSample[p,]
  n = nuSample[p]
  mDraws[,i] = rt_ls(length(m), n, m, s)
  mThetas[i,] = c(mean(m), median(m), s, n)
}

mDraws.t = mDraws
## get the p-values for the test statistics
t1 = apply(mDraws, 2, T1_var)
mChecks['Variance', 'student'] = getPValue(t1, var(ivResp))

## testing for outlier detection
t1 = apply(mDraws, 2, T1_min)
mChecks['Min', 'student'] = getPValue(t1, T1_min(ivResp))

## maximum value
t1 = apply(mDraws, 2, T1_max)
mChecks['Max', 'student'] = getPValue(t1, T1_max(ivResp))

## mean value
t1 = apply(mDraws, 2, T1_mean)
mChecks['Mean', 'student'] = getPValue(t1, T1_mean(ivResp))

## median value
t1 = apply(mDraws, 2, T1_median)
mChecks['Median', 'student'] = getPValue(t1, T1_median(ivResp))

mChecks

## plots of densities
par(p.old)
yresp = density(ivResp)
plot(yresp, xlab='', main='Fitted distribution', ylab='density', lwd=2)
hist(ivResp, prob=T, add=T)
temp = apply(mDraws, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='red', lwd=0.6)
})
lines(yresp, lwd=2)

hist(ivResp, prob=T)
## t samples
temp = apply(mDraws.t, 2, function(x) {x = density(x)
#x$y = x$y/max(x$y)
lines(x, col='red', lwd=0.6)
})

####### end model checks for data distribution
###########################################################

