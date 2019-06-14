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
df = read.csv(n, header=T)

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
save(dfData, file=n2)

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

xyplot(y ~ log(values+0.5) | ind, data=df, type=c('smooth'), pch=19, cex=0.6,
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'))

xyplot(y ~ values | ind, data=df, type=c('smooth'), pch=19, cex=0.6,
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'))

## transformation of the covariates
fit.1 = lm(CD63.Act ~ Peanut.Sp.Act, data=dfData)
fit.2 = lm(CD63.Act ~ log(Peanut.Sp.Act), data=dfData)
summary(fit.1)
summary(fit.2)

s1 = simulate(fit.1, 20)
s2 = simulate(fit.2, 20)
par(mfrow=c(2,3))
plot(density(dfData$CD63.Act))

plot(density(simulate(fit.1, 1)[,1]))
apply(s1, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
lines(density(dfData$CD63.Act), col=2)

plot(density(simulate(fit.2, 1)[,1]))
apply(s2, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
lines(density(dfData$CD63.Act), col=2)

dfData.new = dfData
dfData.new$Peanut.Sp.Act = log(dfData$Peanut.Sp.Act)

## repeat for each covariate
fit.1 = lm(CD63.Act ~ Peanut.SPT, data=dfData)
fit.2 = lm(CD63.Act ~ log(Peanut.SPT), data=dfData)
summary(fit.1)
summary(fit.2)

s1 = simulate(fit.1, 20)
s2 = simulate(fit.2, 20)
par(mfrow=c(2,3))
plot(density(dfData$CD63.Act))

plot(density(simulate(fit.1, 1)[,1]))
apply(s1, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
lines(density(dfData$CD63.Act), col=2)

plot(density(simulate(fit.2, 1)[,1]))
apply(s2, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
lines(density(dfData$CD63.Act), col=2)

###
fit.1 = lm(CD63.Act ~ total.IgE, data=dfData)
fit.2 = lm(CD63.Act ~ log(total.IgE), data=dfData)
summary(fit.1)
summary(fit.2)

s1 = simulate(fit.1, 20)
s2 = simulate(fit.2, 20)
par(mfrow=c(2,3))
plot(density(dfData$CD63.Act))

plot(density(simulate(fit.1, 1)[,1]))
apply(s1, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
lines(density(dfData$CD63.Act), col=2)

plot(density(simulate(fit.2, 1)[,1]))
apply(s2, 2, function(x) lines(density(x), lwd=0.5, col='grey'))
lines(density(dfData$CD63.Act), col=2)

dfData.new$total.IgE = log(dfData$total.IgE)


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

###########################################################
##### Plot Coefficients
###########################################################
## get the coefficient of interest
mCoef = extract(fit.stan.2)$betas
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
######## model checks
###########################################################
mFitted = extract(fit.stan)$mu
fitted = colMeans(mFitted)
# get residuals that is response minus fitted values
iResid = (dfData$CD63.Act - fitted)
plot(fitted, iResid, pch=20, cex=0.5)
lines(lowess(fitted, iResid), col=2, lwd=2)

## calculate standardized residuals
## these are useful to detect non-normality
## see equation 14.7 in Gelman 2013
s = mean(extract(fit.stan)$sigmaPop)
plot(fitted, iResid/s, pch=20, cex=0.5, main='standardized residuals')
lines(lowess(fitted, iResid/s), col=2, lwd=2)

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



####### end model checks
###########################################################



