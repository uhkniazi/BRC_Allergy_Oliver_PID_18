# File: 03_multivariate.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: multivariate analysis of the formatted data set after EDA
# Date: 2/07/2019

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
load(n[2])

# close connection after getting data
dbDisconnect(db)

str(dfData)
dfData = dfData[,-1]

## figures for showing balance and overlap
library(lattice)
summary(dfData)

################# figures for covariate distributions
dfData.bk = dfData
data.frame(colnames(dfData))
df = dfData[,-c(1,17)]
df = stack(df)
df$y = dfData$CD63.Act
df$treatment = dfData$Allergic.Status

xyplot(y ~ log(values+0.5) | ind, data=df, type=c('g', 'p', 'r'), pch=19, cex=0.6, groups=treatment,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'),
                                                   ylab='%CD63 Activation')

xyplot(y ~ values | ind, data=df, type=c('g', 'p', 'r'), pch=19, cex=0.6, groups=treatment,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'),
       ylab='%CD63 Activation')

densityplot(~ log(values+0.5) | ind, data=df, groups=treatment, auto.key=T, scales=list(relation='free'))

par(mfrow=c(2,2))

## code idea from https://stackoverflow.com/questions/3541713/how-to-plot-two-histograms-together-in-r
hist2 = function(x, y, main='', xlab='Values', ylab='Density', legends=c('Value 1', 'Value 2'), legend.pos='topright'){
  p1 = hist(x,plot=FALSE)
  p2 = hist(y,plot=FALSE)
  ## calculate the range of the graph
  xlim = range(p1$breaks,p2$breaks)
  ylim = range(0,p1$density,
                p2$density)
  plot(p1,xlim = xlim, ylim = ylim,
       col = rgb(1,0,0,0.4),xlab = xlab, ylab=ylab,
       freq = FALSE, ## relative, not absolute frequency
       main = main)
  plot(p2,xlim = xlim, ylim = ylim,
       xaxt = 'n', yaxt = 'n', ## don't add axes
       col = rgb(0,0,1,0.4), add = TRUE,
       freq = FALSE)
  legend(legend.pos, legends,
         fill = rgb(1:0,0,0:1,0.4), bty = 'n',
         border = NA)
}

cn = colnames(dfData)[-1]
f = dfData$Allergic.Status

par(mfrow=c(2,2))
for (i in seq_along(cn)){
  x = dfData[,cn[i]]
  hist2(x[f == 'PA'], x[f=='PS'], main=cn[i], legends=c('PA', 'PS'))
}

##### transformation of the covariates
tapply(df$values, df$ind, function(x) quantile(x, 0:10/10))
data.frame(colnames(dfData))
dfData = dfData[,-c(10, 11)]
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

data.frame(colnames(dfData))
df = dfData[,-c(1,15)]
df = stack(df)
df$y = dfData$CD63.Act
df$treatment = dfData$Allergic.Status

densityplot(~ values | ind, data=df, groups=treatment, auto.key=T, scales=list(relation='free'))

xyplot(y ~ values | ind, data=df[df$y > 0,], type=c('g', 'p', 'smooth'), pch=19, cex=0.6, groups=treatment,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'),
       ylab='%CD63 Activation', auto.key=list(columns=2))

xyplot(y ~ values | ind, data=df[df$y > 0,], type=c('g', 'p', 'r'), pch=19, cex=0.6, groups=treatment,
       #index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'),
       ylab='%CD63 Activation', auto.key=list(columns=2))

xyplot(y ~ values | ind, data=df, type=c('smooth'), pch=19, cex=0.6, groups=treatment,
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'))

cn = colnames(dfData)[-1]
f = dfData$Allergic.Status

par(mfrow=c(2,2))
for (i in seq_along(cn)){
  x = dfData[,cn[i]]
  hist2(x[f == 'PA'], x[f=='PS'], main=cn[i], legends=c('PA', 'PS'))
}

## center the variables before modelling
cn = cn[-(length(cn))]
dfData.ps = dfData[dfData$Allergic.Status == 'PS', ]
dfData.pa = dfData[dfData$Allergic.Status == 'PA', ]
for (i in seq_along(cn)){
  x = dfData.pa[,cn[i]]
  dfData.pa[,cn[i]] = x - mean(x)
}

for (i in seq_along(cn)){
  x = dfData.ps[,cn[i]]
  dfData.ps[,cn[i]] = x - mean(x)
}


##### fit a lm
fit.1.pa = lm(CD63.Act ~ ., data=dfData.pa[, -1])
summary(fit.1.pa)

fit.1.ps = lm(CD63.Act ~ ., data=dfData.ps[, -1])
summary(fit.1.ps)

#################### use stan to generate MCMC sample
## censor the response variable to a lower bound of 0.01
dfData.pa$CD63.Act[dfData.pa$CD63.Act <= 0] = 0.01
dfData.ps$CD63.Act[dfData.ps$CD63.Act <= 0] = 0.01

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# stanDso = rstan::stan_model(file='tResponseRegression_censored.stan')
# 
# m = model.matrix(CD63.Act ~ ., data=dfData.pa[dfData.pa$CD63.Act > 0.01, -1])
# m2 = model.matrix(CD63.Act ~ ., data=dfData.pa[dfData.pa$CD63.Act <= 0.01, -1])
# 
# lStanData = list(Ntotal=nrow(m), Ncol=ncol(m), X=m, rLower=0.01,
#                  X2 = m2,
#                  Ncens = nrow(m2),
#                  y=dfData.pa$CD63.Act[dfData.pa$CD63.Act > 0.01])
# 
# fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=2, pars=c('betas', 'mu', 'y_cens', 'mu2', 'sigmaPop', 'nu'),
#                     cores=2)
# print(fit.stan, c('betas', 'sigmaPop', 'nu'), digits=3)
# print(fit.stan, c('mu', 'mu2'))
# print(fit.stan, c('y_cens'))
# # some diagnostics for stan
# traceplot(fit.stan, c('sigmaPop', 'nu'), ncol=1, inc_warmup=F)
# 
# m = extract(fit.stan, 'betas')
# betas = colMeans(m$betas)
# names(betas) = colnames(lStanData$X)
# # compare with lm
# data.frame(coef(fit.1.ps), betas)
# 
# s = cbind(extract(fit.stan)$betas, extract(fit.stan)$sigmaPop)
# colnames(s) = c(colnames(lStanData$X), 'sigmaPop')
# pairs(s, pch=20, cex=0.5, col='grey')

### partial pooling of batches of coefficients
stanDso.2 = rstan::stan_model(file='tResponseRegression_partialPoolingBatches_censored.stan')

## model dataset PA
m = model.matrix(CD63.Act ~ ., data=dfData.pa[dfData.pa$CD63.Act > 0.01, -1])
m2 = model.matrix(CD63.Act ~ ., data=dfData.pa[dfData.pa$CD63.Act <= 0.01, -1])

lStanData = list(Ntotal=nrow(m), Ncol=ncol(m), X=m, rLower=0.01,
                 X2 = m2,
                 Ncens = nrow(m2),
                 NBatchMap = c(1, 2, 3, 4, rep(5, times=4), 6, 5, 5, 7, 7), NscaleBatches=7,
                 y=dfData.pa$CD63.Act[dfData.pa$CD63.Act > 0.01])

fit.stan.pa = sampling(stanDso.2, data=lStanData, iter=5000, chains=4, pars=c('betas', 'mu', 'y_cens', 'mu2', 'sigmaPop', 'nu', 'sigmaRan'),
                      cores=4, control=list(adapt_delta=0.99, max_treedepth = 10))

## model dataset PS
m = model.matrix(CD63.Act ~ ., data=dfData.ps[dfData.ps$CD63.Act > 0.01, -1])
m2 = model.matrix(CD63.Act ~ ., data=dfData.ps[dfData.ps$CD63.Act <= 0.01, -1])

lStanData = list(Ntotal=nrow(m), Ncol=ncol(m), X=m, rLower=0.01,
                 X2 = m2,
                 Ncens = nrow(m2),
                 NBatchMap = c(1, 2, 3, 4, rep(5, times=4), 6, 5, 5, 7, 7), NscaleBatches=7,
                 y=dfData.ps$CD63.Act[dfData.ps$CD63.Act > 0.01])

fit.stan.ps = sampling(stanDso.2, data=lStanData, iter=5000, chains=4, pars=c('betas', 'mu', 'y_cens', 'mu2', 'sigmaPop', 'nu', 'sigmaRan'),
                       cores=4, control=list(adapt_delta=0.99, max_treedepth = 10))

################# model checks
print(fit.stan.pa, c('betas', 'sigmaPop', 'sigmaRan', 'nu'), digits=3)
print(fit.stan.ps, c('betas', 'sigmaPop', 'sigmaRan', 'nu'), digits=3)

traceplot(fit.stan.pa, c('betas'))
traceplot(fit.stan.ps, c('betas'))

betas.pa = colMeans(extract(fit.stan.pa, 'betas')$betas)
betas.ps = colMeans(extract(fit.stan.ps, 'betas')$betas)
names(betas.pa) = colnames(lStanData$X)
names(betas.ps) = colnames(lStanData$X)
# compare with lm 
data.frame(coef(fit.1.pa), betas.pa, coef(fit.1.ps), betas.ps)

s2 = cbind(extract(fit.stan.pa)$betas)
colnames(s2) = c(colnames(lStanData$X))
dim(s2)
pairs(s2[sample(1:nrow(s2), 1000),], pch=20, col='grey', cex=0.5, main='PA')

s2 = cbind(extract(fit.stan.ps)$betas)
colnames(s2) = c(colnames(lStanData$X))
dim(s2)
pairs(s2[sample(1:nrow(s2), 1000),], pch=20, col='grey', cex=0.5, main='PS')


###########################################################
##### Plot Coefficients
###########################################################
## get the coefficient of interest
mCoef = extract(fit.stan.pa)$betas
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

## export the coefficients to results file
av = c(mean(iIntercept), colMeans(mTreatment))
s = c(sd(iIntercept), apply(mTreatment, 2, sd))
r = signif(cbind(av, s), 3)
colnames(r) = c('Coefficient', 'SE')
rownames(r)[1] = 'Intercept'

write.csv(r, file = 'results/cd63Coef_PA.csv')

df.pa = apply(mTreatment, 2, getms)
x.pa = 1:ncol(mTreatment)

## ps data
mCoef = extract(fit.stan.ps)$betas
dim(mCoef)
## get the intercept 
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]
colnames(mCoef) = colnames(lStanData$X)[2:ncol(lStanData$X)]

mTreatment = mCoef[,names(m)]

## export the coefficients to results file
av = c(mean(iIntercept), colMeans(mTreatment))
s = c(sd(iIntercept), apply(mTreatment, 2, sd))
r = signif(cbind(av, s), 3)
colnames(r) = c('Coefficient', 'SE')
rownames(r)[1] = 'Intercept'

write.csv(r, file = 'results/cd63Coef_PS.csv')

df.ps = apply(mTreatment, 2, getms)
x.ps = (1:ncol(mTreatment))+0.5

## line plot
par(p.old)
plot(x.pa, df.pa['m',], ylim=c(min(rbind(df.pa, df.ps)), max(rbind(df.pa, df.ps))), pch=20, xlab='', main='Coefficients for each predictor',
     ylab='Slopes', xaxt='n', xlim=c(0.5, max(x.ps)))
axis(1, at = x.pa+0.25, labels = colnames(mTreatment), las=2, cex.axis=0.7)
for(l in 1:ncol(df.pa)){
  lines(x=c(x.pa[l], x.pa[l]), y=df.pa[c(2,3),l], lwd=0.5)
}
abline(h = 0, col='grey')

points(x.ps, df.ps['m',], ylim=c(min(rbind(df.pa, df.ps)), max(rbind(df.pa, df.ps))), pch=20,
     ylab='Slopes', xaxt='n', col=2)
#axis(1, at = x.ps, labels = colnames(mTreatment), las=2, cex.axis=0.7)
for(l in 1:ncol(df.ps)){
  lines(x=c(x.ps[l], x.ps[l]), y=df.ps[c(2,3),l], lwd=0.5, col=2)
}
legend('bottomleft', c('PA', 'PS'), lty=1, col=c(1,2))

## extract the intercepts
iIntercept.pa = extract(fit.stan.pa)$betas[,1]
iIntercept.ps = extract(fit.stan.ps)$betas[,1]

df = apply(cbind(iIntercept.pa, iIntercept.ps), 2, getms)
x = 1:2

par(p.old)
plot(x, df['m',], ylim=c(min(df), max(df)), pch=20, xlab='', main='Average %CD63 Activation',
     ylab='Intercept', xaxt='n')
axis(1, at = x, labels = c('PA', 'PS'), las=2, cex.axis=1)
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}

######## end plot coefficients
###########################################################

###############################################
#################### plots for slopes
###############################################
# cn = colnames(dfData)[-c(1, 15)]
# dfPlot = rbind(dfData.pa, dfData.ps)
# dfPlot = dfPlot[dfPlot$CD63.Act > 0.01,]
# betas.pa = colMeans(extract(fit.stan.pa)$betas)
# betas.ps = colMeans(extract(fit.stan.ps)$betas)
# names(betas.pa) = colnames(lStanData$X)
# names(betas.ps) = colnames(lStanData$X)
# par(mfrow=c(2,2))
# for (i in seq_along(cn)){
#   #r = range(dfPlot[,cn[i]])
#   #x = seq(r[1], r[2], length.out = 100)
#   #y = betas[1] + betas[cn[i]] * x
#   plot(dfPlot[,cn[i]], dfPlot$CD63.Act, pch=20, xlab=paste(cn[i]), ylab='%CD63 Activation', type='n')
#   #points(dfPlot[dfPlot$Allergic.Status == 'PS',cn[i]], dfPlot$CD63.Act[dfPlot$Allergic.Status == 'PS'], pch=20, col=2)
#   abline(betas.pa[1], betas.pa[cn[i]])
#   abline(betas.ps[1], betas.ps[cn[i]], col=2)
#   #lines(x, y, col=2)
# }
###############################################


###########################################################
######## model checks for residuals
###########################################################
### variance check/residual sd check in each group PS and PA
fitted.pa = colMeans(extract(fit.stan.pa)$mu)
fitted.ps = colMeans(extract(fit.stan.ps)$mu)
# get residuals that is response minus fitted values
iResid.pa = (dfData.pa$CD63.Act[dfData.pa$CD63.Act > 0.01] - fitted.pa)
iResid.ps = (dfData.ps$CD63.Act[dfData.ps$CD63.Act > 0.01] - fitted.ps)
## calculate standardized residuals
## these are useful to detect non-normality
## see equation 14.7 in Gelman 2013
## for t-distribution it is sqrt((scale^2)*(nu/(nu-2)))
s = mean(extract(fit.stan.pa)$sigmaPop)
nu = mean(extract(fit.stan.pa)$nu)
sa = sqrt((s^2)*(nu/(nu-2)))

s = mean(extract(fit.stan.ps)$sigmaPop)
nu = mean(extract(fit.stan.ps)$nu)
ss = sqrt((s^2)*(nu/(nu-2)))

iResid = c(iResid.pa/sa, iResid.ps/ss)
l = c(rep(1, times=length(iResid.pa)), rep(2, times=length(iResid.ps)))
plot(l, abs(iResid), xaxt='n', ylab='Absolute Standardised Residuals', main='Unequal Variance in each group')
axis(1, at = c(1,2), labels = c('PA', 'PS'))

# get residuals that is response minus fitted values
plot(fitted.pa, iResid.pa, pch=20, cex=0.5)
lines(lowess(fitted.pa, iResid.pa), col=2, lwd=2)
plot(dfData.pa$CD63.Act[dfData.pa$CD63.Act > 0.01], fitted.pa, pch=20, xlab='observed')

## checking for non-linearity
n = colnames(dfData[,-1])
n = n[-(length(n))]
par(mfrow=c(2,2))
sapply(n, function(x){
  plot(dfData.pa[dfData.pa$CD63.Act > 0.01 ,x], iResid.pa/sa, main=paste(x))
  lines(lowess(dfData.pa[dfData.pa$CD63.Act > 0.01,x], iResid.pa/sa), col=2)
})

# ## unequal variances
# sapply(n, function(x){
#   plot(dfData.ps[dfData.ps$CD63.Act > 0.01,x], abs(iResid), main=paste(x))
#   lines(lowess(dfData.ps[dfData.ps$CD63.Act > 0.01,x], abs(iResid)), col=2)
# })

### generate some posterior predictive data
## generate random samples from alternative t-distribution parameterization
## see https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
## follow the algorithm in section 14.3 page 363 in Gelman 2013
simulateOne = function(betas, sigma, nu, mModMatrix){
  f = mModMatrix %*% betas
  yrep = rt_ls(length(f), nu, f,  sigma)
  # censor the values below detection limit
  yrep[yrep <= 0.01] = 0.01
  return(yrep)
}

runRegression = function(yrep, mModMatrix){
  f = which(yrep <= 0.01)
  m = mModMatrix[-f,]
  m2 = mModMatrix[f,]
  
  lStanData = list(Ntotal=nrow(m), Ncol=ncol(m), X=m, rLower=0.01,
                   X2 = m2,
                   Ncens = nrow(m2),
                   NBatchMap =c(1, 2, 3, 4, rep(5, times=4), 6, 5, 5, 7, 7), NscaleBatches=7,
                   y=yrep[-f])
  f.s = sampling(stanDso.2, data=lStanData, iter=1000, chains=2, pars=c('betas', 'mu', 'nu', 'sigmaPop'),
                        cores=2)
  iFitted = colMeans(extract(f.s)$mu)
  nu = mean(extract(f.s)$nu)
  sig = mean(extract(f.s)$sigmaPop)
  sig = sqrt((sig^2)*(nu/(nu-2)))
  res = yrep[-f] - iFitted
  return(list(stan=f.s, res=res, res.sd=res/sig, sig=sig, yrep=yrep))
}

### test functions
# l = extract(fit.stan.2)
# names(l)
# s = simulateOne(colMeans(l$betas), 
#                 mean(l$sigmaPop),
#                 mean(l$nu),
#                 model.matrix(CD63.Act ~ ., data=dfData))
# 
# r = runRegression(s, model.matrix(CD63.Act ~ ., data=dfData))

### test works, continue with simulation
## sample n values, 1000 times
mDraws.sim = matrix(NA, nrow = nrow(dfData.pa), ncol=300)
lSims = vector(mode = 'list', length = 300)
l = extract(fit.stan.pa)
for (i in 1:300){
  p = sample(1:nrow(l$betas), 1)
  mDraws.sim[,i] = simulateOne(l$betas[p,], 
                  l$sigmaPop[p],
                  l$nu[p],
                  model.matrix(CD63.Act ~ ., data=dfData.pa[,-1]))
  
  lSims[[i]] = runRegression(mDraws.sim[,i], model.matrix(CD63.Act ~ ., data=dfData.pa[,-1]))
}

## visual distribution check
hist(dfData.pa$CD63.Act, prob=T)
plot(density(dfData.pa$CD63.Act))
apply(mDraws.sim, 2, function(x) lines(density(x)))

iResLimit = ceiling(abs(sa * 3))

## proportion of absolute value residuals that exceed limit iResLimit
T1_proportion = function(Y) {
  return(sum(abs(Y) > iResLimit)/length(Y))
}

iOutliers = sapply(lSims, function(x){
  return(T1_proportion(x$res))
})

iObserved = T1_proportion(iResid.pa)
hist(iOutliers)
points(iObserved, 0, col=2)

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


getPValue(na.omit(iOutliers), iObserved)

## get the p-values for the test statistics
t1 = apply(mDraws.sim, 2, T1_var)
getPValue(t1, var(dfData.pa$CD63.Act))

## testing for outlier detection
t1 = apply(mDraws.sim, 2, T1_min)
getPValue(t1, T1_min(dfData.pa$CD63.Act))

## maximum value
t1 = apply(mDraws.sim, 2, T1_max)
getPValue(t1, T1_max(dfData.pa$CD63.Act))

## mean value
t1 = apply(mDraws.sim, 2, T1_mean)
getPValue(t1, T1_mean(dfData.pa$CD63.Act))

## median value
t1 = apply(mDraws.sim, 2, T1_median)
getPValue(t1, T1_median(dfData.pa$CD63.Act))

## plot the relationship with covariates and simulated data
pdf('results/simulated_vs_covariates.pdf')
par(mfrow=c(2,2))
plot(dfData.pa$Age, dfData.pa$CD63.Act, pch=20, col=2, xlab='Age', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$Age, x), lwd=0.5, col='grey'))

plot(dfData.pa$Peanut.SPT, dfData.pa$CD63.Act, pch=20, col=2, xlab='Peanut.SPT', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$Peanut.SPT, x), lwd=0.5, col='grey'))

plot(dfData.pa$total.IgE, dfData.pa$CD63.Act, pch=20, col=2, xlab='total.IgE', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$total.IgE, x), lwd=0.5, col='grey'))

plot(dfData.pa$f13.Peanut, dfData.pa$CD63.Act, pch=20, col=2, xlab='f13.Peanut', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$f13.Peanut, x), lwd=0.5, col='grey'))

plot(dfData.pa$f422.rAra.h.1, dfData.pa$CD63.Act, pch=20, col=2, xlab='f422.rAra.h.1', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$f422.rAra.h.1, x), lwd=0.5, col='grey'))

plot(dfData.pa$f423.Ara.h.2, dfData.pa$CD63.Act, pch=20, col=2, xlab='f423.Ara.h.2', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$f423.Ara.h.2, x), lwd=0.5, col='grey'))

plot(dfData.pa$f424.rAra.h.3, dfData.pa$CD63.Act, pch=20, col=2, xlab='f424.rAra.h.3', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$f424.rAra.h.3, x), lwd=0.5, col='grey'))

plot(dfData.pa$f423.nAra.H.6, dfData.pa$CD63.Act, pch=20, col=2, xlab='f423.nAra.H.6', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$f423.nAra.H.6, x), lwd=0.5, col='grey'))

plot(dfData.pa$Peanut.Sp.Act, dfData.pa$CD63.Act, pch=20, col=2, xlab='Peanut.Sp.Act', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$Peanut.Sp.Act, x), lwd=0.5, col='grey'))

plot(dfData.pa$Ara.h.2.Sp.Act, dfData.pa$CD63.Act, pch=20, col=2, xlab='Ara.h.2.Sp.Act', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$Ara.h.2.Sp.Act, x), lwd=0.5, col='grey'))

plot(dfData.pa$Ara.h.6.Sp.Act, dfData.pa$CD63.Act, pch=20, col=2, xlab='Ara.h.6.Sp.Act', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$Ara.h.6.Sp.Act, x), lwd=0.5, col='grey'))

plot(dfData.pa$ISAC.Shannon, dfData.pa$CD63.Act, pch=20, col=2, xlab='ISAC.Shannon', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$ISAC.Shannon, x), lwd=0.5, col='grey'))

plot(dfData.pa$Peanut.Shannon, dfData.pa$CD63.Act, pch=20, col=2, xlab='Peanut.Shannon', ylab="%CD63 Activation")
apply(mDraws.sim, 2, function(x) lines(lowess(dfData.pa$Peanut.Shannon, x), lwd=0.5, col='grey'))

dev.off(dev.cur())

#### simulate data for some covariates, holding others at their average value i.e. 0
mDraws.sim = matrix(NA, nrow = nrow(dfData.pa), ncol=500)
l = extract(fit.stan.pa)
m = model.matrix(CD63.Act ~ ., data=dfData.pa[,-1])
i = which(colnames(m) %in% c('(Intercept)', 'Peanut.Sp.Act', 'ISAC.Shannon', 'Peanut.Shannon', 'Age'))
m[,-i] = 0
for (i in 1:500){
  p = sample(1:nrow(l$betas), 1)
  mDraws.sim[,i] = simulateOne(l$betas[p,], 
                               l$sigmaPop[p],
                               l$nu[p],
                               m)
}



####### end model checks residuals
###########################################################

