# File: 05_mechanisticWithAvidity.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: mechanistic EDA with the CD63 activation data and Avidity data
# Date: 7/1/2020


source('header.R')
dfData = read.csv('dataExternal/Mechanistic Data - CD63 Activation Assay with avidity.csv', header=T,
                  stringsAsFactors = F)

## clean up the input data
df = dfData
df$Patient.ID = factor(gsub(' ', '', as.character(df$Patient.ID)))
df$Allergic.Status = gsub(' ', '', as.character(df$Allergic.Status))
df$Allergic.Status[is.na(df$Allergic.Status)] = 'Healthy'
df$Allergic.Status = factor(df$Allergic.Status, levels = c('Healthy', 'PS', 'PA'))
df = droplevels.data.frame(df)
dfData = df

## check where is imputation required
apply(df, 2, function(x) sum(is.na(x)))
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
dfData[c('25', '98'), 'Peanut.SPT'] = 3.55
dfData[c('100'), 'Peanut.SPT'] = 9.69

f = is.na(dfData$Age)
table(f)
tapply(dfData$Age, dfData$Allergic.Status, mean, na.rm=T)
f = which(is.na(dfData$Age))
dfData[f,c('Allergic.Status', 'Age')]
dfData[c('98'), 'Age'] = 5.88
dfData[c('100'), 'Age'] = 8.87

summary(dfData)

f = complete.cases(dfData)
table(f)
## drop missing data, this will be the avidity data
dfData = dfData[f,]

## save this data for further analysis
n = make.names(paste('Imputed data for mechanistic analysis with avidity variable data ID 40 rds'))
n2 = paste0('~/Data/MetaData/', n)
#save(dfData, file=n2)

library(RMySQL)
## note: comment out as this entry has been made in db
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=40, name=n, type='rds', location='~/Data/MetaData/', comment='Olivers allergen data after imputation for PS and PA groups with avidity data added')
#dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)

################# preliminary model fitting
df = stack(dfData[,-c(1, 2, 18)])
df$y = dfData$CD63.Act
df$g = dfData$Allergic.Status

library(lattice)
xyplot(y ~ log(values+1) | ind, data=df, groups=g, type=c('g', 'p', 'r'), pch=19, cex=0.6,
       index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'),
       auto.key=list(columns=2))

xyplot(y ~ log(values+1) | ind, data=df, groups=g, type=c('g', 'p', 'smooth'), pch=19, cex=0.6,
       index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'),
       auto.key=list(columns=2))

xyplot(y ~ values | ind, data=df, groups=g, type=c('g', 'p', 'smooth'), pch=19, cex=0.6,
       index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'),
       auto.key=list(columns=2))


as.data.frame(colnames(dfData))
m = as.matrix(dfData[,c(4:12)])
pairs(log(m[dfData$Allergic.Status == 'PA', ]), pch=20)

tapply(df$values[df$g == 'PA'], df$ind[df$g == 'PA'], function(x) quantile(x, 0:10/10))
dfData.bk = dfData

dfData = dfData[dfData$Allergic.Status == 'PA', ]
dfData = droplevels.data.frame(dfData)
library(rethinking)
library(car)

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
########################################################################################
##### peanut specific ige and avidity
########################################################################################
## see pages 105 to 109 (rethinking book) for extract.samples, link and sim for 
## posterior predictive checks 
str(dfData)
dfData$CD63_Act = dfData$CD63.Act - mean(dfData$CD63.Act)
# dfData$f13_Peanut = dfData$f13.Peanut - mean(dfData$f13.Peanut)
dfData$f13_Peanut = log(dfData$f13.Peanut) 
dfData$f13_Peanut = dfData$f13_Peanut - mean(dfData$f13_Peanut)
dfData$Avidity = dfData$Avidity - mean(dfData$Avidity)

fit.lm = lm(CD63_Act ~ f13_Peanut + Avidity, data=dfData)
summary(fit.lm)

model.1 <- alist(
  # https://github.com/rmcelreath/rethinking/blob/ae3fae963e244c58bbd68dd86478840a87d29f49/man/dstudent.Rd#L37
  CD63_Act ~ dstudent(3, mu, sigmaPop),
  mu <- b0 + b_f13_Peanut * f13_Peanut +
    b_Avidity * Avidity,
  b0 ~ dnorm(0, 2),
  c(b_f13_Peanut, b_Avidity) ~ dnorm(0, 3),
  sigmaPop ~ dexp(1/28)
  # see here for example
  # https://github.com/rmcelreath/rethinking/blob/ae3fae963e244c58bbd68dd86478840a87d29f49/tests/rethinking_tests/test_ulam1.R#L20
  #nu ~ gamma(3, 0.1)
)

fit.1 <- quap(model.1,
              data=dfData,
              start=list(b0=0)
)

summary(fit.1)
plot(coeftab(fit.1), pars=c('b0', 'b_f13_Peanut', 'b_Avidity'))

# fit.1.u = ulam(model.1, data=as.list(dfData[,c('CD63_Act', 'f13_Peanut', 'Avidity')]), log_lik = T,
#                start = list(sigmaPop=14, nu=10), chains = 4, cores=4)
m = extract.samples(fit.1, n=1000)
names(m)
pairs(m, pch=20)
round(apply(m, 2, getDifference),3)

### residual checks
fitted.pa = link(fit.1, data = as.list(dfData), n = 200)
dim(fitted.pa)
fitted.pa = colMeans(fitted.pa)

## calculate standardized residuals
## these are useful to detect non-normality
## see equation 14.7 in Gelman 2013
## for t-distribution it is sqrt((scale^2)*(nu/(nu-2)))
s = coef(fit.1)['sigmaPop']
nu = 3
sa = sqrt((s^2)*(nu/(nu-2)))
# get residuals that is response minus fitted values
iResid.pa = (dfData$CD63_Act - fitted.pa) / sa
plot(fitted.pa, iResid.pa, pch=20)
lines(lowess(fitted.pa, iResid.pa))

## checking for non-linearity
plot(dfData$f13_Peanut, iResid.pa)
lines(lowess(dfData$f13_Peanut, iResid.pa))

plot(dfData$Avidity, iResid.pa)
lines(lowess(dfData$Avidity, iResid.pa))

## unequal variances
plot(dfData$f13_Peanut, abs(iResid.pa))
lines(lowess(dfData$f13_Peanut, abs(iResid.pa)))

plot(dfData$Avidity, abs(iResid.pa))
lines(lowess(dfData$Avidity, abs(iResid.pa)))

### use link and sim functions to get samples of fitted values and posterior predictive data
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
simulateOne = function(betas, sigma, nu, mModMatrix){
  f = mModMatrix %*% betas
  yrep = rt_ls(length(f), nu, f,  sigma)
  return(yrep)
}
mDraws = (sim(fit.1, data=as.list(dfData), n=200))
mMuSim = (link(fit.1, data=as.list(dfData), n=200))
dim(mDraws)
dim(mMuSim)
mParam = extract.samples(fit.1, n=200)
dim(mParam)

mDraws.2 = matrix(NA, nrow = 200, ncol = 36)
for (i in 1:200){
  mDraws.2[i,] = simulateOne(t(mParam[i,c(1,2,3)]), mParam[i,4], 3, model.matrix(CD63_Act ~ f13_Peanut + Avidity, data=dfData)) 
}

## visual checks
plot(density(dfData$CD63_Act))
lines(density(colMeans(mDraws)))
apply(mDraws, 1, function(x) lines(density(x), lwd=0.5, col='grey'))

plot(density(dfData$CD63_Act))
lines(density(colMeans(mDraws.2)))
apply(mDraws.2, 1, function(x) lines(density(x), lwd=0.5, col='grey'))

plot(density(colMeans(mDraws)), xlim=c(-40, 40))
lines(density(colMeans(mDraws.2)))
apply(mDraws, 1, function(x) lines(density(x), lwd=0.5, col='grey'))
apply(mDraws.2, 1, function(x) lines(density(x), lwd=0.5, col='red'))

#### plot covariates vs actual data and fitted values
plot(dfData$f13_Peanut, dfData$CD63_Act, pch=20,
     xlab='log f13_peanut', ylab='%CD63 Act',
     main='Relationship of Predictor to Response')
lines(lowess(dfData$f13_Peanut, colMeans(mMuSim)))
points(dfData$f13_Peanut, colMeans(mMuSim), col=2, pch=20)
lines(lowess(dfData$f13_Peanut, dfData$CD63_Act))

i = range(dfData$f13_Peanut)
iGrid = seq(i[1], i[2], length.out = 50)
## hold second variable at average
coef(fit.1)
mFitted = link(fit.1, data=list(f13_Peanut=iGrid, Avidity=rep(0, times=50)))
## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

### repeat for second variable
plot(dfData$Avidity, dfData$CD63_Act, pch=20,
     xlab='Avidity', ylab='%CD63 Act',
     main='Relationship of Predictor to Response')
lines(lowess(dfData$Avidity, colMeans(mMuSim)))
points(dfData$Avidity, colMeans(mMuSim), col=2, pch=20)
lines(lowess(dfData$Avidity, dfData$CD63_Act))

i = range(dfData$Avidity)
iGrid = seq(i[1], i[2], length.out = 50)
## hold second variable at average
coef(fit.1)
mFitted = link(fit.1, data=list(f13_Peanut=rep(0, times=50), Avidity=iGrid))
## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

######################################################################################

########################################################################################
##### peanut specific ige Ara h2 and avidity
########################################################################################
## see pages 105 to 109 (rethinking book) for extract.samples, link and sim for 
## posterior predictive checks 
str(dfData)
dfData$CD63_Act = dfData$CD63.Act - mean(dfData$CD63.Act)
dfData$f423_Arah2 = log(dfData$f423.Ara.h.2 + 1) 
dfData$f423_Arah2 = dfData$f423_Arah2 - mean(dfData$f423_Arah2)
#dfData$Avidity = dfData$Avidity - mean(dfData$Avidity)

fit.lm = lm(CD63_Act ~ f423_Arah2 + Avidity, data=dfData)
summary(fit.lm)

model.1 <- alist(
  # https://github.com/rmcelreath/rethinking/blob/ae3fae963e244c58bbd68dd86478840a87d29f49/man/dstudent.Rd#L37
  CD63_Act ~ dstudent(3, mu, sigmaPop),
  mu <- b0 + b_f423_Arah2 * f423_Arah2 +
    b_Avidity * Avidity,
  b0 ~ dnorm(0, 2),
  c(b_f423_Arah2, b_Avidity) ~ dnorm(0, 3),
  sigmaPop ~ dexp(1/28)
  # see here for example
  # https://github.com/rmcelreath/rethinking/blob/ae3fae963e244c58bbd68dd86478840a87d29f49/tests/rethinking_tests/test_ulam1.R#L20
  #nu ~ gamma(3, 0.1)
)

fit.1 <- quap(model.1,
              data=dfData,
              start=list(b0=0)
)

summary(fit.1)
plot(coeftab(fit.1), pars=c('b0', 'b_f423_Arah2', 'b_Avidity'))

m = extract.samples(fit.1, n=1000)
names(m)
pairs(m, pch=20)
round(apply(m, 2, getDifference),3)

### residual checks
fitted.pa = link(fit.1, data = as.list(dfData), n = 200)
dim(fitted.pa)
fitted.pa = colMeans(fitted.pa)

## calculate standardized residuals
## these are useful to detect non-normality
## see equation 14.7 in Gelman 2013
## for t-distribution it is sqrt((scale^2)*(nu/(nu-2)))
s = coef(fit.1)['sigmaPop']
nu = 3
sa = sqrt((s^2)*(nu/(nu-2)))
# get residuals that is response minus fitted values
iResid.pa = (dfData$CD63_Act - fitted.pa) / sa
plot(fitted.pa, iResid.pa, pch=20)
lines(lowess(fitted.pa, iResid.pa))

## checking for non-linearity
plot(dfData$f423_Arah2, iResid.pa)
lines(lowess(dfData$f423_Arah2, iResid.pa))

plot(dfData$Avidity, iResid.pa)
lines(lowess(dfData$Avidity, iResid.pa))

## unequal variances
plot(dfData$f423_Arah2, abs(iResid.pa))
lines(lowess(dfData$f423_Arah2, abs(iResid.pa)))

plot(dfData$Avidity, abs(iResid.pa))
lines(lowess(dfData$Avidity, abs(iResid.pa)))

### use link and sim functions to get samples of fitted values and posterior predictive data
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
simulateOne = function(betas, sigma, nu, mModMatrix){
  f = mModMatrix %*% betas
  yrep = rt_ls(length(f), nu, f,  sigma)
  return(yrep)
}
mDraws = (sim(fit.1, data=as.list(dfData), n=200))
mMuSim = (link(fit.1, data=as.list(dfData), n=200))
dim(mDraws)
dim(mMuSim)
mParam = extract.samples(fit.1, n=200)
dim(mParam)

mDraws.2 = matrix(NA, nrow = 200, ncol = 36)
for (i in 1:200){
  mDraws.2[i,] = simulateOne(t(mParam[i,c(1,2,3)]), mParam[i,4], 3, model.matrix(CD63_Act ~ f423_Arah2 + Avidity, data=dfData)) 
}

## visual checks
plot(density(dfData$CD63_Act))
lines(density(colMeans(mDraws)))
apply(mDraws, 1, function(x) lines(density(x), lwd=0.5, col='grey'))

plot(density(dfData$CD63_Act))
lines(density(colMeans(mDraws.2)))
apply(mDraws.2, 1, function(x) lines(density(x), lwd=0.5, col='grey'))

plot(density(colMeans(mDraws)), xlim=c(-40, 40))
lines(density(colMeans(mDraws.2)))
apply(mDraws, 1, function(x) lines(density(x), lwd=0.5, col='grey'))
apply(mDraws.2, 1, function(x) lines(density(x), lwd=0.5, col='red'))

#### plot covariates vs actual data and fitted values
plot(dfData$f423_Arah2, dfData$CD63_Act, pch=20,
     xlab='log f423_Arah2', ylab='%CD63 Act',
     main='Relationship of Predictor to Response')
lines(lowess(dfData$f423_Arah2, colMeans(mMuSim)))
points(dfData$f423_Arah2, colMeans(mMuSim), col=2, pch=20)
lines(lowess(dfData$f423_Arah2, dfData$CD63_Act))

i = range(dfData$f423_Arah2)
iGrid = seq(i[1], i[2], length.out = 50)
## hold second variable at average
coef(fit.1)
mFitted = link(fit.1, data=list(f423_Arah2=iGrid, Avidity=rep(0, times=50)))
## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

### repeat for second variable
plot(dfData$Avidity, dfData$CD63_Act, pch=20,
     xlab='Avidity', ylab='%CD63 Act',
     main='Relationship of Predictor to Response')
lines(lowess(dfData$Avidity, colMeans(mMuSim)))
points(dfData$Avidity, colMeans(mMuSim), col=2, pch=20)
lines(lowess(dfData$Avidity, dfData$CD63_Act))

i = range(dfData$Avidity)
iGrid = seq(i[1], i[2], length.out = 50)
## hold second variable at average
coef(fit.1)
mFitted = link(fit.1, data=list(f423_Arah2=rep(0, times=50), Avidity=iGrid))
## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

######################################################################################

########################################################################################
##### peanut specific ige Ara h3 and avidity
########################################################################################
## see pages 105 to 109 (rethinking book) for extract.samples, link and sim for 
## posterior predictive checks 
str(dfData)
dfData$CD63_Act = dfData$CD63.Act - mean(dfData$CD63.Act)
dfData$f424_rArah3 = log(dfData$f424.rAra.h.3 + 1) 
dfData$f424_rArah3 = dfData$f424_rArah3 - mean(dfData$f424_rArah3)

fit.lm = lm(CD63_Act ~ f424_rArah3 + Avidity, data=dfData)
summary(fit.lm)

model.1 <- alist(
  # https://github.com/rmcelreath/rethinking/blob/ae3fae963e244c58bbd68dd86478840a87d29f49/man/dstudent.Rd#L37
  CD63_Act ~ dstudent(3, mu, sigmaPop),
  mu <- b0 + b_f424_rArah3 * f424_rArah3 +
    b_Avidity * Avidity,
  b0 ~ dnorm(0, 2),
  c(b_f424_rArah3, b_Avidity) ~ dnorm(0, 3),
  sigmaPop ~ dexp(1/28)
  # see here for example
  # https://github.com/rmcelreath/rethinking/blob/ae3fae963e244c58bbd68dd86478840a87d29f49/tests/rethinking_tests/test_ulam1.R#L20
  #nu ~ gamma(3, 0.1)
)

fit.1 <- quap(model.1,
              data=dfData,
              start=list(b0=0)
)

summary(fit.1)
plot(coeftab(fit.1), pars=c('b0', 'b_f424_rArah3', 'b_Avidity'))

m = extract.samples(fit.1, n=1000)
names(m)
pairs(m, pch=20)
round(apply(m, 2, getDifference),3)

### residual checks
fitted.pa = link(fit.1, data = as.list(dfData), n = 200)
dim(fitted.pa)
fitted.pa = colMeans(fitted.pa)

## calculate standardized residuals
## these are useful to detect non-normality
## see equation 14.7 in Gelman 2013
## for t-distribution it is sqrt((scale^2)*(nu/(nu-2)))
s = coef(fit.1)['sigmaPop']
nu = 3
sa = sqrt((s^2)*(nu/(nu-2)))
# get residuals that is response minus fitted values
iResid.pa = (dfData$CD63_Act - fitted.pa) / sa
plot(fitted.pa, iResid.pa, pch=20)
lines(lowess(fitted.pa, iResid.pa))

## checking for non-linearity
plot(dfData$f424_rArah3, iResid.pa)
lines(lowess(dfData$f424_rArah3, iResid.pa))

plot(dfData$Avidity, iResid.pa)
lines(lowess(dfData$Avidity, iResid.pa))

## unequal variances
plot(dfData$f424_rArah3, abs(iResid.pa))
lines(lowess(dfData$f424_rArah3, abs(iResid.pa)))

plot(dfData$Avidity, abs(iResid.pa))
lines(lowess(dfData$Avidity, abs(iResid.pa)))

### use link and sim functions to get samples of fitted values and posterior predictive data
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
simulateOne = function(betas, sigma, nu, mModMatrix){
  f = mModMatrix %*% betas
  yrep = rt_ls(length(f), nu, f,  sigma)
  return(yrep)
}
mDraws = (sim(fit.1, data=as.list(dfData), n=200))
mMuSim = (link(fit.1, data=as.list(dfData), n=200))
dim(mDraws)
dim(mMuSim)
mParam = extract.samples(fit.1, n=200)
dim(mParam)

mDraws.2 = matrix(NA, nrow = 200, ncol = 36)
for (i in 1:200){
  mDraws.2[i,] = simulateOne(t(mParam[i,c(1,2,3)]), mParam[i,4], 3, model.matrix(CD63_Act ~ f424_rArah3 + Avidity, data=dfData)) 
}

## visual checks
plot(density(dfData$CD63_Act))
lines(density(colMeans(mDraws)))
apply(mDraws, 1, function(x) lines(density(x), lwd=0.5, col='grey'))

plot(density(dfData$CD63_Act))
lines(density(colMeans(mDraws.2)))
apply(mDraws.2, 1, function(x) lines(density(x), lwd=0.5, col='grey'))

plot(density(colMeans(mDraws)), xlim=c(-40, 40))
lines(density(colMeans(mDraws.2)))
apply(mDraws, 1, function(x) lines(density(x), lwd=0.5, col='grey'))
apply(mDraws.2, 1, function(x) lines(density(x), lwd=0.5, col='red'))

#### plot covariates vs actual data and fitted values
plot(dfData$f424_rArah3, dfData$CD63_Act, pch=20,
     xlab='log f424_rArah3', ylab='%CD63 Act',
     main='Relationship of Predictor to Response')
lines(lowess(dfData$f424_rArah3, colMeans(mMuSim)))
points(dfData$f424_rArah3, colMeans(mMuSim), col=2, pch=20)
lines(lowess(dfData$f424_rArah3, dfData$CD63_Act))

i = range(dfData$f424_rArah3)
iGrid = seq(i[1], i[2], length.out = 50)
## hold second variable at average
coef(fit.1)
mFitted = link(fit.1, data=list(f424_rArah3=iGrid, Avidity=rep(0, times=50)))
## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

### repeat for second variable
plot(dfData$Avidity, dfData$CD63_Act, pch=20,
     xlab='Avidity', ylab='%CD63 Act',
     main='Relationship of Predictor to Response')
lines(lowess(dfData$Avidity, colMeans(mMuSim)))
points(dfData$Avidity, colMeans(mMuSim), col=2, pch=20)
lines(lowess(dfData$Avidity, dfData$CD63_Act))

i = range(dfData$Avidity)
iGrid = seq(i[1], i[2], length.out = 50)
## hold second variable at average
coef(fit.1)
mFitted = link(fit.1, data=list(f424_rArah3=rep(0, times=50), Avidity=iGrid))
## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

######################################################################################


########################################################################################
##### peanut specific ige Ara h6 and avidity
########################################################################################
## see pages 105 to 109 (rethinking book) for extract.samples, link and sim for 
## posterior predictive checks 
str(dfData)
dfData$CD63_Act = dfData$CD63.Act - mean(dfData$CD63.Act)
dfData$f423_nAraH6 = log(dfData$f423.nAra.H.6 + 1) 
dfData$f423_nAraH6 = dfData$f423_nAraH6 - mean(dfData$f423_nAraH6)

fit.lm = lm(CD63_Act ~ f423_nAraH6 + Avidity, data=dfData)
summary(fit.lm)

model.1 <- alist(
  # https://github.com/rmcelreath/rethinking/blob/ae3fae963e244c58bbd68dd86478840a87d29f49/man/dstudent.Rd#L37
  CD63_Act ~ dstudent(3, mu, sigmaPop),
  mu <- b0 + b_f423_nAraH6 * f423_nAraH6 +
    b_Avidity * Avidity,
  b0 ~ dnorm(0, 2),
  c(b_f423_nAraH6, b_Avidity) ~ dnorm(0, 3),
  sigmaPop ~ dexp(1/28)
  # see here for example
  # https://github.com/rmcelreath/rethinking/blob/ae3fae963e244c58bbd68dd86478840a87d29f49/tests/rethinking_tests/test_ulam1.R#L20
  #nu ~ gamma(3, 0.1)
)

fit.1 <- quap(model.1,
              data=dfData,
              start=list(b0=0)
)

summary(fit.1)
plot(coeftab(fit.1), pars=c('b0', 'b_f423_nAraH6', 'b_Avidity'))

m = extract.samples(fit.1, n=1000)
names(m)
pairs(m, pch=20)
round(apply(m, 2, getDifference),3)

### residual checks
fitted.pa = link(fit.1, data = as.list(dfData), n = 200)
dim(fitted.pa)
fitted.pa = colMeans(fitted.pa)

## calculate standardized residuals
## these are useful to detect non-normality
## see equation 14.7 in Gelman 2013
## for t-distribution it is sqrt((scale^2)*(nu/(nu-2)))
s = coef(fit.1)['sigmaPop']
nu = 3
sa = sqrt((s^2)*(nu/(nu-2)))
# get residuals that is response minus fitted values
iResid.pa = (dfData$CD63_Act - fitted.pa) / sa
plot(fitted.pa, iResid.pa, pch=20)
lines(lowess(fitted.pa, iResid.pa))

## checking for non-linearity
plot(dfData$f423_nAraH6, iResid.pa)
lines(lowess(dfData$f423_nAraH6, iResid.pa))

plot(dfData$Avidity, iResid.pa)
lines(lowess(dfData$Avidity, iResid.pa))

## unequal variances
plot(dfData$f423_nAraH6, abs(iResid.pa))
lines(lowess(dfData$f423_nAraH6, abs(iResid.pa)))

plot(dfData$Avidity, abs(iResid.pa))
lines(lowess(dfData$Avidity, abs(iResid.pa)))

### use link and sim functions to get samples of fitted values and posterior predictive data
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
simulateOne = function(betas, sigma, nu, mModMatrix){
  f = mModMatrix %*% betas
  yrep = rt_ls(length(f), nu, f,  sigma)
  return(yrep)
}
mDraws = (sim(fit.1, data=as.list(dfData), n=200))
mMuSim = (link(fit.1, data=as.list(dfData), n=200))
dim(mDraws)
dim(mMuSim)
mParam = extract.samples(fit.1, n=200)
dim(mParam)

mDraws.2 = matrix(NA, nrow = 200, ncol = 36)
for (i in 1:200){
  mDraws.2[i,] = simulateOne(t(mParam[i,c(1,2,3)]), mParam[i,4], 3, model.matrix(CD63_Act ~ f423_nAraH6 + Avidity, data=dfData)) 
}

## visual checks
plot(density(dfData$CD63_Act))
lines(density(colMeans(mDraws)))
apply(mDraws, 1, function(x) lines(density(x), lwd=0.5, col='grey'))

plot(density(dfData$CD63_Act))
lines(density(colMeans(mDraws.2)))
apply(mDraws.2, 1, function(x) lines(density(x), lwd=0.5, col='grey'))

plot(density(colMeans(mDraws)), xlim=c(-40, 40))
lines(density(colMeans(mDraws.2)))
apply(mDraws, 1, function(x) lines(density(x), lwd=0.5, col='grey'))
apply(mDraws.2, 1, function(x) lines(density(x), lwd=0.5, col='red'))

#### plot covariates vs actual data and fitted values
plot(dfData$f423_nAraH6, dfData$CD63_Act, pch=20,
     xlab='log f423_nAraH6', ylab='%CD63 Act',
     main='Relationship of Predictor to Response')
lines(lowess(dfData$f423_nAraH6, colMeans(mMuSim)))
points(dfData$f423_nAraH6, colMeans(mMuSim), col=2, pch=20)
lines(lowess(dfData$f423_nAraH6, dfData$CD63_Act))

i = range(dfData$f423_nAraH6)
iGrid = seq(i[1], i[2], length.out = 50)
## hold second variable at average
coef(fit.1)
mFitted = link(fit.1, data=list(f423_nAraH6=iGrid, Avidity=rep(0, times=50)))
## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

### repeat for second variable
plot(dfData$Avidity, dfData$CD63_Act, pch=20,
     xlab='Avidity', ylab='%CD63 Act',
     main='Relationship of Predictor to Response')
lines(lowess(dfData$Avidity, colMeans(mMuSim)))
points(dfData$Avidity, colMeans(mMuSim), col=2, pch=20)
lines(lowess(dfData$Avidity, dfData$CD63_Act))

i = range(dfData$Avidity)
iGrid = seq(i[1], i[2], length.out = 50)
## hold second variable at average
coef(fit.1)
mFitted = link(fit.1, data=list(f423_nAraH6=rep(0, times=50), Avidity=iGrid))
## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

######################################################################################

########################################################################################
##### peanut specific ige, peanut shannon diversity and avidity
########################################################################################
## see pages 105 to 109 (rethinking book) for extract.samples, link and sim for 
## posterior predictive checks 
str(dfData)
dfData$CD63_Act = dfData$CD63.Act - mean(dfData$CD63.Act)
dfData$Peanut_Shannon = dfData$Peanut.Shannon - mean(dfData$Peanut.Shannon)

fit.lm = lm(CD63_Act ~ f13_Peanut + Peanut_Shannon + Avidity, data=dfData)
summary(fit.lm)

model.1 <- alist(
  # https://github.com/rmcelreath/rethinking/blob/ae3fae963e244c58bbd68dd86478840a87d29f49/man/dstudent.Rd#L37
  CD63_Act ~ dstudent(3, mu, sigmaPop),
  mu <- b0 + b_f13_Peanut * f13_Peanut +
    b_Peanut_Shannon * Peanut_Shannon +
    b_Avidity * Avidity,
  b0 ~ dnorm(0, 2),
  c(b_f13_Peanut, b_Peanut_Shannon, b_Avidity) ~ dnorm(0, 3),
  sigmaPop ~ dexp(1/28)
  # see here for example
  # https://github.com/rmcelreath/rethinking/blob/ae3fae963e244c58bbd68dd86478840a87d29f49/tests/rethinking_tests/test_ulam1.R#L20
  #nu ~ gamma(3, 0.1)
)

fit.1 <- quap(model.1,
              data=dfData,
              start=list(b0=0)
)

summary(fit.1)
plot(coeftab(fit.1), pars=c('b0', 'b_f13_Peanut', 'b_Peanut_Shannon', 'b_Avidity'))

m = extract.samples(fit.1, n=1000)
names(m)
pairs(m, pch=20)
round(apply(m, 2, getDifference),3)

### residual checks
fitted.pa = link(fit.1, data = as.list(dfData), n = 200)
dim(fitted.pa)
fitted.pa = colMeans(fitted.pa)

## calculate standardized residuals
## these are useful to detect non-normality
## see equation 14.7 in Gelman 2013
## for t-distribution it is sqrt((scale^2)*(nu/(nu-2)))
s = coef(fit.1)['sigmaPop']
nu = 3
sa = sqrt((s^2)*(nu/(nu-2)))
# get residuals that is response minus fitted values
iResid.pa = (dfData$CD63_Act - fitted.pa) / sa
plot(fitted.pa, iResid.pa, pch=20)
lines(lowess(fitted.pa, iResid.pa))

## checking for non-linearity
plot(dfData$f13_Peanut, iResid.pa)
lines(lowess(dfData$f13_Peanut, iResid.pa))

plot(dfData$Peanut_Shannon, iResid.pa)
lines(lowess(dfData$Peanut_Shannon, iResid.pa))

plot(dfData$Avidity, iResid.pa)
lines(lowess(dfData$Avidity, iResid.pa))

## unequal variances
plot(dfData$f13_Peanut, abs(iResid.pa))
lines(lowess(dfData$f13_Peanut, abs(iResid.pa)))

plot(dfData$Peanut_Shannon, abs(iResid.pa))
lines(lowess(dfData$Peanut_Shannon, abs(iResid.pa)))

plot(dfData$Avidity, abs(iResid.pa))
lines(lowess(dfData$Avidity, abs(iResid.pa)))

### use link and sim functions to get samples of fitted values and posterior predictive data
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu
simulateOne = function(betas, sigma, nu, mModMatrix){
  f = mModMatrix %*% betas
  yrep = rt_ls(length(f), nu, f,  sigma)
  return(yrep)
}
mDraws = (sim(fit.1, data=as.list(dfData), n=200))
mMuSim = (link(fit.1, data=as.list(dfData), n=200))
dim(mDraws)
dim(mMuSim)
mParam = extract.samples(fit.1, n=200)
dim(mParam)

mDraws.2 = matrix(NA, nrow = 200, ncol = 36)
for (i in 1:200){
  mDraws.2[i,] = simulateOne(t(mParam[i,c(1,2,3,4)]), mParam[i,5], 3, model.matrix(CD63_Act ~ f13_Peanut + Peanut_Shannon + Avidity, data=dfData)) 
}

## visual checks
plot(density(dfData$CD63_Act))
lines(density(colMeans(mDraws)))
apply(mDraws, 1, function(x) lines(density(x), lwd=0.5, col='grey'))

plot(density(dfData$CD63_Act))
lines(density(colMeans(mDraws.2)))
apply(mDraws.2, 1, function(x) lines(density(x), lwd=0.5, col='grey'))

plot(density(colMeans(mDraws)), xlim=c(-40, 40))
lines(density(colMeans(mDraws.2)))
apply(mDraws, 1, function(x) lines(density(x), lwd=0.5, col='grey'))
apply(mDraws.2, 1, function(x) lines(density(x), lwd=0.5, col='red'))

#### plot covariates vs actual data and fitted values
plot(dfData$f13_Peanut, dfData$CD63_Act, pch=20,
     xlab='log f13_Peanut', ylab='%CD63 Act',
     main='Relationship of Predictor to Response')
lines(lowess(dfData$f13_Peanut, colMeans(mMuSim)))
points(dfData$f13_Peanut, colMeans(mMuSim), col=2, pch=20)
lines(lowess(dfData$f13_Peanut, dfData$CD63_Act))

i = range(dfData$f13_Peanut)
iGrid = seq(i[1], i[2], length.out = 50)
## hold second variable at average
coef(fit.1)
mFitted = link(fit.1, data=list(f13_Peanut=iGrid, 
                                Peanut_Shannon = rep(0, times=50),
                                Avidity=rep(0, times=50)))
## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

#### plot covariates vs actual data and fitted values
## for second variable
plot(dfData$Peanut_Shannon, dfData$CD63_Act, pch=20,
     xlab='Peanut_Shannon', ylab='%CD63 Act',
     main='Relationship of Predictor to Response')
lines(lowess(dfData$Peanut_Shannon, colMeans(mMuSim)))
points(dfData$Peanut_Shannon, colMeans(mMuSim), col=2, pch=20)
lines(lowess(dfData$Peanut_Shannon, dfData$CD63_Act))

i = range(dfData$Peanut_Shannon)
iGrid = seq(i[1], i[2], length.out = 50)
## hold second variable at average
coef(fit.1)
mFitted = link(fit.1, data=list(f13_Peanut=rep(0, times=50), 
                                Peanut_Shannon = iGrid,
                                Avidity=rep(0, times=50)))
## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

### repeat for third variable
plot(dfData$Avidity, dfData$CD63_Act, pch=20,
     xlab='Avidity', ylab='%CD63 Act',
     main='Relationship of Predictor to Response')
lines(lowess(dfData$Avidity, colMeans(mMuSim)))
points(dfData$Avidity, colMeans(mMuSim), col=2, pch=20)
lines(lowess(dfData$Avidity, dfData$CD63_Act))

i = range(dfData$Avidity)
iGrid = seq(i[1], i[2], length.out = 50)
## hold second variable at average
coef(fit.1)
mFitted = link(fit.1, data=list(f13_Peanut=rep(0, times=50), 
                                Peanut_Shannon = rep(0, times=50),
                                Avidity=iGrid))
## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

######################################################################################
