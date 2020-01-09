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
xyplot(y ~ log(values+0.5) | ind, data=df, groups=g, type=c('g', 'p', 'r'), pch=19, cex=0.6,
       index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'),
       auto.key=list(columns=2))

xyplot(y ~ log(values+0.5) | ind, data=df, groups=g, type=c('g', 'p', 'smooth'), pch=19, cex=0.6,
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

dfData$total.IgE = log(dfData$total.IgE)
dfData$f13.Peanut = log(dfData$f13.Peanut)

dfData$f422.rAra.h.1[dfData$f422.rAra.h.1 <= 0] = 1e-2
dfData$f422.rAra.h.1 = log(dfData$f422.rAra.h.1)

dfData$f423.Ara.h.2[dfData$f423.Ara.h.2 <= 0] = 1e-2
dfData$f423.Ara.h.2 = log(dfData$f423.Ara.h.2)

dfData$f424.rAra.h.3[dfData$f424.rAra.h.3 <= 0] = 1e-2
dfData$f424.rAra.h.3 = log(dfData$f424.rAra.h.3)

dfData$f423.nAra.H.6[dfData$f423.nAra.H.6 <= 0] = 1e-2
dfData$f423.nAra.H.6 = log(dfData$f423.nAra.H.6)

dfData$f352.rAra.h.8[dfData$f352.rAra.h.8 <= 0] = 1e-2
dfData$f352.rAra.h.8 = log(dfData$f352.rAra.h.8)

dfData$f427.rAra.h.9[dfData$f427.rAra.h.9 <= 0] = 1e-2
dfData$f427.rAra.h.9 = log(dfData$f427.rAra.h.9)

library(car) # to use logit function
dfData$Peanut.Sp.Act = logit(dfData$Peanut.Sp.Act/100)
dfData$Ara.h.2.Sp.Act = logit(dfData$Ara.h.2.Sp.Act/100)
dfData$Ara.h.6.Sp.Act = logit(dfData$Ara.h.6.Sp.Act/100)
dfData$CD63.Act[dfData$CD63.Act <= 0] = 0
dfData$CD63.Act = logit(dfData$CD63.Act)

as.data.frame(colnames(dfData))
m = as.matrix(dfData[,-c(1:2)])
pairs(m[dfData$Allergic.Status == 'PA', ], pch=20)

df = stack(dfData[,-c(1, 2, 18)])
df$y = dfData$CD63.Act
df$g = dfData$Allergic.Status

xyplot(y ~ values | ind, data=df, groups=g, type=c('g', 'p', 'r'), pch=19, cex=0.6,
       index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'),
       auto.key=list(columns=2))

xyplot(y ~ values | ind, data=df, groups=g, type=c('g', 'p', 'smooth'), pch=19, cex=0.6,
       index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free'),
       auto.key=list(columns=2))

xyplot(y ~ values | ind, data=df[df$g == 'PA',], type=c('g', 'p', 'r'), pch=19, cex=0.6,
       index.cond = function(x,y) coef(lm(y ~ x))[1], #aspect='xy',# layout=c(8,2),
       par.strip.text=list(cex=0.7), scales = list(x=list(rot=45, cex=0.5), relation='free')
       )

#save(dfData, file='dataExternal/dfData_scaled_with_avidity_dataID40.rds')
load('dataExternal/dfData_scaled_with_avidity_dataID40.rds')
setwd('dataID_40/')

str(dfData)
dfData = dfData[dfData$Allergic.Status == 'PA', ]

library(rethinking)

## scale the data 
dfData = data.frame(scale(dfData[,-c(1:3)]))

########################################################################################
#################### model 1 with all covariates
########################################################################################
fit.all <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b_Peanut.SPT * Peanut.SPT + 
      b_total.IgE * total.IgE +
      b_f13.Peanut * f13.Peanut +
      b_f422.rAra.h.1 * f422.rAra.h.1 +
      b_f423.Ara.h.2 * f423.Ara.h.2 +
      b_f424.rAra.h.3 * f424.rAra.h.3 +
      b_f423.nAra.H.6 * f423.nAra.H.6 +
      b_f352.rAra.h.8 * f352.rAra.h.8 +
      b_f427.rAra.h.9 * f427.rAra.h.9 +
      b_Peanut.Sp.Act * Peanut.Sp.Act +
      b_Ara.h.2.Sp.Act * Ara.h.2.Sp.Act +
      b_Ara.h.6.Sp.Act * Ara.h.6.Sp.Act +
      b_ISAC.Shannon * ISAC.Shannon +
      b_Peanut.Shannon + Peanut.Shannon +
      b_Avidity + Avidity,
    b0 ~ dnorm(0, 2),
    c(b_Peanut.SPT, b_total.IgE, b_f13.Peanut, b_f422.rAra.h.1, b_f423.Ara.h.2,
      b_f424.rAra.h.3, b_f423.nAra.H.6, b_f352.rAra.h.8, b_f427.rAra.h.9,
      b_Peanut.Sp.Act, b_Ara.h.2.Sp.Act, b_Ara.h.6.Sp.Act, b_ISAC.Shannon,
      b_Peanut.Shannon, b_Avidity) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData,
  start=list(b0=0)
)

summary(fit.all)
plot(coeftab(fit.all))
m = extract.samples(fit.all, n=200)
names(m)
pairs(m, pch=20)

########### simulate fake data from this model
nsim = nrow(dfData)
df.Sim = data.frame(1:nsim)
df.Sim$Peanut.SPT = rnorm(nsim)
df.Sim$total.IgE = rnorm(nsim)
df.Sim$f13.Peanut = rnorm(nsim)
df.Sim$f422.rAra.h.1 = rnorm(nsim)
df.Sim$f423.Ara.h.2 = rnorm(nsim)
df.Sim$f424.rAra.h.3 = rnorm(nsim)
df.Sim$f423.nAra.H.6 = rnorm(nsim)
df.Sim$f352.rAra.h.8 = rnorm(nsim)
df.Sim$f427.rAra.h.9 = rnorm(nsim)
df.Sim$Peanut.Sp.Act = rnorm(nsim)
df.Sim$Ara.h.2.Sp.Act = rnorm(nsim)
df.Sim$Ara.h.6.Sp.Act = rnorm(nsim)
df.Sim$ISAC.Shannon = rnorm(nsim)
df.Sim$Peanut.Shannon = rnorm(nsim)
df.Sim$Avidity = rnorm(nsim)
colnames(df.Sim)

## model matrix
m = cbind(1, as.matrix(df.Sim[,-1]))
head(m)
coef(fit.all)
m = m %*% coef(fit.all)[-17]
df.Sim$CD63.Act = rnorm(nsim, m)

## refit the model to the simulated data
fit.all.sim <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b_Peanut.SPT * Peanut.SPT + 
      b_total.IgE * total.IgE +
      b_f13.Peanut * f13.Peanut +
      b_f422.rAra.h.1 * f422.rAra.h.1 +
      b_f423.Ara.h.2 * f423.Ara.h.2 +
      b_f424.rAra.h.3 * f424.rAra.h.3 +
      b_f423.nAra.H.6 * f423.nAra.H.6 +
      b_f352.rAra.h.8 * f352.rAra.h.8 +
      b_f427.rAra.h.9 * f427.rAra.h.9 +
      b_Peanut.Sp.Act * Peanut.Sp.Act +
      b_Ara.h.2.Sp.Act * Ara.h.2.Sp.Act +
      b_Ara.h.6.Sp.Act * Ara.h.6.Sp.Act +
      b_ISAC.Shannon * ISAC.Shannon +
      b_Peanut.Shannon + Peanut.Shannon +
      b_Avidity + Avidity,
    b0 ~ dnorm(0, 2),
    c(b_Peanut.SPT, b_total.IgE, b_f13.Peanut, b_f422.rAra.h.1, b_f423.Ara.h.2,
      b_f424.rAra.h.3, b_f423.nAra.H.6, b_f352.rAra.h.8, b_f427.rAra.h.9,
      b_Peanut.Sp.Act, b_Ara.h.2.Sp.Act, b_Ara.h.6.Sp.Act, b_ISAC.Shannon,
      b_Peanut.Shannon, b_Avidity) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim,
  start=list(b0=0)
)

plot(coeftab(fit.all, fit.all.sim))

########################################################################################

########################################################################################
##### peanut specific ige and avidity
########################################################################################
fit.1 <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b_f13.Peanut * f13.Peanut +
      b_Avidity + Avidity,
    b0 ~ dnorm(0, 2),
    c(b_f13.Peanut, b_Avidity) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData,
  start=list(b0=0)
)

summary(fit.1)
plot(coeftab(fit.1))
m = extract.samples(fit.1, n=200)
names(m)
pairs(m, pch=20)

########### simulate fake data from this model
nsim = nrow(dfData)
df.Sim = data.frame(1:nsim)
df.Sim$f13.Peanut = rnorm(nsim)
df.Sim$Avidity = rnorm(nsim)
colnames(df.Sim)

## model matrix
m = cbind(1, as.matrix(df.Sim[,-1]))
head(m)
coef(fit.1)
m = m %*% coef(fit.1)[-4]
df.Sim$CD63.Act = rnorm(nsim, m)

fit.1.sim <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b_f13.Peanut * f13.Peanut +
      b_Avidity + Avidity,
    b0 ~ dnorm(0, 2),
    c(b_f13.Peanut, b_Avidity) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim,
  start=list(b0=0)
)
plot(coeftab(fit.1, fit.1.sim))

#### level 2 of the model, subsets of f13
fit.2 <- quap(
  alist(
    f13.Peanut ~ dnorm(mu, sigmaPop),
    mu <- b0 + 
      b_f422.rAra.h.1 * f422.rAra.h.1 +
      b_f423.Ara.h.2 * f423.Ara.h.2 +
      b_f424.rAra.h.3 * f424.rAra.h.3 +
      b_f423.nAra.H.6 * f423.nAra.H.6 +
      b_f352.rAra.h.8 * f352.rAra.h.8 +
      b_f427.rAra.h.9 * f427.rAra.h.9, 
    b0 ~ dnorm(0, 2),
    c(b_f422.rAra.h.1, b_f423.Ara.h.2,
      b_f424.rAra.h.3, b_f423.nAra.H.6, b_f352.rAra.h.8, b_f427.rAra.h.9
      ) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData,
  start=list(b0=0)
)
summary(fit.2)
plot(coeftab(fit.2))
pairs(extract.samples(fit.2, n=100))

## simulate the data from level 2
nsim = nrow(dfData)
df.Sim = data.frame(1:nsim)
df.Sim$f422.rAra.h.1 = rnorm(nsim)
df.Sim$f423.Ara.h.2 = rnorm(nsim)
df.Sim$f424.rAra.h.3 = rnorm(nsim)
df.Sim$f423.nAra.H.6 = rnorm(nsim)
df.Sim$f352.rAra.h.8 = rnorm(nsim)
df.Sim$f427.rAra.h.9 = rnorm(nsim)
colnames(df.Sim)

## model matrix
m = cbind(1, as.matrix(df.Sim[,-1]))
head(m)
coef(fit.2)
m = m %*% coef(fit.2)[-8]
df.Sim$f13.Peanut = rnorm(nsim, m)

## simulate level 1 of the data
df.Sim$Avidity = rnorm(nsim)
colnames(df.Sim)

## model matrix
m = cbind(1, as.matrix(df.Sim[,c('f13.Peanut', 'Avidity')]))
head(m)
coef(fit.1)
m = m %*% coef(fit.1)[-4]
df.Sim$CD63.Act = rnorm(nsim, m)

fit.1.sim <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b_f13.Peanut * f13.Peanut +
      b_Avidity + Avidity,
    b0 ~ dnorm(0, 2),
    c(b_f13.Peanut, b_Avidity) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim,
  start=list(b0=0)
)
plot(coeftab(fit.1, fit.1.sim))


########################################################################################










