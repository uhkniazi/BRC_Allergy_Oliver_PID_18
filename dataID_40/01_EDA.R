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







