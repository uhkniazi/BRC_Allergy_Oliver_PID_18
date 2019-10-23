# File: 04_mechanistic.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: some tests for mechanistic DAGs
# Date: 14/08/2019


source('header.R')
setwd('dataID_40/')

load('temp/dfData.pa.rds')
str(dfData.pa)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#devtools::install_github("rmcelreath/rethinking", ref = "Experimental" )
library(rethinking)
########################################################
###### univariate models
########################################################
## choose the variable to model
colnames(dfData.pa)
cVar = 'total.IgE'
##### fit a lm
fit.1.pa = lm(CD63.Act ~ ., data=dfData.pa[, c(cVar, 'CD63.Act')])
summary(fit.1.pa)

# fit.1 <- quap(
#   alist(
#     CD63.Act ~ dnorm(mu, sigmaPop),
#     mu <- b0 + b1*Ara.h.2.Sp.Act,
#     b0 ~ dcauchy(0, 2),
#     b1 ~ dnorm(0, 10),
#     sigmaPop ~ dexp(1)
#   ), data=dfData.pa
# )
# summary(fit.1)

stanDso = rstan::stan_model(file='tResponseRegression_censored.stan')

## model dataset PA
m = model.matrix(CD63.Act ~ ., data=dfData.pa[dfData.pa$CD63.Act > 0.01, c(cVar, 'CD63.Act')])
m2 = model.matrix(CD63.Act ~ ., data=dfData.pa[dfData.pa$CD63.Act <= 0.01, c(cVar, 'CD63.Act')])

lStanData = list(Ntotal=nrow(m), Ncol=ncol(m), X=m, rLower=0.01,
                 X2 = m2,
                 Ncens = nrow(m2),
                 y=dfData.pa$CD63.Act[dfData.pa$CD63.Act > 0.01])

fit.stan.pa = sampling(stanDso, data=lStanData, iter=2000, chains=2, 
                       pars=c('betas', 'mu', 'y_cens', 'mu2', 'sigmaPop', 'nu', 'log_lik'),
                       cores=2)

################# model checks
print(fit.stan.pa, c('betas', 'sigmaPop', 'nu'), digits=3)

traceplot(fit.stan.pa, c('betas'))

WAIC(fit.stan.pa)
LOOPk(fit.stan.pa)
plot(LOOPk(fit.stan.pa) ~ WAIC(fit.stan.pa, pointwise = T))

## testing some summary functions from rethinking package
post = extract.samples(fit.stan.pa, n=1000)
dim(post)
length(post)
names(post)
post = as.data.frame(post$betas)
head(post)
precis(post)

betas.pa = extract(fit.stan.pa, 'betas')$betas
colnames(betas.pa) = c('b0', 'b1')
dim(betas.pa)
precis(as.data.frame(betas.pa))

plot(CD63.Act ~ dfData.pa[,cVar], data=dfData.pa, xlab=cVar)
curve(mean(betas.pa[,'b0']) + mean(betas.pa[,'b1']) * x, add=T)
abline(a = mean(betas.pa[,'b0']), b= mean(betas.pa[,'b1']), col='red')

i = range(dfData.pa[,cVar])
iGrid = seq(i[1], i[2], length.out = 50)
mGrid = cbind(1, iGrid)
mFitted = t(mGrid %*% t(betas.pa))
dim(mFitted)

## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

############### new simulated data
###############
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

## sample n values, 1000 times
mDraws.sim = matrix(NA, nrow = length(iGrid), ncol=300)
l = extract(fit.stan.pa)
for (i in 1:300){
  p = sample(1:nrow(l$betas), 1)
  mDraws.sim[,i] = simulateOne(l$betas[p,], 
                               l$sigmaPop[p],
                               l$nu[p],
                               mGrid)
}

dim(mDraws.sim)
mDraws.PI = apply(mDraws.sim, 1, PI)
shade(mDraws.PI, iGrid)
apply(mDraws.sim, 2, function(x) lines(lowess(iGrid, x), lwd=0.5, col='red'))
#apply(mDraws.sim, 2, function(x) points(iGrid, x, cex=0.5, pch=20, col='red'))
###############

# compare with lm 
# coef(fit.1.pa); colMeans(betas.pa); coef(fit.1)
# 
# s2 = cbind(extract(fit.stan.pa)$betas)
# colnames(s2) = c(colnames(lStanData$X))
# dim(s2)
s2 = betas.pa
pairs(s2[sample(1:nrow(s2), 1000),], pch=20, col='grey', cex=0.5, main='PA')

########################################################
# dfData.pa.bk = dfData.pa
# dfData.pa = dfData.pa[dfData.pa$CD63.Act > 0.01,]
# 
# ######## various model sizes
# fit.1 <- quap(
#   alist(
#     CD63.Act ~ dnorm(mu, sigmaPop),
#     mu <- b0 + b1*f13.Peanut,
#     b0 ~ dcauchy(0, 2),
#     b1 ~ dnorm(0, 10),
#     sigmaPop ~ dexp(1)
#   ), data=dfData.pa,
#   start=list(b0=21)
# )
# summary(fit.1)
# 
# fit.2 <- quap(
#   alist(
#     CD63.Act ~ dnorm(mu, sigmaPop),
#     mu <- b0 + b2*ISAC.Shannon,
#     b0 ~ dcauchy(0, 2),
#     b2 ~ dnorm(0, 10),
#     sigmaPop ~ dexp(1)
#   ), data=dfData.pa,
#   start=list(b0=21)
# )
# summary(fit.2)
# 
# fit.3 <- quap(
#   alist(
#     CD63.Act ~ dnorm(mu, sigmaPop),
#     mu <- b0 + b3*Peanut.Sp.Act,
#     b0 ~ dcauchy(0, 2),
#     b3 ~ dnorm(0, 10),
#     sigmaPop ~ dexp(1)
#   ), data=dfData.pa,
#   start=list(b0=21)
# )
# summary(fit.3)
# 
# fit.4 <- quap(
#   alist(
#     CD63.Act ~ dnorm(mu, sigmaPop),
#     mu <- b0 + b1*f13.Peanut + b2*ISAC.Shannon + b3*Peanut.Sp.Act,
#     b0 ~ dcauchy(0, 2),
#     c(b1, b2, b3) ~ dnorm(0, 10),
#     sigmaPop ~ dunif(1, 40)
#   ), data=dfData.pa,
#   start=list(b0=21, b1=0, b2=0, b3=0)
# )
# summary(fit.4)
# 
# fit.5 <- quap(
#   alist(
#     CD63.Act ~ dnorm(mu, sigmaPop),
#     mu <- b0 + b1*Ara.h.2.Sp.Act + b2*ISAC.Shannon, #+ b3*Peanut.Sp.Act,
#     b0 ~ dcauchy(0, 2),
#     c(b1, b2) ~ dnorm(0, 10),
#     sigmaPop ~ dunif(1, 40)
#   ), data=dfData.pa,
#   start=list(b0=21, b1=0, b2=0)
# )
# summary(fit.5)
# 
# fit.6 <- quap(
#   alist(
#     CD63.Act ~ dnorm(mu, sigmaPop),
#     mu <- b0 + b1*Ara.h.2.Sp.Act + b3*Peanut.Sp.Act,
#     b0 ~ dcauchy(0, 2),
#     c(b1, b3) ~ dnorm(0, 10),
#     sigmaPop ~ dunif(1, 40)
#   ), data=dfData.pa,
#   start=list(b0=21, b1=0, b3=0)
# )
# summary(fit.6)
# 
# fit.7 <- quap(
#   alist(
#     CD63.Act ~ dnorm(mu, sigmaPop),
#     mu <- b0 + b2*ISAC.Shannon + b3*Peanut.Sp.Act,
#     b0 ~ dcauchy(0, 2),
#     c(b2, b3) ~ dnorm(0, 10),
#     sigmaPop ~ dunif(1, 40)
#   ), data=dfData.pa,
#   start=list(b0=21, b2=0, b3=0)
# )
# summary(fit.7)
# 
# 
# plot(coeftab(fit.1, fit.2, fit.3, fit.4), pars=c('b1', 'b2', 'b3'))
# 
# cov2cor(vcov(fit.4))
# sapply(list(fit.1, fit.2, fit.3, fit.4), function(m) sum(lppd(m)))
# compare(fit.1, fit.2, fit.3, fit.4)
# dfData.pa = dfData.pa.bk

###################################
####### multivariate censored model
## model dataset PA
cVar = c(cVar, 'f13.Peanut')
m = model.matrix(CD63.Act ~ ., data=dfData.pa[dfData.pa$CD63.Act > 0.01, c(cVar, 'CD63.Act')])
m2 = model.matrix(CD63.Act ~ ., data=dfData.pa[dfData.pa$CD63.Act <= 0.01, c(cVar, 'CD63.Act')])

lStanData = list(Ntotal=nrow(m), Ncol=ncol(m), X=m, rLower=0.01,
                 X2 = m2,
                 Ncens = nrow(m2),
                 y=dfData.pa$CD63.Act[dfData.pa$CD63.Act > 0.01])

fit.stan.2 = sampling(stanDso, data=lStanData, iter=2000, chains=2, 
                      pars=c('betas', 'mu', 'y_cens', 'mu2', 'sigmaPop', 'nu', 'log_lik'),
                       cores=2)

################# model checks
print(fit.stan.2, c('betas', 'sigmaPop', 'nu'), digits=3)

compare(fit.stan.pa, fit.stan.2)
plot(compare(fit.stan.pa, fit.stan.2))
compare(fit.stan.pa, fit.stan.2, func = LOO)

ct = coeftab(fit.stan.pa, fit.stan.2)
rn = rownames(ct@coefs)
i = grep('betas', rn)
i = i[-1]
rownames(ct@coefs)[i] = cVar
rownames(ct@se)[i] = cVar
plot(ct, pars=cVar)

# traceplot(fit.stan.pa, c('betas'))
# 
betas.pa = extract(fit.stan.2)$betas
colnames(betas.pa) = c('b0', cVar)
# dim(betas.pa)
# precis(as.data.frame(betas.pa))

## set index of cVar to 1 or 2 depending on which variable is being modelled
## the other can be zero as the data has been centered
plot(CD63.Act ~ dfData.pa[,cVar[2]], data=dfData.pa, xlab=cVar[2])
curve(mean(betas.pa[,'b0']) + mean(betas.pa[,cVar[2]]) * x, add=T)
abline(a = mean(betas.pa[,'b0']), b= mean(betas.pa[,cVar[2]]), col='red')

i = range(dfData.pa[,cVar[2]])
iGrid = seq(i[1], i[2], length.out = 50)
mGrid = cbind(1, 0, iGrid)
dim(mGrid)
mFitted = t(mGrid %*% t(betas.pa))
dim(mFitted)

## posterior predictive values for fitted
lines(iGrid, colMeans(mFitted), col='green')
mu.hpdi = apply(mFitted, 2, HPDI)
shade(mu.hpdi, iGrid)

############### new simulated data
###############
### generate some posterior predictive data

## sample n values, 1000 times
mDraws.sim = matrix(NA, nrow = length(iGrid), ncol=300)
l = extract(fit.stan.2)
for (i in 1:300){
  p = sample(1:nrow(l$betas), 1)
  mDraws.sim[,i] = simulateOne(l$betas[p,], 
                               l$sigmaPop[p],
                               l$nu[p],
                               mGrid)
}

dim(mDraws.sim)
mDraws.PI = apply(mDraws.sim, 1, PI)
shade(mDraws.PI, iGrid)
apply(mDraws.sim, 2, function(x) lines(lowess(iGrid, x), lwd=0.5, col='red'))

### repeat for the second covariate
# plot(CD63.Act ~ dfData.pa[,cVar[2]], data=dfData.pa, xlab=cVar[2])
# curve(mean(betas.pa[,'b0']) + mean(betas.pa[,'b2']) * x, add=T)
# abline(a = mean(betas.pa[,'b0']), b= mean(betas.pa[,'b2']), col='red')
# 
# i = range(dfData.pa[,cVar[2]])
# iGrid = seq(i[1], i[2], length.out = 50)
# mGrid = cbind(1, 0, iGrid)
# mFitted = t(mGrid %*% t(betas.pa))
# dim(mFitted)
# 
# ## posterior predictive values for fitted
# lines(iGrid, colMeans(mFitted), col='green')
# mu.hpdi = apply(mFitted, 2, HPDI)
# shade(mu.hpdi, iGrid)
# 
# ############### new simulated data
# 
# ## sample n values, 1000 times
# mDraws.sim = matrix(NA, nrow = length(iGrid), ncol=300)
# l = extract(fit.stan.pa)
# for (i in 1:300){
#   p = sample(1:nrow(l$betas), 1)
#   mDraws.sim[,i] = simulateOne(l$betas[p,], 
#                                l$sigmaPop[p],
#                                l$nu[p],
#                                mGrid)
# }
# 
# dim(mDraws.sim)
# mDraws.PI = apply(mDraws.sim, 1, PI)
# shade(mDraws.PI, iGrid)
# apply(mDraws.sim, 2, function(x) lines(lowess(iGrid, x), lwd=0.5, col='red'))

s2 = cbind(extract(fit.stan.2)$betas)
colnames(s2) = c(colnames(lStanData$X))
dim(s2)
pairs(s2[sample(1:nrow(s2), 1000),], pch=20, col='grey', cex=0.5, main='PA')

## posterior predictive values vs actual data
mGrid = model.matrix(CD63.Act ~ ., data=dfData.pa[, c(cVar, 'CD63.Act')])
mFitted = t(mGrid %*% t(betas.pa))
dim(mFitted)
mu.1 = mFitted
mu.1.pi = apply(mu.1, 2, PI)
plot(colMeans(mu.1) ~ dfData.pa$CD63.Act, ylim=range(mu.1.pi))
abline(0, 1)
for (i in 1:nrow(dfData.pa)) lines(rep(dfData.pa$CD63.Act[i], 2), mu.1.pi[,i])
r = dfData.pa$CD63.Act - colMeans(mu.1)
plot(r ~ colMeans(mu.1))
lines(lowess(r ~ colMeans(mu.1)))

## how does the model see the data in terms of individual observations P 231 rethinking plot
w1 = WAIC(fit.stan.pa, pointwise = T)
w2 = WAIC(fit.stan.2, pointwise = T)
w = as.numeric(w1 - w2)
plot(lStanData$y ~ w)


##############################################################################################
############################################### check behaviour of other subset of predictors
##############################################################################################
dfData.bk = dfData.pa
dfData.pa = dfData.pa[dfData.pa$CD63.Act > 0.01,]
colnames(dfData.pa)
dfData.pa = dfData.pa[,-1]
dfData.pa = apply(dfData.pa, 2, scale)
str(dfData.pa)
dfData.pa = data.frame(dfData.pa)
pairs(dfData.pa)


fit.1 <- quap(
  alist(
    f13.Peanut ~ dnorm(mu, sigmaPop),
    mu <- b0 + b1*total.IgE,
    b0 ~ dcauchy(0, 2),
    c(b1) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData.pa,
  start=list(b0=0)
)
summary(fit.1)
precis(fit.1)
post = extract.samples(fit.1, n=1000)
dim(post)
head(post)
precis(post)
plot(precis(post))
#plot(coeftab(fit.1), pars=c('b1', 'b2', 'b3', 'b4'))
cov2cor(vcov(fit.1))
pairs(post)

fit.1.cd <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b1*f13.Peanut + b2*total.IgE,
    b0 ~ dcauchy(0, 2),
    c(b1, b2) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData.pa,
  start=list(b0=0)
)
summary(fit.1.cd)

fit.2.cd <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b1*f13.Peanut,
    b0 ~ dcauchy(0, 2),
    c(b1) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData.pa,
  start=list(b0=0)
)
summary(fit.2.cd)

fit.3.cd <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b2*total.IgE,
    b0 ~ dcauchy(0, 2),
    c(b2) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData.pa,
  start=list(b0=0)
)
summary(fit.3.cd)

plot(compare(fit.1.cd, fit.2.cd, fit.3.cd))
plot(coeftab(fit.1.cd, fit.2.cd, fit.3.cd), pars=c('b1', 'b2'))
compare(fit.1.cd, fit.2.cd, fit.3.cd)
compare(fit.1.cd, fit.2.cd, fit.3.cd, func=LOO)

# ## second model
# mu.1 = link(fit.2)
# mu.1.pi = apply(mu.1, 2, PI)
# plot(colMeans(mu.1) ~ dfData.pa$Peanut.Sp.Act, ylim=range(mu.1.pi))
# abline(0, 1)
# for (i in 1:nrow(dfData.pa)) lines(rep(dfData.pa$Peanut.Sp.Act[i], 2), mu.1.pi[,i])
# r = dfData.pa$Peanut.Sp.Act - colMeans(mu.1)
# plot(r ~ colMeans(mu.1))
# lines(lowess(r ~ colMeans(mu.1)))
# 
# ## third model
# mu.1 = link(fit.3)
# mu.1.pi = apply(mu.1, 2, PI)
# plot(colMeans(mu.1) ~ dfData.pa$Peanut.Sp.Act, ylim=range(mu.1.pi))
# abline(0, 1)
# for (i in 1:nrow(dfData.pa)) lines(rep(dfData.pa$Peanut.Sp.Act[i], 2), mu.1.pi[,i])
# r = dfData.pa$Peanut.Sp.Act - colMeans(mu.1)
# plot(r ~ colMeans(mu.1))
# lines(lowess(r ~ colMeans(mu.1)))

# ## new model with actual response vs predictor
# fit.4 <- quap(
#   alist(
#     CD63.Act ~ dnorm(mu, sigmaPop),
#     mu <- b0 + b1*f13.Peanut,
#     b0 ~ dcauchy(0, 2),
#     c(b1) ~ dnorm(0, 1),
#     sigmaPop ~ dexp(1)
#   ), data=dfData.pa,
#   start=list(b0=0)
# )
# summary(fit.4)
# plot(coeftab(fit.4))

## simulate the DAGs
# totalIge -> F13.Peanut
nsim = 1000
df.Sim = data.frame(1:nsim)
df.Sim$total.IgE = rnorm(nsim)
# df.Sim$f422.rAra.h.1 = rnorm(nsim)
# df.Sim$f423.Ara.h.2 = rnorm(nsim)
# df.Sim$f424.rAra.h.3 = rnorm(nsim)
# df.Sim$f423.nAra.H.6 = rnorm(nsim)
# s = sim(fit.1, data = df.Sim, n = 100)
# df.Sim$Peanut.Sp.Act = colMeans(s)
m = cbind(1, as.matrix(df.Sim[,2]))
head(m)
coef(fit.1)
m = m %*% coef(fit.1)[-3]
df.Sim$f13.Peanut = rnorm(nsim, m)

# fit.sim <- quap(
#   alist(
#     Peanut.Sp.Act ~ dnorm(mu, sigmaPop),
#     mu <- b0 + b1*Ara.h.2.Sp.Act + b2*Ara.h.6.Sp.Act,
#     b0 ~ dcauchy(0, 2),
#     c(b1, b2) ~ dnorm(0, 1),
#     sigmaPop ~ dexp(1)
#   ), data=df.Sim,
#   start=list(b0=0)
# )

## recover coefficients using simulated data fit
fit.1.sim <- quap(
  alist(
    f13.Peanut ~ dnorm(mu, sigmaPop),
    mu <- b0 + b1*total.IgE,
    b0 ~ dcauchy(0, 2),
    c(b1) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim,
  start=list(b0=0)
)

plot(coeftab(fit.1, fit.1.sim))
## simulate new response variable
str(df.Sim)
coef(fit.1.cd)
m = cbind(1, df.Sim$f13.Peanut, df.Sim$total.IgE)
head(m)
m = m %*% coef(fit.1.cd)[-4]
dim(m)
df.Sim$CD63.Act = rnorm(nsim, m[,1])

## now fit the simulated data in the same model
fit.1.cd.sim <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b1*f13.Peanut + b2*total.IgE,
    b0 ~ dcauchy(0, 2),
    c(b1, b2) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim,
  start=list(b0=0)
)
summary(fit.1.cd.sim)
plot(coeftab(fit.1.cd, fit.1.cd.sim), pars=c('b0', 'b1', 'b2'))

plot(density(dfData.pa$CD63.Act))
plot(density(df.Sim$CD63.Act))
# fit.sim <- quap(
#   alist(
#     CD63.Act ~ dnorm(mu, sigmaPop),
#     mu <- b0 + b1*Peanut.Sp.Act,
#     b0 ~ dcauchy(0, 2),
#     c(b1) ~ dnorm(0, 1),
#     sigmaPop ~ dexp(1)
#   ), data=df.Sim,
#   start=list(b0=0)
# )
# summary(fit.sim)
fit.4.sim <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b1*f13.Peanut,
    b0 ~ dcauchy(0, 2),
    c(b1) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim,
  start=list(b0=0)
)


plot(coeftab(fit.4, fit.4.sim), pars=c('b1'))
plot(coeftab(fit.4, fit.4.sim))
###################### add another variable to the simulation
## find the coefficient first for new variable
fit.5 <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b2*ISAC.Shannon,
    b0 ~ dcauchy(0, 2),
    c(b2) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData.pa,
  start=list(b0=0)
)
summary(fit.5)
df.Sim$ISAC.Shannon = rnorm(nsim)
df.Sim$CD63.Act = rnorm(nsim, 0.58*df.Sim$Peanut.Sp.Act + 0.38*df.Sim$ISAC.Shannon)

fit.sim <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b1*Peanut.Sp.Act + b2*ISAC.Shannon,
    b0 ~ dcauchy(0, 2),
    c(b1, b2) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=df.Sim,
  start=list(b0=0)
)
summary(fit.sim)

plot(coeftab(fit.4, fit.5, fit.sim), pars=c('b1', 'b2'))

## multivariate model on actual data
fit.6 <- quap(
  alist(
    CD63.Act ~ dnorm(mu, sigmaPop),
    mu <- b0 + b2*ISAC.Shannon + b1*Peanut.Sp.Act,
    b0 ~ dcauchy(0, 2),
    c(b1, b2) ~ dnorm(0, 1),
    sigmaPop ~ dexp(1)
  ), data=dfData.pa,
  start=list(b0=0)
)
summary(fit.6)

plot(coeftab(fit.4, fit.5, fit.sim, fit.6), pars=c('b1', 'b2'))
# fit.5 <- quap(
#   alist(
#     CD63.Act ~ dnorm(mu, sigmaPop),
#     mu <- b0 + b1*Peanut.Shannon + b2*ISAC.Shannon + b3*Peanut.Sp.Act,
#     b0 ~ dcauchy(0, 2),
#     c(b1, b2, b3) ~ dnorm(0, 10),
#     sigmaPop ~ dunif(1, 40)
#   ), data=df.Sim,
#   start=list(b0=21, b1=0, b2=0, b3=0)
# )
summary(fit.5)
dfData.pa = df.Sim
dfData.pa = dfData.pa.bk
