# Name: 03_variableSelection.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 6/11/2018
# Desc: try some variable selection and dimension reduction workflows

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
df$Allergic.Status = factor(gsub(' ', '', as.character(df$Allergic.Status)))
df = droplevels.data.frame(df)
## make count matrix
mData = as.matrix(df[,-c(1:2)])
dfSample = df[,1:2]
rownames(mData) = as.character(dfSample$Patient)
str(dfSample)
## remove NAs and 0s by adding a jitter
table(mData == 0)
mData[mData == 0] = round(runif(8436, 1e-3, 0.02), 3)
dim(na.omit(mData))
dim(mData)

lData.train = list(data=mData, covariates=dfSample)
rm(df)
############ end data loading

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

########################## perform a random forest step
dfData = data.frame(lData.train$data)
fGroups = lData.train$covariates$Allergic.Status

oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100)
save(oVar.r, file='temp/oVar.r_reverse.rds')

plot.var.selection(oVar.r)


######################## Stan section for binomial regression approach
dfData = data.frame(lData.train$data)
dim(dfData)
dfData$fGroups = fGroups

lData = list(resp=ifelse(dfData$fGroups == 'PS', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('tau', 'betas2'), init=initf, 
                    control=list(adapt_delta=0.99, max_treedepth = 13))

save(fit.stan, file='temp/fit.stan.binom.rds')

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
# ## get the intercept at population level
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]

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
names(m) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
text(colMeans(mCoef), ivPval, names(m), pos=1)
m = abs(m)
m = sort(m, decreasing = T)
cvTopGenes.binomial = names(m)[1:8]

p.old = par(mar=c(6,3,4,2)+0.1)
l2 = barplot(m[1:20], 
             las=2, xaxt='n', col='grey', main='Top Variables')
axis(1, at = l2, labels = names(m)[1:20], tick = F, las=2, cex.axis=0.7 )

#################### test performance of both results, from Random Forest and Binomial Regression
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# select the top 30 variables
cvTopGenes = rownames(dfRF)[1:8]

# use the top genes to find top combinations of genes
dfData = data.frame(lData.train$data[, cvTopGenes])

oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

# use the top binomial genes
dfData = data.frame(lData.train$data[, cvTopGenes.binomial])

oVar.sub2 = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub2)

# print variable combinations
for (i in 1:7){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

for (i in 1:7){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub2, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

### try a combination of top genes from both models
cvTopGenes.comb = NULL;
for (i in 1:8){
  cvTopGenes.comb = append(cvTopGenes.comb, CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)) 
  cat(i)
}
cvTopGenes.comb = unique(cvTopGenes.comb)

for (i in 1:8){
  cvTopGenes.comb = append(cvTopGenes.comb, CVariableSelection.ReduceModel.getMinModel(oVar.sub2, i))
  cat(i)
}

cvTopGenes.comb = unique(cvTopGenes.comb)
length(cvTopGenes.comb)

# use these combined variables to find top combinations of genes
dfData = data.frame(lData.train$data[,cvTopGenes.comb])

oVar.subComb = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.subComb)

for (i in 1:10){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.subComb, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

lModelReduction = list(rf=oVar.sub, bin=oVar.sub2, combined=oVar.subComb)
save(lModelReduction, file='temp/lModelReduction.rds')

################################################# binomial regression with mixture model section
################ fit a binomial model on the chosen model size based on previous results
## this can be another classifier as well e.g. LDA. Using this model check how is the performance 
## and using this make some calibration curves to select decision boundary

#cvTopGenes.comb = CVariableSelection.ReduceModel.getMinModel(lModelReduction$rf, 2)
cvTopGenes.comb = rownames(CVariableSelection.RandomForest.getVariables(oVar.r))[1:2]

library(LearnBayes)
logit.inv = function(p) {exp(p)/(exp(p)+1) }
dfData = data.frame(lData.train$data[, cvTopGenes.comb])
colnames(dfData) = cvTopGenes.comb
## binomial prediction
mypred = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  mModMatrix = data$mModMatrix
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  return(iFitted)
}

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
  iFitted = logit.inv(iFitted)
  # write the priors and likelihood
  lp = dnorm(betas[1], 0, 10, log=T) + sum(dnorm(betas[-1], 0, 10, log=T))
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lik + lp
  return(val)
}

# dfData = data.frame(dfData[ , cvTopGenes.comb])
dim(dfData)
head(dfData)
dfData = data.frame(dfData, fGroups=fGroups)

lData = list(resp=ifelse(dfData$fGroups == 'PS', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))
start = c(rep(0, times=ncol(lData$mModMatrix)))
mylogpost(start, lData)

fit.2 = laplace(mylogpost, start, lData)
fit.2

fit.1 = glm(fGroups ~ ., data=dfData, family='binomial')
data.frame(coef(fit.1), fit.2$mode)

stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('tau', 'betas2'), init=initf, 
                    control=list(adapt_delta=0.99, max_treedepth = 13))

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')
traceplot(fit.stan, 'betas2')
## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', cvTopGenes.comb)
pairs(mCoef, pch=20)


### once we have results from the classifier we can make some plots to see
### the performance
library(lattice)
library(car)
## get the predicted values
dfData.new = dfData
str(dfData.new)
## create model matrix
X = as.matrix(cbind(rep(1, times=nrow(dfData.new)), dfData.new[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict = mypred(colMeans(mCoef), list(mModMatrix=X))[,1]
xyplot(ivPredict ~ fGroups, xlab='Actual Group', ylab='Predicted Probability of Being PS (1)')
xyplot(ivPredict ~ lData.train$covariates$Allergic.Status, xlab='Actual Group', ylab='Predicted Probability of Being PS (1)',
       main='Predicted scores vs Actual groups')
densityplot(~ ivPredict, data=dfData, type='n')
densityplot(~ ivPredict | fGroups, data=dfData, type='n', xlab='Predicted Score', main='Actual Scale')
densityplot(~ ivPredict, groups=fGroups, data=dfData, type='n', 
            xlab='Predicted Score', main='Actual Scale', auto.key = list(columns=2))

## lets check on a different scale of the score
densityplot(~ logit(ivPredict), data=dfData)
xyplot(logit(ivPredict) ~ lData.train$covariates$Allergic.Status, xlab='Actual Group', ylab='Predicted Probability of Being PS (1)')
densityplot(~ logit(ivPredict), groups=fGroups, data=dfData, type='n', 
            xlab='Predicted Score', main='Logit Scale', auto.key = list(columns=2))

# convert to logit scale for model fitting
ivPredict = logit(ivPredict)