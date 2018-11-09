# File: 02_EDA.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the data
# Date: 05/11/2018

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
mData = t(mData)
mData.bk = mData
str(dfSample)
## remove NAs and 0s by adding a jitter
table(mData < 0.3)
mData[mData < 0.3] = runif(8619, 1e-3, 1e-1)
dim(na.omit(mData))
dim(mData)

summary(mData)
# # drop the samples where average across rows is less than 1
# i = rowMeans(mData)
# table( i < 1)
# mData = mData[!(i< 1),]
# dim(mData)

ivProb = apply(mData, 1, function(inData) {
  inData[inData < 1] = 0  
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

hist(ivProb)

table(ivProb < 0.1)
# mData = mData[ivProb > 0.1, ]
dim(mData)
## load CDiagnostics and test
## compare the normalised and raw data
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(mData, 'Original')
oDiag.2 = CDiagnosticPlots(log(mData), 'Log')
# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
str(dfSample)
fBatch = dfSample$Allergic.Status
levels(fBatch)
#pdf('results/matrixClustering.pdf')

## check using various plots
par(mfrow=c(1,2))
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7, legend.pos = 'topright')
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7, legend.pos = 'topright')

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)
plot.missing.summary(oDiag.2, fBatch, axis.label.cex = 0.7, cex.main=1)

plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.PCA(oDiag.2, fBatch, cex.main=1)
# 
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8, cex.main=0.7)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)
par(mfrow=c(1,1))
plot.PCA(oDiag.1, fBatch)
plot.PCA(oDiag.2, fBatch)

plot.dendogram(oDiag.1, fBatch, labels_cex = 0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.7)

plot(oDiag.1@lData$PCA$sdev)
###############################################################################################
################ clusters and covariates explaining variance
mPC = oDiag.1@lData$PCA$x[,1:4]
i = which(mPC[,1] > 0)
fNewBatch = rep(1, times=length(dfSample$Patient))
fNewBatch[i] = 2
fNewBatch = factor(fNewBatch)

## check if batch assigned correctly
plot.PCA(oDiag.1, fNewBatch)
## try a linear mixed effect model to account for varince
library(lme4)

dfData = data.frame(mPC)
dfData = stack(dfData)
dfData$fBatch = factor(dfSample$Allergic.Status)
dfData$fAdjust1 = factor(dfSample$Patient)
dfData$fAdjust2 = fNewBatch

dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData$Coef.adj1 = dfData$fAdjust1 #factor(dfData$fAdjust1:dfData$ind)
dfData$Coef.adj2 = factor(dfData$fAdjust2:dfData$ind)

str(dfData)

fit.lme1 = lmer(values ~ 1  + (1 | Coef) + (1 | Coef.adj1) + (1 | Coef.adj2), data=dfData)
summary(fit.lme1)

fit.lme2 = lmer(values ~ 1  + (1 | Coef) + (1 | Coef.adj1), data=dfData)
summary(fit.lme2)

fit.lme3 = lmer(values ~ 1  + (1 | Coef) + (1 | Coef.adj2), data=dfData)
summary(fit.lme3)

anova(fit.lme1, fit.lme2, fit.lme3)

par(p.old)
plot((fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme1)), resid(fit.lme1)), col=2)
hist(dfData$values, prob=T)
lines(density(fitted(fit.lme1)))

plot((fitted(fit.lme2)), resid(fit.lme2), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme2)), resid(fit.lme2)), col=2)
hist(dfData$values, prob=T)
lines(density(fitted(fit.lme2)))
lines(density(fitted(fit.lme1)), col=2)

plot((fitted(fit.lme3)), resid(fit.lme3), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme3)), resid(fit.lme3)), col=2)
hist(dfData$values, prob=T)
lines(density(fitted(fit.lme2)))
lines(density(fitted(fit.lme1)), col=2)
lines(density(fitted(fit.lme3)), col=3)


library(flexmix)

fit.flex = flexmix(values ~ 1 + Coef, data=dfData, k=2)
summary(fit.flex)

fit.flex.2 = flexmix(values ~ 1 + Coef.adj2, data=dfData, k=2)
summary(fit.flex.2)

lines(density(predict(fit.flex, aggregate=T)[[1]][,1]))
lines(density(predict(fit.flex.2, aggregate=T)[[1]][,1]))

### heatmap
if (!require(NMF)) stop('R package NMF needs to be installed.')
fBatch = dfSample$Allergic.Status
i = order(fBatch)
mData = mData[, i]
fBatch = fBatch[i]

## remove NAs and 0s and convert other values to 1 by adding a jitter
f = mData < 0.3
table(f)
mData[f] = 0
mData[!f] = 1

# standardize the variables
s = apply(mData, 1, sd)
## remove any variables with sd 0
f = s <= 0
table(f)
# s = s[!f]
# mData = mData[!f,]
# mData = t(scale(t(mData)))
# cluster the samples
hc = hclust(dist(t(mData), method='binary'))
# cluster the variables
hcV = hclust(dist(mData, method = 'binary'))
colnames(mData) = as.character(fBatch)
# draw the heatmap  
pdf('results/figures/heatmaps.pdf')
aheatmap(mData, color=c('grey', 'black'), breaks=0.5, scale='none', Rowv = hcV, annColors=NA, Colv=hc)
aheatmap(mData, color=c('grey', 'black'), breaks=0.5, scale='none', Rowv = hcV, annColors=NA, Colv=NA)
dev.off(dev.cur())
