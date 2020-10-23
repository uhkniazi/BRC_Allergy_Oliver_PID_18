# Name: 01_modelAndCV.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 13/10/2020
# Desc: model fitting with cross validation and roc curves

############ data loading and cleaning
source('header.R')

## load the training set
df = read.csv(file.choose(), header=T)
dim(df)

# remove NA patients
f = is.na(df$Allergic.Status)
table(f)
df = df[!f,]
dim(df)
# remove any additi
df = na.omit(df)
dim(df)
## remove white space
colnames(df)

df$Patient.ID = factor(gsub(' ', '', as.character(df$Patient.ID)))
df$Allergic.Status = factor(gsub(' ', '', as.character(df$Allergic.Status)), levels = c('PS', 'PA'))
df = droplevels.data.frame(df)
## make count matrix
mData = as.matrix(df[,-c(1:2)])
dfSample = df[,1:2]
rownames(mData) = as.character(dfSample$Patient.ID)
str(dfSample)
colnames(mData) = 'Arah2'
dim(na.omit(mData))
dim(mData)
lData.train = list(data=mData, covariates=dfSample)
rm(df)

############ end data loading

####################### cross validation under LDA and binomial models
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

## setup input data
dfData.train = data.frame(lData.train$data)
fGroups.train = lData.train$covariates$Allergic.Status

dim(dfData.train)

# create the cross validation object
url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/bernoulli.stan'
download(url, 'bernoulli.stan')

oCV.s = CCrossValidation.StanBern(train.dat = dfData.train, 
                                  test.dat = dfData.train, 
                                  test.groups = fGroups.train, 
                                  train.groups = fGroups.train,
                                  level.predict = 'PA',
                                  boot.num = 10, k.fold = 10, 
                                  ncores = 2, nchains = 2) 


plot.cv.performance(oCV.s)
unlink('bernoulli.rds')
unlink('bernoulli.stan')
