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

