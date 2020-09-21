# scratch.R
# some random bits of rough work

### generate shannon diversity
source('header.R')
df = read.csv(file.choose(), header=T, stringsAsFactors = F, row.names = 1)
df2 = df[-1,]
m = as.matrix(df2)
dim(m)
colnames(m)
rownames(m)
m = apply(m, 2, as.numeric)
dim(m)

s = apply(m, 2, shannon)
s
hist(s)
dfResult = data.frame(Whole.ISAC.Shannon.Diversity=s)
rownames(df2)
i = grep('Ara', rownames(df2), ignore.case = T)
rownames(df2)[i]
s = apply(m[i,], 2, shannon)
s
hist(s)
hist(log(s+1e-5))
dfResult$Peanut.ISAC.Shannon.Diversity = s
identical(rownames(dfResult), colnames(df))
dfResult$Allergic.Status = as.character(df[1,])

boxplot(dfResult$Whole.ISAC.Shannon.Diversity ~ factor(dfResult$Allergic.Status))
boxplot(dfResult$Peanut.ISAC.Shannon.Diversity ~ factor(dfResult$Allergic.Status))

write.csv(dfResult, file='dataExternal/External Validation ISAC Data to Diversity.xls')
