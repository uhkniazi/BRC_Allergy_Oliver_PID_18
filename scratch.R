# scratch.R
# some random bits of rough work

### generate shannon diversity
df2 = read.csv('dataExternal/Mechanistic Data for Umar May 2019 OH .csv', header=T, stringsAsFactors = F)
df$Patient = as.character(df$Patient)
identical(df2$Patient.ID, df$Patient)

colnames(df)
m = as.matrix(df[,-c(1,2)])
colnames(m)
dim(m)

s = apply(m, 1, shannon)
s
hist(s)

df2$Whole.ISAC.Diversity..Shannons. = s
colnames(m)
i = grep('Ara', colnames(m))
colnames(m)[i]
s = apply(m[,i], 1, shannon)
s
hist(s)
hist(log(s+1e-5))
df2$Peanut.ISAC.Diversity..Shannon. = s
write.csv(df2, file='dataExternal/Mechanistic Data for Umar May 2019 OH with diversity.csv')
