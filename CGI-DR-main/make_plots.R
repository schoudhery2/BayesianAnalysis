args = commandArgs(trailingOnly=TRUE)

ifile = args[1]
gene = args[2]
data = read.table(ifile,head=T,sep='\t')

subset = data[data$gene==gene,]
print(dim(subset))
conc0 = subset[,6:8]
conc0vals = as.matrix(conc0)
hist(log10(conc0vals),breaks=20,main="depletion (abundance at conc=0 relative to -ATC)")

#####################
# noise level

#####################

#summary(data$betaE)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.71179 -0.30520 -0.06266 -0.20121 -0.01427  0.84667 

data$betaE = data$betaE-min(data$betaE)+0.01 # shift positive

subset = data[data$gene==gene,]

# squashing functions:
# 1/(1+exp(-x))
# (1-exp(-x))/(1+exp(-x))
# (1-exp(-2x))/(1+exp(-2x))
# Wolfram Alpha: plot (1-exp(-2x))/(1+exp(-2x)) from -1 to 3

for (i in 6:17)
{
  subset[,i] = (1-exp(-2*subset[,i]))/(1+exp(-2*subset[,i]))
}

conc0 = subset[,6:8]
conc0vals = as.matrix(conc0)
hist(log10(conc0vals),breaks=20,main="squashed abund0")

subset$conc0mean = apply(conc0,1,mean)

plot(subset$betaE,log10(subset$conc0mean),ylim=c(-4,1),main=gene)

###################

logsigmoid = function(x) { log10(x/(1-x)) }

hist(logsigmoid(conc0vals),breaks=20,main="log(sigmoid(abund0))")

plot(log10(subset$betaE),logsigmoid(subset$conc0mean),ylim=c(-5,5),main=gene)

###################

slopes = NULL
for (i in 1:nrow(subset))
{
  row = subset[i,]
  vals = row[6:17]
  PC = 0.01
  vals = as.numeric(lapply(vals,function(x) { max(PC,x) }))
  logconcs = c(0,0,0,1,1,1,2,2,2,3,3,3) # already on pretend log-scale
  df = data.frame(t(rbind(logconc=logconcs,abund=vals)))
  LOG = T
  if (LOG)
  {
    #mod = lm(log10(abund)~logconc,data=df)
    mod = lm(logsigmoid(abund)~logconc,data=df)
    plot(logconcs,predict(mod),type="l",ylim=c(-5,1),main=gene,ylab="predicted logsigmoid(abund)")
  }
  else
  {
    mod = lm(abund~logconc,data=df)
    plot(logconcs,predict(mod),type="l",ylim=c(0,3),main=gene)
  }
  slope = coefficients(mod)[2]
  slopes = cbind(slopes,slope)
  par(new=T)
}

par(new=F)
#plot(subset$betaE,slopes,xlim=c(-1.5,0.5),ylim=c(-1,1))
plot(log10(subset$betaE),slopes,ylim=c(-1,1),main=gene)
abline(h=0)

