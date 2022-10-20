args = commandArgs(trailingOnly=TRUE)

ifile = args[1]
data = read.table(ifile,head=T,sep='\t')

logsigmoid = function(x) { log10(x/(1-x)) }

#####################
# noise level

SCV = NULL
rows = sample(1:nrow(data),5000)
temp = data[rows,]
for (i in 1:nrow(temp))
{
  row = temp[i,]
  for (i in 0:3)
  {
    colstart = 6+i*3
    vals = as.numeric(row[colstart:(colstart+2)])
    m = mean(vals)
    s = sd(vals)
    scv = s/m
    m2 = mean(logsigmoid(vals))
    v2 = var(logsigmoid(vals))
    SCV = rbind(SCV,c(m=m,s=s,logm=log10(m),scv=scv,lsig_mean=m2,lsig_var=v2))
  }
}
SCV = as.data.frame(SCV)
nzdata = SCV[SCV$m>0 & SCV$s>0,]
plot(scv~logm,data=nzdata,main="noise among replicates",ylim=c(0,1.5),xlab="log10(mean fractional abundance)",ylab="SCV=stdev/mean (among replicates)")
fit = loess(scv~logm,data=nzdata)
ypred = predict(fit,nzdata)
par(new=T)
plot(nzdata$logm,ypred,col=2,xlab="",ylab="",ylim=c(0,1.5),pch=20)

plot(lsig_var~lsig_mean,data=nzdata,main="noise among replicates",xlab="mean logsigmoid (of fractional abundance among replicates)",ylab="var(logsigmoid)",ylim=c(0,2.5))
fit = loess(lsig_var~lsig_mean,data=nzdata)
ypred = predict(fit,nzdata)
par(new=T)
plot(nzdata$lsig_mean,ypred,col=3,xlab="",ylab="",ylim=c(0,2.5),pch=20)
