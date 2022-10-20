args = commandArgs(trailingOnly=TRUE)
options(width=200)

ifile = args[1]
data = read.table(ifile,head=T,sep='\t')

ORFs = sort(unique(data$orf))
THEGENE = "???"
if (length(args)>1) 
{
  #ORFs = c(args[2])
  if (length(args)>1) THEGENE = c(args[2])
  #temp = data[data$gene==args[2],]
  #ORFs = c(temp$orf[1])
}

logsigmoid = function(x) { log10(x/(1-x)) }

######################################
# estimate noise level
#  (apparently, this has to be done before squashing)

NOISE = F
if (NOISE)
{
 print("fitting variances...")

 SCV = NULL
 info = NULL
 rows = sample(1:nrow(data),5000) # sampling
 temp = data[rows,]
 for (i in 1:nrow(temp))
 {
  row = temp[i,]
  for (i in 0:3) # same as log2 of "fake" concs below
  {
    colstart = 6+i*3
    id = as.character(row[3])
    vals = as.numeric(row[colstart:(colstart+2)])
    m = mean(vals)
    s = sd(vals)
    scv = s/m
    m2 = mean(logsigmoid(vals))
    v2 = var(logsigmoid(vals))
    SCV = rbind(SCV,c(m=m,s=s,logm=log10(m),scv=scv,lsig_mean=m2,lsig_var=v2))
    info = rbind(info,c(id=id,logconc=i))
  }
 }
 SCV = as.data.frame(SCV)
 nzdata = SCV[SCV$m>0 & SCV$s>0,]
 fit = loess(lsig_var~lsig_mean,data=nzdata)
 ypred = predict(fit,nzdata)

 #plot(nzdata$lsig_mean,nzdata$lsig_var,ylim=c(0,2.5))
 #par(new=T)
 #plot(nzdata$lsig_mean,ypred,col=2,ylim=c(0,2.5))

 # interpolate the LOESS fit so I can look up predicted var based on mean of logsigmoid
 lsig_mean = as.numeric(seq(-3,2,0.1))
 temp = as.data.frame(lsig_mean)
 lsig_predvar = predict(fit,temp)
 noise = data.frame(lsig_mean,predvar=lsig_predvar) # not needed
 print("noise:")
 print(noise)
 #print(approx(lsig_mean,lsig_predvar,xout=-2)) # 0.1543382

 noisefn = function(X) { approx(lsig_mean,lsig_predvar,xout=X)$y } # interpolating a lookup table
}

###################################
# squashing function:
#   see Wolfram Alpha: plot (1-exp(-2x))/(1+exp(-2x)) from -1 to 3

PC = 0.01
for (i in 6:17)
{
  data[,i] = PC+(1-PC)*(1-exp(-2*data[,i]))/(1+exp(-2*data[,i]))
}

#summary(data$betaE)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.71179 -0.30520 -0.06266 -0.20121 -0.01427  0.84667 

#data$betaEpos = -(data$betaE-max(data$betaE))
data$betaEpos = data$betaE-min(data$betaE)+0.01

######################################

require(reshape2)
require(lmtest)

for (orf in ORFs)
{
 subset = data[data$orf==orf,]
 nobs = dim(subset)[1]
 gene = as.character(subset[1,"gene"])
 if (THEGENE!="???" & gene!=THEGENE) { next }

 cat("----------------------\n")
 cat(sprintf("%s %s %s\n",orf,gene,nobs))
 if (gene=="ligC") { print("skipping"); next } # fit nearly perfect because some frac abune >30

 # note: although concs are like 0, 0.0625, 0.125, 0.25 (different range for each drug)
 #  treat them as 1, 2, 4, 8 for simplicity (it doesn't matter, as long as they are 2 fold dilutions)
 #  this also assumes "0" is half of the lowest concentration 

 temp = cbind(subset[,c("id","betaEpos")],subset[,6:17])
# print("temp:")
# print(dim(temp))
 colnames(temp)=c("id","betaEpos","a1","a2","a3","b1","b2","b3","c1","c2","c3","d1","d2","d3")
 melted = melt(temp,id.vars=c("id","betaEpos"),variable.name="conc",value.name="abund")
# print("melted:")
# print(dim(melted))

 melted$newconc = 0
 melted$newconc[melted$conc=="a1"] = 1
 melted$newconc[melted$conc=="a2"] = 1
 melted$newconc[melted$conc=="a3"] = 1
 melted$newconc[melted$conc=="b1"] = 2
 melted$newconc[melted$conc=="b2"] = 2
 melted$newconc[melted$conc=="b3"] = 2
 melted$newconc[melted$conc=="c1"] = 4
 melted$newconc[melted$conc=="c2"] = 4
 melted$newconc[melted$conc=="c3"] = 4
 melted$newconc[melted$conc=="d1"] = 8
 melted$newconc[melted$conc=="d2"] = 8
 melted$newconc[melted$conc=="d3"] = 8

 means = aggregate(logsigmoid(abund)~id+newconc,data=melted,FUN=mean)
 a = colnames(means); a[3] = "mean"; colnames(means) = a
 vars = aggregate(logsigmoid(abund)~id+newconc,data=melted,FUN=var)
 a = colnames(vars); a[3] = "var"; colnames(vars) = a

# print("subset")
# print(head(subset))
# print("means")
# print(dim(means))
# print(head(means))
# print("variances")
# print(head(vars))

 X = round(means[,3],1)
 if (NOISE) 
 { 
   predvars = noisefn(X)
   weights = 1/predvars

   # print("X:")
   # print(X)
   # print("predvars:")
   # print(head(predvars))
   # print("weights:")
   # print(head(weights))

   temp = means
   temp$var = vars$var
   temp$predvar = predvars
   temp$weight = weights
   # print("temp")
   # print(temp)
 }
 else { weights = 1 }

 #w = temp[temp$id==melted$id & temp$newconc==melted$newconc,"weight"] # why won't this work?!
 f = function(row) { temp[temp$id==row[1] & temp$newconc==row[2],"weight"] }
 w = apply(melted[,c("id","newconc")],1,f)
 melted$weight = w

 write.table(melted,"melted.txt",sep='\t',quote=F); print("writing melted.txt")

 tryCatch( 
  {
   #mod = lm(logsigmoid(abund)~log10(betaEpos)+log2(newconc),data=melted,weights=w) # weighted regression
   #print(length(w))
   #print(length(melted$newconc))
   mod = lm(logsigmoid(abund)~log10(betaEpos)+log2(newconc),data=melted)
   print(summary(mod))
   summ = summary(mod)

   # LRT gives essentially same P-value as Wald test
   #m0 =  lm(logsigmoid(abund)~log10(betaEpos),data=melted) # reduced model, without conc term
   #print(summary(m0))
   #print(lrtest(mod,m0))
 
   coeffs = summ$coefficients[,1]
   pvals = summ$coefficients[,4]
   cat(sprintf("result: %s %s %s %s %s %s %s %s %s\n",orf,gene,nobs,coeffs[1],coeffs[2],coeffs[3],pvals[1],pvals[2],pvals[3]))
  },
  error = function (e) { print("skipping due to error with lm")  } 
 )
}
