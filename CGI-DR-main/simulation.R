#########################
# simulation

logsigmoid = function(x) { log10(x/(1-x)) }

concs = c(0.0625,0.125,0.25,0.5,1)
betaEs = seq(-2,0,0.1)

Kb = -1 # midpoint of betaE range
Hb = -4

Kc = 0.25 # midpoint of conc range
Hc = -0.5

vals = NULL
for (b in betaEs)
{
  for (c in concs)
  {
    theta = 1/(1+((Kb/b)**Hb)*(Kc/c)**Hc)
    vals = rbind(vals,c(conc=c,betaE=b,theta=theta))
  }
}
print(vals)
vals = as.data.frame(vals)
par(mfrow=c(1,2))
for (b in betaEs)
{
  plot(theta~log(conc),data=vals[vals$betaE==b,],type='l',ylim=c(0,1))
  par(new=T)
}
par(new=F)
for (b in betaEs)
{
  plot(logsigmoid(theta)~log(conc),data=vals[vals$betaE==b,],type='l',ylim=c(-2,2))
  par(new=T)
}
