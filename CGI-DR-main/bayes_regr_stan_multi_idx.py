# use python3 
# based on PyStan2

# https://mc-stan.org/docs/2_18/stan-users-guide/linear-regression.html
# https://datascienceplus.com/bayesian-regression-with-stan-part-1-normal-regression/

#https://num.pyro.ai/en/0.7.1/tutorials/bayesian_hierarchical_linear_regression.html  (lung capacity in patients)
#note: I am assuming there is no covariance between features, hence sigmas are independent; otherwise use cov_matrix and wishart...

#https://m-clark.github.io/bayesian-basics/enhancements.html  (robust regr, uses cauchy for sig)
#alpha ~ normal(0,10);
#for (k in 1:K) { beta[k] ~ normal(mu_b[k],var_b[k]); }
#sigma ~ cauchy(0,5);

import sys,math,numpy
import pystan

def logsigmoid(x): return math.log(x/(1-x),10)

def HDI(vals,conf=0.95):
  n = len(vals)
  a = int(n*conf)
  b = n-a # e.g. 5%
  temp = sorted(vals)
  besti,minwidth = -1,-1
  for i in range(b):
    width = temp[i+a]-temp[i]
    if i==0 or width<minwidth: besti,minwidth = i,width
  return (temp[besti],temp[besti+a])

#############################################################

# here's the model in R notation:
#   mod = lm(logsigmoid(abund)~log10(betaEpos)+log2(newconc),data=melted)

model = """
data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of covariates (columns in X)
  int<lower=0> G;   // number of genes
  vector[N] Y;
  matrix[N,K] X; 
  int<lower=1,upper=G> Z[N]; // gene indexes (1..G)
}
parameters {
  real alpha[G];           // intercept
  vector[K] beta[G];       // slopes
  real<lower=0> var_e[G];  // error variances (one for each gene)

  real mu_b[K]; // distribution of slopes for each covariate
  real<lower=0> var_b[K]; // variances of slope distns
}
model 
{
  alpha ~ normal(0,10); // for each gene
  var_e ~ inv_gamma(1,1); // error variance for each gene
  mu_b ~ normal(0,10); // prior distribution of slopes for each covariate
  var_b ~ inv_gamma(1,1);
  for (k in 1:K) { 
    beta[k] ~ normal(mu_b,var_b); // slopes (vector over covars) for each gene, drawn from hierarchical distribution
  }
  for (n in 1:N) {
    Y[n] ~ normal(alpha[Z[n]] + X[n]*beta[Z[n]], sqrt(var_e[Z[n]])); // see examples in Sec 5.8 of stan-reference-2.6.2-1.pdf
  }
}
"""

#############################################################

genes = {} # index by ints 1..G
data = []
skip = 0
for line in open(sys.argv[1]):
  if skip>0: skip -= 1; continue
  if line.startswith("#"): 
    w = line.rstrip().split()
    i = int(w[1])
    genes[i] = w[2:] # orf, gene, num_sgRNAs
    continue
  w = line.rstrip().split('\t')
  data.append(w)

logsigmoid_abund = [float(w[0]) for w in data]
log_betaEpos =     [float(w[1]) for w in data]
log_newconc =      [float(w[2]) for w in data]
#Z = [[int(x) for x in w[3:]] for w in data]
Z = [int(w[3]) for w in data]

N = len(data)
G = max([int(x[-1]) for x in data]) # last column is gene indexes
K = 2
print("N=%s, K=%s, G=%s" % (N,K,G))

data = { "N":N, "K":K, "G":G }
data["Y"] = logsigmoid_abund
data["X"] = [ [b,c] for (b,c) in zip(log_betaEpos,log_newconc)]
data["Z"] = Z

#############################################################

# for PyStan2
sm = pystan.StanModel(model_code=model)
fit = sm.sampling(data=data, iter=1000, chains=4, seed=1)
print(fit)

samples = fit.extract()
beta = samples['beta'] # 2000xGx2
N = beta.shape[0]
for i in range(G):
  slopes = beta[:,i,1]
  neg,pos = numpy.sum(slopes<0)/float(N),numpy.sum(slopes>0)/float(N)
  orf,gene,num_sgRNAs = genes[i+1]
  lb,ub = HDI(slopes)
  flag = "."
  if numpy.mean(slopes)<0 and pos<0.05: flag = '*'
  if numpy.mean(slopes)>0 and neg<0.05: flag = '*'
  overlap = lb<0 and ub>0 # overlap 0 or ROPE?
  flag2 = "."
  if overlap==False: flag2 = "$"
  print("%s %-10s %-10s %4s %8.3f %8.3f %8.3f %8.3f %8.3f %5.3f %5.3f %s %s" % (i+1,orf,gene,num_sgRNAs,numpy.min(slopes),lb,numpy.mean(slopes),ub,numpy.max(slopes),neg,pos,flag,flag2))
