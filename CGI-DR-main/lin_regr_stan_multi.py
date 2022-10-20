# use python3 
# based on PyStan2

# https://mc-stan.org/docs/2_18/stan-users-guide/linear-regression.html
# https://datascienceplus.com/bayesian-regression-with-stan-part-1-normal-regression/

import sys,math
import pystan

def logsigmoid(x): return math.log(x/(1-x),10)


#############################################################

# mod = lm(logsigmoid(abund)~log10(betaEpos)+log2(newconc),data=melted)

model = """
data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of covariates (columns in X)
  int<lower=0> G;   // number of genes
  vector[N] Y;
  matrix[N,K] X; // covariates (betaE, conc)
  matrix[N,G] Z; // bit vectors for genes
}
parameters {
  vector[G] alpha;     // intercept
  matrix[G,K] beta;    // slopes
  vector[G] sigma;     // std deviations (one for each gene)
}
model {
  // Z is NxG; alpha is Gx1; Zxalpha is Nx1; X is NxK; beta is GxK; Z*beta is NxK
  Y ~ normal(Z*alpha + rows_dot_product(X,Z*beta), Z*sigma); 
}
"""

#############################################################

data = []
skip = 0
for line in open(sys.argv[1]):
  if skip>0: skip -= 1; continue
  w = line.rstrip().split('\t')
  data.append(w)

logsigmoid_abund = [float(w[0]) for w in data]
log_betaEpos =     [float(w[1]) for w in data]
log_newconc =      [float(w[2]) for w in data]
Z = [[int(x) for x in w[3:]] for w in data]

N = len(data)
G = len(data[0])-3

data = { "N":N, "K":2, "G":G }
data["Y"] = logsigmoid_abund
data["X"] = [ [b,c] for (b,c) in zip(log_betaEpos,log_newconc)]
data["Z"] = Z

#############################################################

# for PyStan2
sm = pystan.StanModel(model_code=model)
fit = sm.sampling(data=data, iter=1000, chains=4, seed=1)
print(fit)


