# use python3 
# based on PyStan2

# https://mc-stan.org/docs/2_18/stan-users-guide/linear-regression.html

import sys,math
import pystan

def logsigmoid(x): return math.log(x/(1-x),10)


#############################################################

# mod = lm(logsigmoid(abund)~log10(betaEpos)+log2(newconc),data=melted)

model = """
data {
  int<lower=0> N;   // number of data items
  vector[N] logsigmoid_abund;
  vector[N] log_betaEpos;
  vector[N] log_newconc;
}
parameters {
  real alpha;           // intercept
  real beta1;
  real beta2;
  real<lower=0> sigma;  // error scale
}
model {
  logsigmoid_abund ~ normal(alpha + beta1 * log_betaEpos + beta2 * log_newconc, sigma);  // likelihood
}
"""

#############################################################

#schools_data = {"J": 8,
#                "y": [28,  8, -3,  7, -1,  1, 18, 12],
#                "sigma": [15, 10, 16, 11,  9, 11, 10, 18]}

data = []
skip = 1
for line in open(sys.argv[1]):
  if skip>0: skip -= 1; continue
  w = line.rstrip().split('\t')
  data.append(w)

betaEpos = [float(w[2]) for w in data]
abund =    [float(w[4]) for w in data]
newconc =  [float(w[5]) for w in data]

data = { "N":3696, "K":2 }
data["logsigmoid_abund"] = [logsigmoid(a) for a in abund]
data["log_betaEpos"] = [math.log(b,10) for b in betaEpos]
data["log_newconc"] = [math.log(c,2) for c in newconc]


#############################################################

# for PyStan2
sm = pystan.StanModel(model_code=model)
fit = sm.sampling(data=data, iter=1000, chains=4, seed=1)
print(fit)


