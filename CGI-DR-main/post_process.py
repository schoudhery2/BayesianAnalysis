import sys
from statsmodels.stats.multitest import fdrcorrection

data = []
for line in open(sys.argv[1]):
  if line.startswith("result:"):
    w = line.split()
    pval = float(w[-1]) if w[-1]!="NA" else 1
    data.append((pval,w))

pvals = [x[0] for x in data]
qvals = fdrcorrection(pvals)[1] # I assume this is Benjamini-Hochberg method
data = [(x[0],x[1]+[q]) for x,q in zip(data,qvals)]

data.sort()
print '\t'.join("rank orf gene num_sgRNAs coeff_intercept coeff_log_betaE coeff_log_conc Pval_intercept Pval_log_betaE Pval_log_conc Qval_log_conc".split())
rank = 1
for (pval,w) in data:
  vals = [rank]+w[1:]
  print '\t'.join([str(x) for x in vals])
  rank += 1
