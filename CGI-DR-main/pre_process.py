import sys,math

PC = 0.01

def squashing(x): return PC+(1-PC)*(1-math.exp(-2*x))/(1+math.exp(-2*x))

def logsigmoid(x): return math.log(x/(1-x),10)

newconcs = [1,1,1,2,2,2,4,4,4,8,8,8]

data = []
skip = 1
for line in open(sys.argv[1]):
  if skip>0: skip -= 1; continue
  w = line.rstrip().split('\t')
  data.append(w)

betaE = [float(w[4]) for w in data]
m = min(betaE)
betaEpos = [x-m+0.01 for x in betaE]

if len(sys.argv)<3: sys.stderr.write("error: give one or more genes on command line\n"); sys.exit(0)

for i in range(len(data)):
  bEp = betaEpos[i]
  w = data[i]
  gene = w[1]
  Z = [1 if w[1]==g else 0 for g in sys.argv[2:]]
  if sum(Z)==1:
    for j in range(12):
      newconc = newconcs[j]
      abund = float(w[5+j])
      logsigmoid_abund = logsigmoid(squashing(abund))
      vals = [logsigmoid_abund, math.log(bEp,10), math.log(newconc,2)]+Z
      print '\t'.join([str(x) for x in vals])
  
