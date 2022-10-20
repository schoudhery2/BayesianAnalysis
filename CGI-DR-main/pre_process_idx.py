import sys,math

PC = 0.01

def squashing(x): return PC+(1-PC)*(1-math.exp(-2*x))/(1+math.exp(-2*x))

def logsigmoid(x): return math.log(x/(1-x),10)

newconcs = [1,1,1,2,2,2,4,4,4,8,8,8]

data = []
orfs = {}
genes = {} # orf->gene
skip = 1
for line in open(sys.argv[1]):
  if skip>0: skip -= 1; continue
  if len(sys.argv)>2:
    keywords = sys.argv[2:]
    found = False
    for k in keywords:
      if k in line: found = True
    if found==False: continue
  w = line.rstrip().split('\t')
  orf,gene = w[0],w[1]
  data.append(w)
  if orf not in orfs: orfs[orf] = 0
  orfs[orf] += 1
  genes[orf] =gene

orfids = sorted(orfs.keys())

for i,orf in enumerate(orfids):
  print("# %s %s %s %s" % (i+1,orf,genes[orf],orfs[orf]))

betaE = [float(w[4]) for w in data]
m = min(betaE)
betaEpos = [x-m+0.01 for x in betaE]

for i in range(len(data)):
  bEp = betaEpos[i]
  w = data[i]
  orf,gene = w[0],w[1]
  #Z = [1 if gene==g else 0 for g in sys.argv[2:]]
  #if sum(Z)==1:
  
  Z = 1+orfids.index(orf)
  for j in range(12):
      newconc = newconcs[j]
      abund = float(w[5+j])
      logsigmoid_abund = logsigmoid(squashing(abund))
      vals = [logsigmoid_abund, math.log(bEp,10), math.log(newconc,2),Z]
      print '\t'.join([str(x) for x in vals])
  
