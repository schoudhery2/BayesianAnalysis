import sys
from Spreadsheet import *

betaE = {}
skip = 1
for line in open("Bosch21_TableS2.txt"):
  if skip>0: skip -= 1; continue
  w = line.rstrip().split('\t')
  if len(w[9])<2: continue # some betaE entries are empty
  id,b = w[0],float(w[9])
  id = id[:id.rfind("_")] # strip off v4PAMscore...
  betaE[id] = b

#####################

drug,days = sys.argv[1],sys.argv[2]

metadata = Spreadsheet("ShiquiCGI_metadata.txt")

files = []
concs = []
for i in range(metadata.nrows):
  if metadata.get(i,"drug")==drug and metadata.get(i,"days_predep")==days:
    files.append(metadata.get(i,"filename"))
    concs.append(float(metadata.get(i,"conc_xMIC")))

if days=="1": files.append("counts_1962_DMSO_D1.txt"); concs.append(0) # what about different pools?
if days=="5": files.append("counts_1972_DMSO_D5.txt"); concs.append(0)
if days=="10": files.append("counts_1952_DMSO_D10.txt"); concs.append(0)

pairs = sorted(zip(concs,files))
concs,files = [x[0] for x in pairs],[x[1] for x in pairs]

for i in range(len(files)):
  sys.stderr.write("%s\t%s\n" % (concs[i],files[i]))

#####################

no_dep = {}
IDs = []
for line in open("no_depletion_abundances.txt"):
  w = line.rstrip().split('\t')
  id,abund = w[0],float(w[1])
  id = id[:id.rfind("_")]
  no_dep[id] = abund
  IDs.append(id)

#####################

# also need to get DMSO to represent 0xMIC

Abund = []
Concs = [] # one for each replicate, include those for DMSO (0)
for file,conc in zip(files,concs):
  counts = {} # rows of 3, indexed by id
  skip = 1
  for line in open("data/counts/%s" % file):
    if skip>0: skip -= 1; continue
    if "Negative" in line or "Empty" in line: continue
    w = line.rstrip().split('\t')
    id = w[0]
    id = id[:id.rfind("_")]
    vals = [int(x) for x in w[1:]]
    counts[id] = vals # rows of 3
  for i in range(3): # assume there are 3 replicates in each file
    cnts = [x[i] for x in counts.values()] # a column
    tot = sum(cnts)
    abund = []
    for id in IDs: abund.append(counts[id][i]/float(tot)) # a column for fracs, parallel to IDs
    Abund.append(abund)
    Concs.append(conc)

print('\t'.join("orf gene id abund betaE".split()+[str(x) for x in Concs]))
a,b = 0,0
for i,id in enumerate(IDs):
  #vals = [id]+["%0.8f" % x[i] for x in Abund]
  PC = 0.00000001
  if id not in betaE: a += 1; continue; #sys.stderr.write("%s not in betaE\n" % id); continue
  if id not in no_dep: b += 1; continue # sys.stderr.write("%s not in no_dep\n" % id); continue
  temp = id[:id.find("_")].split(":") # "RVBD00067:rpoB"
  orf,gene = temp[0],temp[1]
  vals = [orf,gene,id,"%0.8f" % (no_dep[id]),str(betaE[id])]+["%0.6f" % ((x[i]+PC)/(no_dep[id]+PC)) for x in Abund]
  print('\t'.join(vals))

sys.stderr.write("warning: betaE values not found for %s gRNAs\n" % a)
sys.stderr.write("warning: no_dep values not found for %s gRNAs\n" % b)

