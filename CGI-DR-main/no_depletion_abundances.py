import sys

data = []
skip = 1
for line in open("counts_616_H37Rv_pass_INPUT-plATc.txt"):
  if skip>0: skip -= 1; continue
  if "Negative" in line or "Empty" in line: continue
  w = line.rstrip().split('\t')
  data.append(w)

A = [float(x[1]) for x in data]
B = [float(x[2]) for x in data]
C = [float(x[3]) for x in data]

totA = sum(A)
totB = sum(B)
totC = sum(C)

fracA = [x/float(totA) for x in A]
fracB = [x/float(totB) for x in B]
fracC = [x/float(totC) for x in C]

for i in range(len(data)):
  id = data[i][0]
  mn = (fracA[i]+fracB[i]+fracC[i])/3.0
  print "%s\t%0.8f" % (id,mn)
