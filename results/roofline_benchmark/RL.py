file = open('sequential/FP_OPS/naive-parallel-one-sided.data') 

i = 0
f =  file.readlines()
m = len(f)/10
fp_ops = [0]*m

for n in f:
    fp_ops[i/10] += int(n)
    i += 1

fp_ops = [ x / 10 for x in fp_ops ]

file = open('sequential/TOT_CYC/naive-parallel-one-sided.data') 

f =  file.readlines()
m = len(f)/10
tot_cyc = [0]*m

for n in range(0,m):
    print min([ int(x) for x in f[n*10:(n+1)*10] ])
    tot_cyc[n] = min([ int(x) for x in f[n*10:(n+1)*10] ])


file = open('sequential/L3_TCM/naive-parallel-one-sided.data') 

i = 0
f =  file.readlines()
m = len(f)/10
l3_tcm = [0]*m

for n in f:
    l3_tcm[i/10] += int(n)
    i += 1

l3_tcm = [ x / 10 for x in l3_tcm ]


X = [ fp_ops[x]/float(l3_tcm[x]) for x in range(0, m) ]
Y = [ fp_ops[x]/float(tot_cyc[x])*2500000000 for x in range(0,m) ]

file = open('naive-sequential-one-sided.plotData', 'w')
for i in range(0,m):
    file.write(str(X[i]) + " " + str(Y[i]) + "\n")
file.close()
