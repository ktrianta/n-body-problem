file = open('sequential/FP_OPS/naive-parallel-one-sided.data') 

i = 0
f =  file.readlines()
m = len(f)/10
fp_ops = [0]*m

for n in f:
    fp_ops[i/10] += int(n)
    i += 1
fp_ops = [ x / 10 for x in fp_ops ]

#for n in range(0,m):
#    print min([ int(x) for x in f[n*10:(n+1)*10] ])
#    fp_ops[n] = min([ int(x) for x in f[n*10:(n+1)*10] ])

file = open('sequential/TOT_CYC/naive-parallel-one-sided.data') 

i = 0
f =  file.readlines()
m = len(f)/10
tot_cyc = [0]*m

#for n in f:
#    tot_cyc[i/10] += int(n)
#    i += 1
#tot_cyc = [ x / 10 for x in tot_cyc ]

for n in range(0,m):
    tot_cyc[n] = min([ int(x) for x in f[n*10:(n+1)*10] ])


file = open('sequential/L3_TCM/naive-parallel-one-sided.data') 

i = 0
f =  file.readlines()
m = len(f)/10
l3_tcm = [0]*m

#for n in f:
#    l3_tcm[i/10] += int(n)
#    i += 1
#l3_tcm = [ x / 10 for x in l3_tcm ]

for n in range(0,m):
    l3_tcm[n] = min([ int(x) for x in f[n*10:(n+1)*10] ])

X = [ fp_ops[x]/float(l3_tcm[x]) for x in range(0, m) ]
Y = [ fp_ops[x]/float(tot_cyc[x])*2500000000 for x in range(0,m) ]

file = open('naive-sequential-one-sided.plotData', 'w')
for i in range(0,m):
    file.write(str(X[i]) + " " + str(Y[i]) + "\n")
file.close()

#barnes-one-sided
file = open('sequential/FP_OPS/barnes-parallel-one-sided.data') 

i = 0
f =  file.readlines()
m = len(f)/10
fp_ops = [0]*m

for n in f:
    fp_ops[i/10] += int(n)
    i += 1
fp_ops = [ x / 10 for x in fp_ops ]

#for n in range(0,m):
#    print min([ int(x) for x in f[n*10:(n+1)*10] ])
#    fp_ops[n] = min([ int(x) for x in f[n*10:(n+1)*10] ])

file = open('sequential/TOT_CYC/barnes-parallel-one-sided.data') 

i = 0
f =  file.readlines()
m = len(f)/10
tot_cyc = [0]*m

#for n in f:
#    tot_cyc[i/10] += int(n)
#    i += 1
#tot_cyc = [ x / 10 for x in tot_cyc ]

for n in range(0,m):
    tot_cyc[n] = min([ int(x) for x in f[n*10:(n+1)*10] ])


file = open('sequential/L3_TCM/barnes-parallel-one-sided.data') 

i = 0
f =  file.readlines()
m = len(f)/10
l3_tcm = [0]*m

#for n in f:
#    l3_tcm[i/10] += int(n)
#    i += 1
#l3_tcm = [ x / 10 for x in l3_tcm ]

for n in range(0,m):
    l3_tcm[n] = min([ int(x) for x in f[n*10:(n+1)*10] ])

X = [ fp_ops[x]/float(l3_tcm[x]) for x in range(0, m) ]
Y = [ fp_ops[x]/float(tot_cyc[x])*2500000000 for x in range(0,m) ]

file = open('barnes-sequential-one-sided.plotData', 'w')
for i in range(0,m):
    file.write(str(X[i]) + " " + str(Y[i]) + "\n")
file.close()

#barnes-balanced
file = open('sequential/FP_OPS/barnes-parallel-balanced.data') 

i = 0
f =  file.readlines()
m = len(f)/10
fp_ops = [0]*m

for n in f:
    fp_ops[i/10] += int(n)
    i += 1
fp_ops = [ x / 10 for x in fp_ops ]

#for n in range(0,m):
#    print min([ int(x) for x in f[n*10:(n+1)*10] ])
#    fp_ops[n] = min([ int(x) for x in f[n*10:(n+1)*10] ])

file = open('sequential/TOT_CYC/barnes-parallel-balanced.data') 

i = 0
f =  file.readlines()
m = len(f)/10
tot_cyc = [0]*m

#for n in f:
#    tot_cyc[i/10] += int(n)
#    i += 1
#tot_cyc = [ x / 10 for x in tot_cyc ]

for n in range(0,m):
    tot_cyc[n] = min([ int(x) for x in f[n*10:(n+1)*10] ])


file = open('sequential/L3_TCM/barnes-parallel-balanced.data') 

i = 0
f =  file.readlines()
m = len(f)/10
l3_tcm = [0]*m

#for n in f:
#    l3_tcm[i/10] += int(n)
#    i += 1
#l3_tcm = [ x / 10 for x in l3_tcm ]

for n in range(0,m):
    l3_tcm[n] = min([ int(x) for x in f[n*10:(n+1)*10] ])

print [ fp_ops[x]/float(tot_cyc[x]) for x in range(0,m) ]
X = [ fp_ops[x]/float(l3_tcm[x]) for x in range(0, m) ]
Y = [ fp_ops[x]/float(tot_cyc[x])*2500000000 for x in range(0,m) ]

file = open('barnes-sequential-balanced.plotData', 'w')
for i in range(0,m):
    file.write(str(X[i]) + " " + str(Y[i]) + "\n")
file.close()
