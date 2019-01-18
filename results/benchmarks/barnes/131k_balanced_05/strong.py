from __future__ import division
from numpy import *
from pylab import *
from decimal import *
import math
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker
from scipy.stats import mstats


# python ../eulerRuns/naive/naive_parallelOneSided.txt ./eulerRuns/naive/naive_parallelAvx256.txt plot_time.py
# ----- Font
rc('font',size=11)
rc('font',family='sans serif')
rc('axes',labelsize=14)
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams['mathtext.rm'] = 'cm'


# ----- Make data.

NUM_MEAS = 50
SQNUM_MEAS = np.sqrt(NUM_MEAS)
z_alpha0verTwo = 3
index_low = int(math.floor( (NUM_MEAS - z_alpha0verTwo*SQNUM_MEAS)/2.0 ))
index_up = int( math.ceil( 1.0 + (NUM_MEAS + z_alpha0verTwo*SQNUM_MEAS)/2.0 ) )
NUM_PR = np.zeros(7)
NUM_PR = [1,2,4,8,16,32,48]

med = np.zeros(NUM_MEAS)
global_median_b05 = np.zeros( len(NUM_PR) )
global_median_b025 = np.zeros( len(NUM_PR) )
global_median_p05 = np.zeros( len(NUM_PR) )
global_median_p025 = np.zeros( len(NUM_PR) )
global_median_n = np.zeros( len(NUM_PR) )
global_median_avx = np.zeros( len(NUM_PR) )

errorCI_total = np.zeros((2,7))
errorCI_comm = np.zeros((2,7))
errorCI_comp = np.zeros((2,7))
errorCI_gather = np.zeros((2,7))
errorCI_tree= np.zeros((2,7))

j = 0
for p in NUM_PR:
    sample = np.zeros( (p,NUM_MEAS) )
    for i in range(0,p):
        sample[:][i] = genfromtxt(str(p)+"/plot"+str(i)+".txt",usecols = 0,unpack=True,delimiter=',',skip_header=1,skip_footer=0)

    medians = np.median(sample,axis = 0)
    global_median_b025[j] = np.median(medians)
    s = np.sort(medians)
    j = j+1

j = 0
for p in NUM_PR:
    sample = np.zeros( (p,NUM_MEAS) )
    for i in range(0,p):
        sample[:][i] = genfromtxt("../131k_balanced_05/"+str(p)+"/plot"+str(i)+".txt",usecols = 0,unpack=True,delimiter=',',skip_header=1,skip_footer=0)
    medians = np.median(sample,axis = 0)
    global_median_b05[j] = np.median(medians)
    s = np.sort(medians)
    j = j+1


j = 0
for p in NUM_PR:
    sample = np.zeros( (p,NUM_MEAS) )
    for i in range(0,p):
        sample[:][i] = genfromtxt("../131k_parallel_05/"+str(p)+"/plot"+str(i)+".txt",usecols = 0,unpack=True,delimiter=',',skip_header=1,skip_footer=0)
    medians = np.median(sample,axis = 0)
    global_median_p05[j] = np.median(medians)
    s = np.sort(medians)
    j = j+1

j = 0
for p in NUM_PR:
    sample = np.zeros( (p,NUM_MEAS) )
    for i in range(0,p):
        sample[:][i] = genfromtxt("../131k_parallel_025/"+str(p)+"/plot"+str(i)+".txt",usecols = 0,unpack=True,delimiter=',',skip_header=1,skip_footer=0)
    medians = np.median(sample,axis = 0)
    global_median_p025[j] = np.median(medians)
    s = np.sort(medians)
    j = j+1

global_n = genfromtxt("../../naive/131k_parallel.txt",usecols = 0, unpack=True,delimiter=',',skip_header=0,skip_footer=0)

INDEX = [0,50,100,150,200,250,300]
j = 0
for i in INDEX:

    global_median_n[j] = np.median(global_n[i:NUM_MEAS+i])
    
    j = j+1
global_avx = genfromtxt("../../naive/131k_parallel_avx_256.txt",usecols = 0, unpack=True,delimiter=',',skip_header=0,skip_footer=0)

INDEX = [0,50,100,150,200,250,300]
j = 0
for i in INDEX:

    global_median_avx[j] = np.median(global_avx[i:NUM_MEAS+i])
    
    j = j+1
fig,ax = plt.subplots()
#ax.plot(PR,avg_total_oneSided,'--o',PR,avg_total_par256,'--o',PR,avg_comm_oneSided,'--o',PR,avg_comm_par256,'--o')
#ax.plot(PR,median_total_oneSided,'--o',PR,median_total_par256,'--o',PR,median_comm_oneSided,'--o',PR,median_comm_par256,'--o')
#ax.plot(NUM_PR[1:],global_median_total[1:],'--o',NUM_PR[1:],global_median_comm[1:],'--o',NUM_PR[1:],global_median_comp[1:],'--o',NUM_PR[1:],global_median_gather[1:],'--o',NUM_PR[1:],global_median_tree[1:],'--o')
#f = 5/100
#Y = [global_median_total[0]*( (1.0 - f) / n + f) for n in NUM_PR]
#ax.plot(NUM_PR[1:],Y[1:],'--s')

ax.plot(NUM_PR[1:],global_median_n[1:],'--o')
ax.plot(NUM_PR[1:],global_median_avx[1:],'--s')
ax.plot(NUM_PR[1:],global_median_p05[1:],'--x')
#ax.plot(NUM_PR[1:],global_median_p025[1:],'--o')
#ax.plot(NUM_PR[1:],global_median_b025[1:],'--o')
ax.plot(NUM_PR[1:],global_median_b05[1:],'--^')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\#(\mathrm{Processors}) \ p$')
ax.set_ylabel(r'$\mathrm{Time}(s)$')
#plt.tight_layout()
#plt.title(r'$N = 200 \ \mathrm{particules}$')
ax.legend([r'$\mathrm{Naive \ Implementation}$',r'$\mathrm{Naive-AVX256}$',r'$\mathrm{Equal \ \# \ of \ Particlesi, \ \theta=0.5 }$',r'$\mathrm{Equal \ Amount \ of \ Work, \ \theta=0.5}$'],frameon=False);
ax.set_xticks([2,4,8,16,32,48])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

textstr = r' $N=130k \ \mathrm{particles}$'
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.05,0.05, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='bottom', bbox=props)


for axis in [ax.xaxis, ax.yaxis]:
    formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
    axis.set_major_formatter(formatter)
ax.grid(True,linestyle='--', linewidth=0.2)
plt.tight_layout()
#plt.ylim(0.001, 700)
#plt.show()
plt.savefig("130kComparison.eps",transparent=True)
