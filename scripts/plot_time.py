from __future__ import division
from matplotlib import ticker
from numpy import *
from pylab import *
from decimal import *
import math
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker as mtick
#  python plot_time.py naive_parallelOneSided.txt naive_parallelAvx256.txt
filename1 = sys.argv[1]
filename2 = sys.argv[2]

# ----- Font
rc('font',size=11)
rc('font',family='sans serif')
rc('axes',labelsize=14)
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams['mathtext.rm'] = 'cm'


# ----- Make data.
# Parallel one Sided Naive
total_oneSided, compute_oneSided, IO_oneSided, comm_oneSided = genfromtxt("../eulerRuns/naive/"+ filename1,
                                                                unpack=True,delimiter=',',skip_header=0,skip_footer=0)
total_oneSided = total_oneSided - IO_oneSided
sort_total_oneSided = (np.sort(total_oneSided))[::-1]

# Parallel Avx256 Naive
total_par256, compute_par256, IO_par256, comm_par256 = genfromtxt("../eulerRuns/naive/"+filename2, 
                                                                unpack=True,delimiter=',',skip_header=0,skip_footer=0)
total_par256 = total_par256 - IO_par256
sort_total_par256 = (np.sort(total_par256))[::-1]

avg_total_oneSided = np.zeros(7)
avg_compute_oneSided = np.zeros(7)
avg_comm_oneSided = np.zeros(7)

avg_total_par256 = np.zeros(7)
avg_compute_par256 = np.zeros(7)
avg_comm_par256 = np.zeros(7)

median_total_oneSided = np.zeros(7)
median_compute_oneSided = np.zeros(7)
median_comm_oneSided = np.zeros(7)

median_total_par256 = np.zeros(7)
median_compute_par256 = np.zeros(7)
median_comm_par256 = np.zeros(7)

errorCI_total_oneSided = np.zeros( (2,7) )
errorCI_total_par256 = np.zeros( (2,7) )

alpha = 0.05
RUNS = 50
SQRUNS = math.sqrt(RUNS)
z_alpha0verTwo = 1.96 
index_low = int( math.floor( (RUNS - z_alpha0verTwo*SQRUNS)/2.0 ) )
index_up = int( math.ceil( 1.0 + (RUNS + z_alpha0verTwo*SQRUNS)/2.0 ) )


INDEX = [0,50,100,150,200,250,300]
j = 0
for i in INDEX:
    avg_total_oneSided[j] = np.mean(total_oneSided[i:RUNS+i])
    avg_compute_oneSided[j] = np.mean(compute_oneSided[i:RUNS+i])
    avg_comm_oneSided[j] = np.mean(comm_oneSided[i:RUNS+i])
    avg_total_par256[j] = np.mean(total_par256[i:RUNS+i])
    avg_comm_par256[j] = np.mean(comm_par256[i:RUNS+i])
    avg_compute_par256[j] = np.mean(compute_par256[i:RUNS+i])

    median_comm_oneSided[j] = np.median(comm_oneSided[i:RUNS+i])
    median_total_oneSided[j] = np.median(total_oneSided[i:RUNS+i])
    median_compute_oneSided[j] = np.median(compute_oneSided[i:RUNS+i]) 
    median_comm_par256[j] = np.median(comm_par256[i:RUNS+i])
    median_total_par256[j] = np.median(total_par256[i:RUNS+i])
    median_compute_par256[j] = np.median(compute_par256[i:RUNS+i])
    
    errorCI_total_oneSided[0][j] = -sort_total_oneSided[i+index_low] + median_total_oneSided[j] #error below median
    errorCI_total_oneSided[1][j] = sort_total_oneSided[i+index_up] -  median_total_oneSided[j] #error above median

    errorCI_total_par256[0][j] = -sort_total_par256[i+index_low] + median_total_par256[j]
    errorCI_total_par256[1][j] = sort_total_par256[i+index_up] - median_total_par256[j]
    
    j = j+1

PR = [1,2,4,8,16,32,48]

speedup_oneSided = np.zeros(6)

speedup_par256 = np.zeros(6)

for i in range(0,6):
    speedup_oneSided[i] = median_total_oneSided[0]  /  median_total_oneSided[i+1]
    speedup_par256[i] = median_total_par256[0]  / median_total_par256[i+1] 


fig,ax = plt.subplots()
#ax.plot(PR,avg_total_oneSided,'--o',PR,avg_total_par256,'--o',PR,avg_comm_oneSided,'--o',PR,avg_comm_par256,'--o')
#ax.plot(PR[1:],speedup_oneSided,'-o',PR[1:],speedup_par256,'-o')
#ax.plot(PR[1:],PR[1:],'-s')
#ax.plot(PR[1:],median_comm_oneSided[1:],'--o',PR[1:],median_comm_par256[1:],'--o')
ax.plot(PR[1:],median_total_oneSided[1:],'-o',PR[1:],median_total_par256[1:],'-o')
#ax.errorbar(PR,median_total_par256,errorCI_total_par256,capsize = 4,marker='s',linestyle="None")
#ax.errorbar(PR,median_total_oneSided,errorCI_total_oneSided,capsize = 4,marker='s',linestyle="None" )
ax.plot(PR[1:],median_total_oneSided[0]/PR[1:],'-s')
ax.plot(PR[1:],median_total_par256[0]/PR[1:],'-s')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'#$\mathrm{(Processors)}p$')
ax.set_ylabel(r'$\mathrm{Time}(s)$')
#ax.set_ylabel(r'$\mathrm{Speed \ up}$')
#ax.legend([r'$\mathrm{Parallel}$', r'$\mathrm{Parallel \ AVX256}$',r'$\mathrm{Theoretical \ speed \  up}$'],frameon=False,loc=2);
ax.legend([r'$\mathrm{Parallel}$',r'$\mathrm{Parallel \ AVX256}$',r'$\mathrm{Parallel \ (Theoretical)}$',r'$\mathrm{Parallel \ AVX256 \ (Theoretical)}$'],frameon=False);
ax.set_xticks([2,4,8,16,32,48])
#ax.set_xticks([1,2,4,8,16,32,48])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

for axis in [ax.xaxis, ax.yaxis]:
    formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
    axis.set_major_formatter(formatter)
ax.grid(True,linestyle='--', linewidth=0.2)

plt.tight_layout()
#plt.title(r'$N = 200 \ \mathrm{particules}$')
#formatter = ticker.ScalarFormatter(useMathText=True)
#formatter.set_scientific(True) 
#formatter.set_powerlimits((10,10)) 
#ax.yaxis.set_major_formatter(formatter) 
#plt.show()
#ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
save_results_to = '../scripts/plots/'
plt.savefig(save_results_to + "total.eps",transparent=True)
