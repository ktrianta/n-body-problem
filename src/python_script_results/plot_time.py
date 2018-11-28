from __future__ import division
from numpy import *
from pylab import *
from decimal import *
import math
import matplotlib.pyplot as plt
from matplotlib import rc

filename1 = sys.argv[1]

# ----- Font
rc('font',size=11)
rc('font',family='sans serif')
rc('axes',labelsize=14)
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams['mathtext.rm'] = 'cm'

# ----- Make data.
NP,time1,time2,time3 = genfromtxt(filename1, unpack=True,delimiter='',skip_header=2)

plt.plot(NP,time1,'--o',NP,time2,'--x',NP,time3,'--^')

plt.xlabel(r'$\#(\mathrm{Processors}) \ p$',labelpad = 2 )
plt.ylabel(r'$\mathrm{Time}(s)$')
plt.xticks((1,2,3,4,5,6,7),(1,2,4,8,16,32,64))
#plt.tight_layout()
plt.title(r'$N = 200 \ \mathrm{particules}$')
plt.legend([r'$\mathrm{Force \  Computation}$',r'$\mathrm{Tree \  Construction}$',r'$\mathrm{Tree \  Sending}$'],frameon=False);
save_results_to = '../python_script_results/plots/'
plt.savefig(save_results_to + "time.eps")
