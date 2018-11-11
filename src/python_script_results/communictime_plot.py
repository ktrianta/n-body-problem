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
numProcessors,execTime = genfromtxt(filename1, unpack=True,delimiter='')

plt.plot(numProcessors,execTime,'g^',color='red')
plt.xlabel(r'$\#(\mathrm{Processors}) \ p$')
plt.ylabel(r'$\mathrm{Communication \ Time}(s)$')
plt.tight_layout()
save_results_to = '../python_script_results/plots/'
plt.savefig(save_results_to + 'communictime.pdf')
