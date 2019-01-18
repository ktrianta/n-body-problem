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
x1,y1 = genfromtxt("naive-sequential-one-sided.plotData",delimiter=' ',unpack=True)
x2,y2 = genfromtxt("naive-parallel-one-sided.plotData",delimiter=' ',unpack=True)
x3,y3 = genfromtxt("barnes-sequential-one-sided.plotData",delimiter=' ',unpack=True)
x4,y4 = genfromtxt("barnes-parallel-one-sided.plotData",delimiter=' ',unpack=True)
x5,y5 = genfromtxt("barnes-parallel-balanced.plotData",delimiter=' ',unpack=True)
x6,y6 = genfromtxt("barnes-sequential-balanced.plotData",delimiter=' ',unpack=True)





fig,ax = plt.subplots()
plt.text(2**10,10**12,r'$\pi_{Par} = 8*\pi_{\mathrm{XeonE5_{2680v3}}}$')
plt.text(2**10,10**11,r'$\pi_{\mathrm{XeonE5_{2680v3}}}$')
plt.text(2**10,9.5**15,r'$\beta_{\mathrm{XeonE5_{2680v3}}} I$',rotation=38)


ax.set_xscale('log',basex=2)
ax.set_yscale('log',basey=2)
ax.plot(x5,y5,'kv')
ax.plot(x6,y6,'yo')
ax.plot(x1,y1,'ro')
ax.plot(x2,y2,'cs')
ax.plot(x3,y3,'b^')
ax.plot(x4,y4,'g>')

ax.legend([r'$\mathrm{Barnes \ Hut \ balanced, \  parallel} $', r'$\mathrm{Barnes \ Hut \ balanced, \  sequential} $', r'$\mathrm{Naive \  algorithm, \ sequential} $', r'$\mathrm{Naive \  algorithm, \ parallel} $', r'$\mathrm{Barnes \ Hut, \  sequential} $', r'$\mathrm{Barnes\ Hut, \  parallel} $'],frameon=False)

x = np.arange(0,100000000,1000)
l1 = [8*16*2500000000]*len(x)
l2 = 8*68000000000*x
l3 = [16*2500000000]*len(x)
l4 = 68000000000*x

ax.plot(x,l1,'b')
ax.plot(x,l2,'b')
ax.plot(x,l3,'r')
ax.plot(x,l4,'r')

ax.set_xlabel(r'$I \ \mathrm{[FLOPs/byte]}$')
ax.set_ylabel(r'$P \ \mathrm{[FLOPs/cycle]}$')

plt.ylim(10000000, 1000000000000000000)
plt.xlim(0.001,2**23)

#plt.show()
plt.savefig("roofline.eps",transparent=True)
