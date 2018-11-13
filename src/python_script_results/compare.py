from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import rc
import sys
import argparse
import csv
import itertools

# Sample execution command: python script.py ../../test/resources/tab1024.txt ../../test/resources/tab1024.txt 7 7 3

filename_me = sys.argv[1]
filename2_other = sys.argv[2]
num_of_columns = int(sys.argv[3])
epsilon  = float(sys.argv[4])

if num_of_columns == 3:
    rx_me,ry_me,rz_me = genfromtxt(filename_me,unpack=True,delimiter=' ')
    rx_other,ry_other,rz_other = genfromtxt(filename_other,unpack=True,delimiter=' ')
elif num_of_columns == 6:
    rx_me,ry_me,rz_me,vx_me,vy_me,vz_me = genfromtxt(filename_me,unpack=True,delimiter=' ')
    rx_other,ry_other,rz_other,vx_other,vy_other,vz_other = genfromtxt(filename_other,unpack=True,delimiter=' ')

#For accuracy, use np.square(), for speed use **2
r_me = sqrt(rx_me**2+ry_me**2+rz_me**2)
r_other = sqrt(rx_other**2+ry_other**2+rz_other**2)
v_me = sqrt(vx_me**2+vy_me**2+vz_me**2)
v_other = sqrt(vx_other**2+vy_other**2+vz_other**2)

difference_error_r = 2 * abs(r_me - r_other) / ( r_me + r_other )
L2_error_r = LA.norm(difference_error_r, 2)     # sqrt(sum(e_i),1<=i<=num_particules)
LInf_error_r = LA.norm(difference_error_r, Inf) # max [e_i], 1<=i<=num_particules

if ( L2_error_r > epsilon ) :
    print("No convergence")
else:
    print("The solution is accurate")
