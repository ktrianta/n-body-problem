from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import rc
import sys
import argparse

filename1 = sys.argv[1]
filename2 = sys.argv[2]
num_of_columns = int(sys.argv[3])
theta  = float(sys.argv[4])

if num_of_columns == 3: 
    rx_me,ry_me,rz_me = genfromtxt(filename1,unpack=True,delimiter=' ',skip_header=0,skip_footer=0)
    rx_other,ry_other,rz_other = genfromtxt(filename2,unpack=True,delimiter=' ',skip_header=0,skip_footer=0)
else if num_of_colums == 6:
    rx_me,ry_me,rz_me,vx_me,vy_me,vz_me = genfromtxt(filename1,unpack=True,delimiter=' ',skip_header=0,skip_footer=0)
    rx_other,ry_other,rz_other,vx_other,vy_other,vz_other = genfromtxt(filename2,unpack=True,delimiter=' ',skip_header=0,skip_footer=0)

error = abs(X-Y)

L2_error = LA.norm(error,2)     # sqrt(sum(e_i),1<=i<=num_particules)
LInf_error = LA.norm(error,Inf) # max [e_i], 1<=i<=num_particules

if ( L2_error > theta ) :
    print("No convergence")
else:
    print("The solution is accurate")
