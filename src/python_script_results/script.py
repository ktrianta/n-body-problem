from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import rc
import sys

filename1 = sys.argv[1]
filename2 = sys.argv[2]
num_of_columns = int(sys.argv[3])
epsilon = float(sys.argv[4])

X,Y = genfromtxt(filename1,unpack=True,delimiter=' ',skip_header=0,skip_footer=0)
#X,Y = genfromtxt(filename2,unpack=True,delimiter=' ',skip_header=0,skip_footer=0)

print(X[0],Y[0])
print(X[1],Y[1])

error = abs(X-Y)

print(error)

print(LA.norm(error,2))
print(LA.norm(error,inf))


if ( LA.norm(error,2) > epsilon) :
    print("No convergence")
