from numpy import *
from numpy import square
from numpy import linalg as LA
import sys

# Sample execution command: python compare.py filename_me filename_other epsilon
# Note that the file inputs should be in the same directory as the compare.py file.
# Otherwise, provide the path of each file in the terminal. 

filename_me = sys.argv[1]
filename_other = sys.argv[2]
epsilon = float(sys.argv[3])

rx_me, ry_me, rz_me = genfromtxt(filename_me,unpack=True,delimiter='    ',usecols=(0,1,2))
rx_other, ry_other, rz_other = genfromtxt(filename_other,unpack=True,delimiter='    ',usecols=(0,1,2))

# For accuracy, use np.square(), for speed use **2
r_me = sqrt(square(rx_me) + square(ry_me) + square(rz_me))
r_other = sqrt(square(rx_other) + square(ry_other) + square(rz_other))

abs_difference_error_r = zeros(len(r_me))
rel_difference_error_r = zeros(len(r_me))

for i in range(0,len(r_me)):
    diffr = r_me[i] + r_other[i]
    if ( diffr == 0 ):
        abs_difference_error_r[i] = 0
        rel_difference_error_r[i] = 0
    else: 
        abs_difference_error_r[i] = abs(r_me[i] - r_other[i]) 
        rel_difference_error_r[i] = 2 * abs(r_me[i] - r_other[i]) / diffr

L2_error_a = LA.norm(abs_difference_error_r, 2)
LInf_error_a = LA.norm(abs_difference_error_r, Inf)

print "Maximum absolute error difference:", LInf_error_a
print "Total absolute error difference:", L2_error_a

L2_error_r = LA.norm(rel_difference_error_r, 2)
LInf_error_r = LA.norm(rel_difference_error_r, Inf)

print "Maximum relative error difference:", LInf_error_r
print "Total relative error difference:", L2_error_r

if ( L2_error_r > epsilon ) :
    print("No convergence")
else:
    print("The solution is accurate")
