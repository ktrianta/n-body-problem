from numpy import *
from numpy import linalg as LA
import sys

# Sample execution command: python compare.py filename_me filename_other num_of_columns epsilon
# Note that the file inputs should be in the same directory as the compare.py file. Otherwise, provide the path 
# of each file in the terminal. 

filename_me = sys.argv[1]
filename_other = sys.argv[2]
num_of_columns = int(sys.argv[3])
epsilon  = float(sys.argv[4])

rx_me,ry_me,rz_me,vx_me,vy_me,vz_me = genfromtxt(filename_me,unpack=True,delimiter='    ')
rx_other,ry_other,rz_other,vx_other,vy_other,vz_other = genfromtxt(filename_other,unpack=True,delimiter='    ')

#For accuracy, use np.square(), for speed use **2
r_me = sqrt(rx_me**2+ry_me**2+rz_me**2)
r_other = sqrt(rx_other**2+ry_other**2+rz_other**2)
v_me = sqrt(vx_me**2+vy_me**2+vz_me**2)
v_other = sqrt(vx_other**2+vy_other**2+vz_other**2)

difference_error_r =  100 * 2 * abs(r_me - r_other) / ( r_me + r_other )
L2_error_r = LA.norm(difference_error_r, 2)     # sqrt(sum(e_i),1<=i<=num_particules)
LInf_error_r = LA.norm(difference_error_r, Inf) # max [e_i], 1<=i<=num_particules

print("The error difference is: ")
print(LInf_error_r)

if ( L2_error_r > epsilon ) :
    print("No convergence")
else:
    print("The solution is accurate")
