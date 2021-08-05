#!/usr/bin/python
import numpy as np
from collections import Counter
import argparse
import sys

'''
v' == New row vector (atom position vector)
v == Old row vector (atom position vector)
B == Lattice vector matrix
A == Bassis vector matix

   v'.(B.A) = v.A
=> v' = (v.A).(B.A)^(-1)
'''
parser = argparse.ArgumentParser(prog='sqs2poscarMyPy', description='This script converts ATAT sqs.out file to VASP POSCAR',
                                        epilog='Have fun!')
parser.add_argument('filename', metavar='FILENAME', type=argparse.FileType('r'), help='The sqs file path e.g. /home/mondal/bestsqs.out')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)
#try:
#       args = parser.parse_args()
#except:
#       #parser.print_help()
#       sys.exit(0)

args = parser.parse_args()

fname = args.filename.name

f = open("POSCAR","w+")

#Load data from a text file using np.genfromtxt
#Each row in the text file must have the same number of values.

# Basis vectors
B = np.genfromtxt(fname, dtype=float, unpack=False, max_rows=3)
# Lattice vectors
L = np.genfromtxt(fname, dtype=float, unpack=False, skip_header=3, max_rows=3) 
# Calculate L.B
LB = np.dot(L,B)
BL_inv = np.linalg.inv(LB)
# Fractional positions w.r.t above basis vectors, read from file
d = np.genfromtxt(fname, dtype=float,usecols=(0,1,2), unpack=False, skip_header=6)
# Modified Fractional positions w.r.t new basis vectors; v'=(v.B).(B.L)^(-1)
vf = np.array([np.dot(np.dot(v,B),BL_inv) for v in d]) 
# Read the elements from file
e = np.genfromtxt(fname, dtype=str, usecols=3, unpack=False, skip_header=6)
# Sort elements
sort_elements = np.argsort(e)
elements_types = Counter(e[sort_elements])
elements = [list(elements_types.keys())]
elements_number = [list(elements_types.values())]
# Sort final position vectors according to elements sort
vfs = vf[sort_elements]

# Prepare POSCAR
f.write("POSCAR from "+fname+'\n'+'1.00000000\n')
np.savetxt(f,LB,fmt='%.16f')
np.savetxt(f,elements,fmt='%s')
np.savetxt(f,elements_number,fmt='%d')
f.write("Direct\n")
np.savetxt(f,vfs,fmt='%.16f')
f.close()
