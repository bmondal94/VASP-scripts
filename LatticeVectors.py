#!/usr/bin/python
import numpy as np
import argparse
import sys
# -------------------------------- Parsers ------------------------------------------
parser = argparse.ArgumentParser(prog='LatticeVectors', description='This script calculates the different cell parameters from VASP file (POSCAR/CONTCAR format)',
					epilog='Have fun!')
parser.add_argument('filename', metavar='FILENAME', type=argparse.FileType('r'), help='The file path e.g. /home/mondal/VASP/test.vasp')
parser.add_argument('-UC', action='store_true', default=False, help='Is the cell already a conventional unit cell? True/False (default: False)')
parser.add_argument('-SC', action='store_true', default=False, help='Is the cell a supercell? True/False (default: False)')
parser.add_argument('-SF', nargs='+', type=int, default=[6,6,6], help='Supercell dimensions (must be int) in a, b and c lattice vector directions respectively. (default: [6,6,6])')
parser.add_argument('-CS', type=str, default='cubic', help='Crystal Structure e.g. cubic, zincblende, strechedzincblende etc. (default: cubic)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

# ---------------------------------Functions -------------------------------------------
def vec(a):
	return (np.sqrt(np.sum(a**2)))

def anglefn(v1, v2, acute=True):
# v1 is your firsr vector
# v2 is your second vector
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    if (acute == True):
        return angle
    else:
        return 2 * np.pi - angle

"""
For streched zincblende:
x, y, z == Lattice parameter of streched zincblende structure.
xy, yz, and xz == Diagonals of streched zincblende or 2 times lattice vactors
Note1: Angle b/w x, y, z are always 90 degree in (streched)zincblende structure
	Angle b/w xy, xz and yz might be differ from 60 in streched case
Note2: 6x6 primitive supercell is equivalent to 3x3 unit cell for zincblende

(2*xy)**2 = x**2 + y**2; (2*yz)**2 = y**2 + z**2; (2*xz)**2 = x**2 + z**2

Solving above 3 equations we get;
	y**2 = xy**2 + yz**2 - xz**2
	x**2 = xy**2 + xz**2 - yz**2
	z**2 = yz**2 + xz**2 - xy**2
"""
def lattice_parameter(vector,crys_str='cubic'):
	if(crys_str == 'cubic'):
		lp = vector
	elif(crys_str == 'zincblende'):
		lp = vector * np.sqrt(2)
	elif(crys_str == 'strechedzincblende'): 
		lp = np.array([0.,0.,0.])
		nv = vector**2; factor = np.sqrt(2)
		lp[1] =  factor * np.sqrt(nv[0] + nv[1] - nv[2])
		lp[0] =  factor * np.sqrt(nv[0] + nv[2] - nv[1])
		lp[2] =  factor * np.sqrt(nv[1] + nv[2] - nv[0])
	print('Crystal Structure (conventional unit cell)\t=\t%s\n'%crys_str)
	return (lp)

# ----------------------------------------------------------------------------------------

if __name__=="__main__":
	# ************************ Read parsers ***************************************
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	#try:
	#	args = parser.parse_args()
	#except:
	#	#parser.print_help()
	#	sys.exit(0)

	args = parser.parse_args()

	filename = args.filename
	UNIT_CELL = args.UC
	SuperCell = args.SC
	S_factor = np.asarray(args.SF, dtype=np.float)
	cst = args.CS

	# ********************** Read coordinate file ************************************
	print('\nFILENAME\t=\t%s\nSupercell\t=\t%s'%(filename.name,SuperCell),'\nSupercell Dimensions\t=\t',S_factor,'\n')

	lines = filename.readlines()
	factor = float(lines[1].split()[0])
	a = np.asarray(lines[2].split()).astype(np.float)
	b = np.asarray(lines[3].split()).astype(np.float)
	c = np.asarray(lines[4].split()).astype(np.float)
	N_atom = np.sum(np.asarray(lines[6].split()).astype(np.int))
	
	# ********************** Primitive --> Unit cell *********************************
	vector = np.array([vec(a),vec(b),vec(c)])
	vector *= factor

	angle = np.array([anglefn(a,b), anglefn(b,c), anglefn(a,c)])

	C_Volm = np.prod(vector)*np.sqrt(1. - np.sum(np.cos(angle)**2) + 2.*np.prod(np.cos(angle)))

	angle = np.degrees(angle)

	print('Cell parameters (Lattice vectors and angles):\n')
	print('Total number of atoms (/ions) in the cell\t=\t%d\n'%N_atom)
	print('a\t=\t%7.5f\nb\t=\t%7.5f\nc\t=\t%7.5f\n'%(vector[0],vector[1],vector[2]))
	print('alpha\t=\t%7.5f\nbeta\t=\t%7.5f\ngamma\t=\t%7.5f\n'%(angle[0],angle[1],angle[2]))
	print('Cell Volume\t=\t%8.4f\n'%C_Volm)

	if (SuperCell):
		primitive_vector = np.divide(vector,S_factor)  #S_factor=Supercell factor
	else:
		primitive_vector = vector
	if  UNIT_CELL:
		print("** Note: The lattice is already a conventional unit cell")
		lattice = primitive_vector
	else:
		print('Primitive cell lattice vectors:\n')
		print('a\t=\t%7.5f\nb\t=\t%7.5f\nc\t=\t%7.5f\n'%(primitive_vector[0],primitive_vector[1],primitive_vector[2]))

		print("Note: The convention used here for lattice parameters calculation:\n\t 'a' should be in xy-plane \n\t 'b' should be in yz-plane \n\t 'c' hsould be in xz-plane.\n")
		print("** WARNING: If you don't follow above convention then depending on how you define lattice vectors the dummy index x, y, z in L may change\n")
		print("** You might need to rotate the coordinate system to follow the convention. Please check the lattice vectors in input file (e.g. POSCAR).\n")

		rotate = input("If you need the rotation please provide here the number of rotation you need to get to the above condition. Else press 'ENTER'.\n") 

		if(rotate):
			primitive_vector = np.roll(primitive_vector,int(rotate))

		lattice = lattice_parameter(primitive_vector, crys_str=cst)
		#lattice = np.divide(np.array([a[0],b[1],c[2]]),S_factor)*(factor*2.0) # Use the projection of lattice vectors to get lattice constants
	print('(Conventional) Unit cell lattice constants: \n')
	print('Lx\t=\t%7.5f\nLy\t=\t%7.5f\nLz\t=\t%7.5f\n'%(lattice[0],lattice[1],lattice[2]))
