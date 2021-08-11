#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 13:05:27 2021

@author: bmondal
"""

import numpy as np
from scipy.spatial.transform import Rotation
from numpy.linalg import norm
from scipy import linalg
import argparse
import sys, os

#%%
parser = argparse.ArgumentParser(prog='MakeVacuum', description='This script creates new vasp cordinate file with vacuum. Input: POSCAR/CONTCAR/.vasp file or supply 3 lattice parameters)',
                                        epilog='Have fun!')
parser.add_argument('-filename', type=argparse.FileType('r'), default='POSCAR', help='The file path e.g. /home/mondal/VASP/test.vasp')
parser.add_argument('-lp', nargs='+', type=float, help='Lattice parameters in x, y, and z directions respectively. ')
parser.add_argument('-SC', action='store_true', default=False, help='Is the cell a supercell? True/False (default: False)')
parser.add_argument('-SF', nargs='+', type=int, default=[6,6,6], help='Supercell dimensions (must be int) in a, b and c lattice vector directions respectively. (default: [6,6,6])')
parser.add_argument('-CS', type=str, default='zincblende', help='Crystal Structure e.g. cubic, zincblende, strechedzincblende etc. (default: zincblende)')
parser.add_argument('-plane', type=int, default=110, help='Lattice plane over which the vaccum will be. (default: 110)')
parser.add_argument('-layer', type=int, default=1, help='How many layer you want to repeat along vacuum direction? (default: 1)')
parser.add_argument('-vac', type=float, default=0, help='Vacuum layer length (default: 0). Default unit in same unit of lattice vector unit. Use -Lvac for unit in layer number.')
parser.add_argument('-Lvac', action='store_true',default=False, help='To switch the vaccum length unit in no. of layer. It can be fraction.')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)

args = parser.parse_args()

filename = args.filename
lp_ary = np.asarray(args.lp, dtype=float)
SuperCell = args.SC
S_factor = np.asarray(args.SF, dtype=float)
cst = args.CS
plane110 = args.plane
repeat_layer = args.layer
vacuum_layer = args.vac
Layer_vacuum = args.Lvac


nfname = filename.name+'_new.vasp'
try:
    os.remove(nfname)
except OSError:
    pass
fnamm = open(nfname,'a')

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
def lattice_parameter(vector,crys_str='zincblende'):
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

if lp_ary.size > 1:
    lattice = np.copy(lp_ary)
    system = ['New Structure']
    sname = ['XX   XX']
else:
    print('\nFILENAME\t=\t%s\nSupercell\t=\t%s'%(filename.name,SuperCell),'\nSupercell Dimensions\t=\t',S_factor,'\n')

    lines = filename.readlines()
    system=[lines[0].rstrip()]
    factor = float(lines[1].split()[0])
    a = np.asarray(lines[2].split()).astype(float)
    b = np.asarray(lines[3].split()).astype(float)
    c = np.asarray(lines[4].split()).astype(float)
    sname = [lines[5].rstrip()]
    N_atom = np.sum(np.asarray(lines[6].split()).astype(int))
    
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
    print('Lx\t=\t%7.5f\nLy\t=\t%7.5f\nLz\t=\t%7.5f'%(lattice[0],lattice[1],lattice[2]))
#%%

def unitcell_coordinate_zincblende(plane):
    # unitcell contains the minimum no. of atom considering periodicity in the two perpendicular
    # directions of the vacuum plane. 
    # Format: The atom coordinates should be grouped together
    if plane==110:
        # 1st two coordinates are fractional coordinates for Ga, next teo are for As
        unitcell = np.array([[0.000000000,         0.000000000,         0.000000000],
                             [0.250000000,         0.250000000,         0.500000000],
                             [0.250000000,         0.000000000,         0.250000000],
                             [0.000000000,         0.250000000,         0.750000000]])
        natom = 2 # No. of Ga atom (= No. of As atom)
        new_lvec = np.array([[1,1,0],[-1,1,0],[0,0,1]])
        lattice_rescale2 = np.array([0.5, 0.5, 1])
    if plane==100:
        unitcell = np.array([[0.000000000,         0.000000000,         0.000000000],
                             [0.000000000,         0.500000000,         0.500000000],
                             [0.500000000,         0.000000000,         0.500000000],
                             [0.500000000,         0.500000000,         0.000000000],
                             [0.250000000,         0.250000000,         0.250000000],
                             [0.750000000,         0.750000000,         0.250000000],
                             [0.750000000,         0.250000000,         0.750000000],
                             [0.250000000,         0.750000000,         0.750000000]])
        natom = 4 # No. of Ga atom (= No. of As atom)
        new_lvec = np.array([[1,0,0],[0,1,0],[0,0,1]])
        lattice_rescale2 = np.array([1, 1, 1])
    '''
    # This part is wrong.
    if plane==111:
        unitcell = np.array([[0.000000000,         0.000000000,         0.000000000],
                             [0.333333343,         0.250000000,         0.083333336],
                             [0.666666687,         0.000000000,         0.166666642],
                             [0.250000000,         0.000000000,         0.000000000],
                             [0.583333313,         0.250000000,         0.083333336],
                             [0.916666687,         0.000000000,         0.166666642]])
        natom = 3 # No. of Ga atom (= No. of As atom)
        new_lvec = np.array([[1,1,1],[-1,1,0],[-1,-1,2]])
        lattice_rescale2 = np.array([1, 1/3, 0.25])
    '''    
    lattice_rescale1 = norm(new_lvec,axis=1)
    
    return unitcell, natom, lattice_rescale1,lattice_rescale2

def rotated_new_pos_zincblende(lattice, plane):
    unitcell, natom, lattice_rescale1,lattice_rescale2 = unitcell_coordinate_zincblende(plane)    
    lattice *= lattice_rescale1
    unitcell *= lattice
    posGa = unitcell[:natom]
    posAs = unitcell[natom:]
    lattice *= lattice_rescale2
    return lattice, posGa, posAs
    

if cst == 'zincblende':
    lattice, posGa, posAs = rotated_new_pos_zincblende(lattice, plane110)
    add_array = np.array([lattice[0], 0, 0])
    # #%%
    # IDENTITY = np.identity(3)
    # axis = [0, 0, 1]
    # theta = np.deg2rad(45)
    # axis = np.array(axis, dtype=float)
    # axis = axis / norm(axis)  # normalize the rotation vector first
    # rot = Rotation.from_rotvec(theta * axis)
    
    # new_v = rot.apply(IDENTITY) 
    # axis = [0, 1, 0]
    # theta = np.deg2rad(45)
    # axis = np.array(axis, dtype=float)
    # axis = axis / norm(axis)  # normalize the rotation vector first
    # rot = Rotation.from_rotvec(theta * axis)
    # new_v = rot.apply(new_v)
    # print(new_v)
    
    #%%

    # Create the layers    
    unitcellGa_core = np.copy(posGa)
    unitcellAs_core = np.copy(posAs)
    for I in range(repeat_layer-1):
        unitcellGa_core += add_array
        posGa = np.concatenate((posGa,unitcellGa_core), axis=0)
        
        unitcellAs_core += add_array
        posAs = np.concatenate((posAs,unitcellAs_core), axis=0)
    
    # New no. of Ga and As
    nGa = len(posGa)
    nAs = len(posAs) 

    if Layer_vacuum:
        vacuum_layer *= lattice[0]

    # New lattice vectors
    lattice[0] *= repeat_layer
    lattice[0] += vacuum_layer
    lattice_new_vec = np.zeros((3, 3))
    np.fill_diagonal(lattice_new_vec, lattice)
         
print("\n----------------- Creating final file -----------------------")   
print("* Please change the line 6 and 7 in the '_new.vasp' file if needed. For diamond structure add the 2 species in line 6 and 7. \n")

np.savetxt(fnamm,system,fmt='%s')
np.savetxt(fnamm,[1])
np.savetxt(fnamm,lattice_new_vec)
np.savetxt(fnamm,sname,fmt='%s')
np.savetxt(fnamm,[[nGa, nAs]],fmt='%d')
np.savetxt(fnamm,['Cartesian'],fmt='%s')
np.savetxt(fnamm,posGa)
np.savetxt(fnamm,posAs)

fnamm.close()
