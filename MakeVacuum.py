#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 13:05:27 2021

@author: bmondal
"""

import numpy as np
from numpy.linalg import norm
import argparse
import sys, os

#%% --------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(prog='MakeVacuum', description='This script creates new vasp cordinate file with vacuum. Input: POSCAR/CONTCAR/.vasp file or supply 3 lattice parameters',
                                        epilog='Have fun!')
parser.add_argument('-filename', type=str, default='POSCAR', help='The file path e.g. /home/mondal/VASP/test.vasp')
parser.add_argument('-lp', nargs='+', type=float, default=[],help='Lattice parameters in x, y, and z directions respectively. If a single number is given then that no. will be applied to all 3 lattice parameters.')
parser.add_argument('-UC', action='store_true', default=False, help='Is the cell already a conventional unit cell? True/False (default: False)')
parser.add_argument('-SC', action='store_true', default=False, help='Is the cell a supercell? True/False (default: False)')
parser.add_argument('-SF', nargs='+', type=int, default=[6,6,6], help='Supercell dimensions (must be int) in a, b and c lattice vector directions respectively. (default: [6,6,6])')
parser.add_argument('-CS', type=str, default='zincblende', help='Crystal Structure e.g. cubic, zincblende, strechedzincblende etc. (default: zincblende)')
parser.add_argument('-plane', type=str, default=110, help='Lattice plane over which the vaccum will be. (default: 110)')
parser.add_argument('-layer', nargs='+', type=int, default=[1,1,1], help='How many layer you want to repeat? If single no. is given then repetation will be applied along z-axis only. If 2 no.s are given repetation will be in xy-plane. (default: no repetation)')
parser.add_argument('-vac', type=float, default=0, help='Vacuum layer (along z-axis) length (default: 0). Default unit in same unit of lattice vector unit. Use -Lvac for unit in layer number.')
parser.add_argument('-Lvac', action='store_true',default=False, help='To switch the vaccum length unit in no. of layer. It can be fraction.')
parser.add_argument('--version', action='version', version='%(prog)s 2.0')

#%% ----------------------------------------------------------------------------------------------------
## All function  definitions
def vec(a):
        return (np.sqrt(np.sum(a**2)))

def anglefn(v1, v2, acute=True):
# v1 is your firsr vector
# v2 is your second vector
    angle = np.arccos(np.dot(v1, v2) / (norm(v1) * norm(v2)))
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

def read_file(filename, SuperCell, S_factor, UNIT_CELL):
    print("* Reading lattice vectors from file.")
    print("*****************************************")
    print('FILENAME\t=\t%s\nSupercell\t=\t%s'%(filename,SuperCell),'\nSupercell Dimensions\t=\t',S_factor,'\n')
    ffname = open(filename,'r')
    lines = ffname.readlines()
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

    print('(Conventional) Unit cell lattice constants: \n')
    print('Lx\t=\t%7.5f\nLy\t=\t%7.5f\nLz\t=\t%7.5f'%(lattice[0],lattice[1],lattice[2]))
    ffname.close()
    print("*******************************************")
    return system, sname, lattice

def unitcell_coordinate_zincblende(plane):
    # unitcell contains the minimum no. of atom considering periodicity in the two perpendicular
    # directions of the vacuum plane. 
    # Format: The atom coordinates should be grouped together
    print("* Surface plane is %s"%(plane))
    if plane=='011' or plane=='101' or plane=='110':
        # 1st two coordinates are fractional coordinates for Ga, next teo are for As
        unitcell = np.array([[0.000000000,         0.000000000,         0.000000000],
                             [0.500000000,         0.250000000,         0.250000000],
                             [0.250000000,         0.000000000,         0.250000000],
                             [0.750000000,         0.250000000,         0.000000000]])
        natom = 2 # No. of Ga atom (= No. of As atom)
        new_lvec = np.array([[1,0,0],[0,1,-1],[0,1,1]])
        lattice_rescale2 = np.array([1, 0.5, 0.5])
        #if plane=='101':
        #    unitcell[:,[0,1]] = unitcell[:,[1,0]]
        #    new_lvec = np.array([[1,0,-1],[0,1,0],[1,0,1]])
        #    lattice_rescale2[[0,1]] = lattice_rescale2[[1,0]]
        #else:
        #    unitcell[:,[0,2]] = unitcell[:,[2,0]]
        #    new_lvec = np.array([[1,1,0],[-1,1,0],[0,0,1]])
        #    lattice_rescale2[[0,2]] = lattice_rescale2[[2,0]]
    elif plane=='100' or plane=='010' or plane=='001':
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
    elif plane=='111':
        unitcell = np.array([[0.000000000,         0.000000000,         0.000000000],
                             [0.166666667,         0.000000000,         0.333333333],
                             [0.083333333,         0.250000000,         0.666666666],
                             [0.250000000,         0.250000000,         0.000000000],
                             [0.416666667,         0.250000000,         0.333333333],
                             [0.333333333,         0.000000000,         0.666666667],
                             [0.000000000,         0.000000000,         0.250000000],
                             [0.166666667,         0.000000000,         0.583333333],
                             [0.083333333,         0.250000000,         0.916666667],
                             [0.250000000,         0.250000000,         0.250000000],
                             [0.416666667,         0.250000000,         0.583333333],
                             [0.333333333,         0.000000000,         0.916666667]])
        natom = 6 # No. of Ga atom (= No. of As atom)
        new_lvec = np.array([[1,1,-2],[-1,1,0],[1,1,1]])
        lattice_rescale2 = np.array([0.5, 0.5, 1])
    else:
        sys.exit("Error: The plane you are asking for is not implemented/supported. Please use one of the below:\n 100/010/001, 110/101/011, 111\n")
        
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

#%% ------------------------------------------------------------------------------------------

if __name__=="__main__":
    ## Print help  msg if no arguments are given
    if len(sys.argv)==1:
            parser.print_help()
            sys.exit(0)
    
    ## Parse the arguments
    args = parser.parse_args()
    
    filename = args.filename
    lp_ary = np.asarray(args.lp, dtype=float)
    SuperCell = args.SC
    S_factor = np.asarray(args.SF, dtype=float)
    UNIT_CELL = args.UC
    cst = args.CS
    plane110 = args.plane
    repeat_layer = args.layer
    vacuum_layer = args.vac
    Layer_vacuum = args.Lvac
    lllp = len(repeat_layer)
    
    print(" ")
    ## Obvious error check
    if len(lp_ary) > 3:
        sys.exit("Error: Only 1 or 3 lattice parameters are allowed.\n")
    
    ## Repeation unit
    if lllp == 1:
        print("* Only repetation of %d unit along z-axis.\n"%(repeat_layer[0]))
        repeat_layer = [1, 1, repeat_layer[0]]
    elif lllp == 2:
        print("* Repetation of %d, %d unit along x and y-axis, respectively.\n"%(repeat_layer[0], repeat_layer[1]))
        repeat_layer += [1]
    elif lllp == 3:
        print("* Repetation of %d, %d, %d unit along x, y and z-axis, respectively.\n"%(repeat_layer[0], repeat_layer[1], repeat_layer[2]))
    else:
        sys.exit("Error: Only 1, 2 or 3 directions for the repeations are allowed.\n")
    
    print("* If vacuum is given, it will always be applied along z/c-axis.\n")
    
    ## Read lattice parameters
    if len(lp_ary) > 0:
        if len(lp_ary) == 1:
            lattice = lp_ary.repeat(3)
        else:
            lattice = np.copy(lp_ary)
        system = ['New Structure']
        sname = ['XX   XX']
    elif os.path.isfile(filename):
        system, sname, lattice = read_file(filename, SuperCell, S_factor, UNIT_CELL)
    else:
        sys.exit("Error: Neither lattice parameter nor any filename was supplied. Default 'POSCAR' file is also not found in current directory.\n")
    
    ## Create the final unit
    if cst == 'zincblende' or cst=='cubic' or cst=='strechedzincblende':
        lattice, posGa, posAs = rotated_new_pos_zincblende(lattice, plane110)
        lattice_new_vec = np.zeros((3, 3))
        np.fill_diagonal(lattice_new_vec, lattice)
    
        for J in range(lllp):
            unitcellGa_core = np.copy(posGa)
            unitcellAs_core = np.copy(posAs)
            for I in range(repeat_layer[J]-1):
                unitcellGa_core += lattice_new_vec[J]
                posGa = np.concatenate((posGa,unitcellGa_core), axis=0)
                
                unitcellAs_core += lattice_new_vec[J]
                posAs = np.concatenate((posAs,unitcellAs_core), axis=0)
    
        # New no. of Ga and As
        nGa = len(posGa)
        nAs = len(posAs) 
    
        if Layer_vacuum:
            vacuum_layer *= lattice[2]
    
        # New lattice vectors
        lattice *= repeat_layer
        lattice[2] += vacuum_layer
        np.fill_diagonal(lattice_new_vec, lattice)
    else:
        sys.exit("Error: The crystal you asked for is not supported in this script. It only supports onle one of the below:\n zinc blende, streched zinc blende, cubic")
             
    #%% ----------------------------------------------------------------------------------------------------
    print("\n----------------- Creating final file -----------------------")   
    ## Open output file
    nfname = filename+'_new.vasp'
    try:
        os.remove(nfname)
    except OSError:
        pass
    fnamm = open(nfname,'a')
    
    ## Save output
    np.savetxt(fnamm,system,fmt='%s')
    np.savetxt(fnamm,[1])
    np.savetxt(fnamm,lattice_new_vec)
    np.savetxt(fnamm,sname,fmt='%s')
    np.savetxt(fnamm,[[nGa, nAs]],fmt='%d')
    np.savetxt(fnamm,['Cartesian'],fmt='%s')
    np.savetxt(fnamm,posGa)
    np.savetxt(fnamm,posAs)
    
    fnamm.close()
    print("* Finish: Please change the line 6 and 7 in the '%s' file if needed. For diamond structure add the 2 species in line 6 and 7. \n"%(nfname))
