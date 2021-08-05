#!/home/bmondal/anaconda3/bin/python
import os
import time
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import argparse
import sys

parser = argparse.ArgumentParser(prog='BandStructure', description='This script plots the band structure using VASP EIGENVAL, PROCAR and POSCAR file', epilog='Have fun!')
parser.add_argument('-d', metavar='DIRECTORYNAME', default=".", help='The file path for EIGENVAL, PROCAR and POSCAR (default: current directory). e.g. /home/mondal/VASP/test/')
parser.add_argument('-EFermi', type=float, help='Fermi energy (default: Efermi from OUTCAR). All band energies will be scaled by Fermi energy such that Fermi level = at 0.0 eV.')
parser.add_argument('-EMIN', type=float, help='Minium energy during plotting (default: EMIN from DOSCAR).')
parser.add_argument('-EMAX', type=float, help='Maximum energy during plotting (default: EMAX from DOSCAR).')
parser.add_argument('-NSKPTS', type=int, default=5, help='No. of high symmetry KPOINTS (default=5). e.g 5 for L-G-X-U K-G. U and K together considered to be one.')
parser.add_argument('-NPTS', type=int, default=40, help='No. of intermediate points in each symmetry line section (default: 40)')
parser.add_argument('-NELEC', type=int, default=8, help='No. of total electron in the system (default: 8)')
parser.add_argument('-LNONCOL', action='store_true', default=False, help='LNONCOLLINEAR=True or False (default: False)')
parser.add_argument('-ISPIN', type=int, default=1, help='ISPIN value. If LNONCOLLINEAR is True this is not needed. (default: 1)')
parser.add_argument('-Dcircle', action='store_true', default=False, help='If you want to draw circle on special points (default: False)')
parser.add_argument('-OnlyCB', action='store_true', default=False, help='If you want to plot only conduction band (default: False)')
parser.add_argument('-SaveFig', action='store_true', default=False, help='If you want to save the figure (default: False)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

#if len(sys.argv)==1:
#	parser.print_help()
#	sys.exit(0)
#args = parser.parse_args()

try:
	args = parser.parse_args()
except:
	#parser.print_help()
	sys.exit(0)


print("")
print("*-*-* IT IS RECOMMENDED TO USE 'IPYTHON' FOR THE INTERACTIVE PLOTS *-*-*")
print("* Ipython recipe: \n 1. Go to the folder where you have EIGENVAL, PROCAR and POSCAR files.\n 2. $ipython \n 3. $%run ~/script/BandStructure")
print("")
#input("++++  Are you ready to go?\n")
print("")
dirname = args.d
NsymmetryKPTS = args.NSKPTS 
NPOINTS = args.NPTS
ispin = args.ISPIN
LNONCOLLINEAR = args.LNONCOL
NELECTRON = args.NELEC
drawcircle = args.Dcircle
onlycb=args.OnlyCB
savefig=args.SaveFig

poscarfile = dirname+"/POSCAR"
eigenvalfile = dirname+"/EIGENVAL"
procarfile = dirname+"/PROCAR"
doscar = dirname+"/DOSCAR"
basic_info_system = np.genfromtxt(doscar,skip_header=5,max_rows=1,dtype=float)
EMAX = basic_info_system[0]
EMIN = basic_info_system[1]
NEDOS = basic_info_system[2]
EFERMI = basic_info_system[3]
if (args.EFermi): EFERMI = args.EFermi
if (args.EMIN): EMIN = args.EMIN
if (args.EMAX): EMAX = args.EMAX

if (ispin==2) or (LNONCOLLINEAR==True): 
	VB1_index = NELECTRON - 2 # -2 because python array starts from index 0
	indexfactor = 2
else:
	VB1_index = NELECTRON//2 - 1 # -1 because python array starts from index 0
	indexfactor = 1
CB1_index = VB1_index+indexfactor
CB2_index = CB1_index+indexfactor

NKPTS = (NsymmetryKPTS - 1)*NPOINTS  # N symmetry kpoints ==> (N-1) symmetry line section

if (savefig):
	file_fig = "/home/bmondal/MyFolder/VASP/DailyReport/SupplementaryMovies/movie20/"
	file_Eg_fig = file_fig+str(dirname)+".png"
	if (os.path.exists(file_Eg_fig)):
		print("\nWarning:: Fig %s already exists"%file_Eg_fig)
		#input("Do you want to re-evaluate the fig.?\n Press 'Enter' if yes, else quit.\n")
	else:
		pass

def read_file(filen, NKPTS):
    with open(filen, 'r') as f:
    	data_list = [l for l in (line.strip() for line in f) if l]
    system_info = list(map(int,data_list[5].split()))
    #kpoints = np.array([np.array(data.split()).astype(np.float) for data in data_list[6::19]])
    #Bands = np.array([np.array(data.split()).astype(np.float) for I,data in enumerate(data_list[6:]) if not (I%(system_info[2]+1)==0)])
    skip = system_info[2]+1
    t = NKPTS*(1+system_info[2])
    kpoints, Bands= [], []
    for I,data in enumerate(data_list[-t:]):
    	arr = np.array(data.split(), dtype=float)
    	if (I%skip==0):
    		kpoints.append(arr[:4])
    	else:
    		Bands.append(arr[1])
    BandsArray = np.split(np.array(Bands),NKPTS)
    
    print ("Actual no. of Electrons, KPOINTS and Bands = ",system_info)
    print ("No. of KPOINTS read(from backward) = %d, Bands read = %d"%(len(kpoints), len(BandsArray[-1])))
    if ((len(kpoints) != NKPTS) or (len(BandsArray[-1]) != system_info[2])):
    	print ("Actual number of KPOINTS are = %d, but only %d KPOINTS has been read"%(system_info[1],len(kpoints)))
    	print ("Actual number of Bands are = %d, but only %d Bands has been read"%(system_info[2],len(BandsArray[-1])))
    	exit()
    return system_info, kpoints, BandsArray

def ReadPoscar(poscar):
	filename=open(poscar,"r")
	lines=filename.readlines()
	filename.close()
	factor=float(lines[1].split()[0])
	a=np.asarray(lines[2].split()).astype(np.float)
	b=np.asarray(lines[3].split()).astype(np.float)
	c=np.asarray(lines[4].split()).astype(np.float)
	vector=np.array([a,b,c])
	vector=vector*factor
	ion = np.array(lines[5].split(),dtype=str)
	ionnumber = np.array(lines[6].split(),dtype=int)
	return vector, ion, ionnumber

def LatticeConstant(vectors,strain=False):
	def vec(a):
	        return (np.sqrt(np.sum(a**2)))
	vector=np.array([vec(vectors[0]),vec(vectors[1]),vec(vectors[2])])
	S_factor=np.array([1., 1., 1.])
	primitive_vectors=np.divide(vector,S_factor)
	nv=primitive_vectors**2
	lp=  np.sqrt(2) * np.sqrt(nv[0] + nv[1] - nv[2])
	#lp=  np.sqrt(2) * np.sqrt(nv[0] + nv[2] - nv[1])
	#lp= np.sqrt(2) * np.sqrt(nv[1] + nv[2] - nv[0])
	if (strain):
		eqm = 5.47393 #4.55138 #GaP=5.47393 #GaAs=5.689
		print("\nThe Lattice constant used for strain calculation is = %6.3f angstrom"%eqm)
		lp = (lp - eqm)/eqm *100.

	return lp

system_info, kpoints, bands = read_file(eigenvalfile, NKPTS)

# Read the elements from POSCAR file
B, iontype, iontypenumber = ReadPoscar(poscarfile)
#reciprocal_lattice = (np.roll(np.cross(B,np.roll(B,1,axis=1)),1,axis=1))/np.linalg.det(B)
reciprocal_lattice = np.linalg.inv(B).T
print("Reciprocal lattice vectors: \n", reciprocal_lattice)
print("Elements: ", iontype)
print("* Please provide total cation and anion numbers as in order given above. \n e.g. For [Ga, As, N]=[216, 196, 20] the cation-anion array will be [216 216 216]\n")
#CA = [input("No. of cation/anion:  ") for i in range(len(iontype))]
CA = [1, 1]
print("Cation/anion array:  ",CA)
print("")
ionconc = iontypenumber/np.array(CA, dtype=float)*100.
icontitle = np.array2string(ionconc, precision=2,separator=', ',suppress_small=True)

#x = np.linspace(1,NKPTS,NKPTS)
KP = np.dot(np.array(kpoints)[:,:3],reciprocal_lattice)
vkpts = np.split(KP,NsymmetryKPTS - 1)

vkpt_diff = []; shift_x= 0; SpecialKPOINTS = [0]
for ii in range(NsymmetryKPTS-1):
	kp_distants = distance.cdist(np.array([vkpts[ii][0]]), vkpts[ii])
	vkpt_diff.append(kp_distants+shift_x)
	shift_x += kp_distants[0][-1]
	SpecialKPOINTS.append(shift_x)
x = np.array(vkpt_diff).flatten()
#print(x)
# or x = np.cumsum(np.insert(np.linalg.norm(np.diff(vkpts,axis=1),axis=2),0,0,axis=1))
# SpecialKPOINTS = np.append(x[::NKPOINTS],x[-1])

ypre = np.array(bands) - EFERMI  # Scaling of energy axis by E-fermi
y=np.copy(ypre)
EMIN = max(EMIN,np.amin(ypre))
EMAX = min(EMAX,np.amax(ypre))
#EMINband = np.searchsorted(ypre[0],EMIN)
#EMAXband = np.searchsorted(ypre[0],EMAX)
#y = ypre[:, EMINband:EMAXband]

#SpecialKPOINTS = np.linspace(1,NKPTS,NsymmetryKPTS,endpoint=True,dtype=int)

fig, ax = plt.subplots(num=1,figsize=(10,8))

CB_min = np.amin(ypre[:,CB1_index])
if (onlycb):
        # Only for conduction band plot
        # Use compression4 from csc
        SpecialKPTS = ["L","$\Lambda_{1,2}$","$\Gamma$","$\Delta_{1,2}$","X","S$_1$"]
        EMIN = CB_min-0.1; EMAX = np.amax(ypre[:,8])+0.1
        x1 = [40,120] # Position of discontinuties
        xcb = np.insert(x,x1,np.nan) # NaN will ne inserted before the positions given in x1 array
        ycb = np.insert(ypre[:,CB1_index:CB2_index],x1,[np.nan,np.nan],axis=0)
        ax.plot(xcb,ycb,linewidth=2)
        ax.set_xticks(SpecialKPOINTS)
        ax.set_xticklabels(SpecialKPTS,fontsize=28)
        ax.tick_params(axis='x',which='both',length=0,labelsize=28)
else:
	SpecialKPTS = ["L","$\Gamma$","X","U,K","$\Gamma$"]
	ax.plot(x, y, "k-",linewidth=2)
	# Plot CBM and VBM with color
	ax.plot(x,ypre[:,VB1_index:CB1_index],linewidth=2)
	ax.plot(x,ypre[:,CB1_index:CB2_index],linewidth=2)
	ax.set_xticks(SpecialKPOINTS)
	ax.set_xticklabels(SpecialKPTS,fontsize=28)

if (drawcircle):
        from matplotlib.patches import Ellipse
        # For GaAs isotropic expansion
        cc = 'r'
        #ax.add_artist(Ellipse((SpecialKPOINTS[2],ypre[2*NPOINTS][CB1_index]),0.04,0.7,color=cc,linewidth=3,fill=False))
        #ax.add_artist(Ellipse((SpecialKPOINTS[1],ypre[1*NPOINTS][VB1_index-1]),0.04,0.7,color=cc,linewidth=3,fill=False))
        # For GaAs isotropic compression
        ax.add_artist(Ellipse((SpecialKPOINTS[1],ypre[1*NPOINTS][CB2_index+1]),0.1,2.2,color=cc,linewidth=3,fill=False))

ax.set_ylabel('E (eV)',fontsize=28)
ax.set_xlabel('k-points',fontsize=28)
# Adding extra line at CBM and gamma point.
# 40: gamma k-point number;  8: CBM band number   # Based on spin-orbit calculation
GKNP = (SpecialKPTS.index("$\Gamma$")) * NPOINTS
CB_min_gamma = ypre[GKNP][CB1_index]
ax.hlines(CB_min_gamma,0,SpecialKPOINTS[-1],color='b',linestyle='-',linewidth=3)
ax.hlines(CB_min,0,SpecialKPOINTS[-1],color='m',linestyles='-',linewidth=3)

ax.vlines(SpecialKPOINTS,EMIN,EMAX,linestyles='--',color='k')
ax.tick_params(axis='y',which='both',labelsize=28)
#ax.set_title(str(iontype)+"  "+icontitle+"  "+str(dirname),fontsize=18)
ax.set_title(''.join(iontype)+",   "+"Strain = "+np.array2string(LatticeConstant(B,strain=True),precision=2,suppress_small=True)+" $\%$",fontsize=18)
#ax.set_title("GaAs,   "+"$a = $"+np.array2string(LatticeConstant(B),precision=4)+r" $\AA$",fontsize=18)
ax.set_ylim(bottom=EMIN, top=EMAX)
ax.set_xlim(left=0, right=SpecialKPOINTS[-1])
#plt.show()

#OrbitalContributionAnalysis = input("Q: Do you want to do different orbital contribution analysis? Press 'y', else 'ENTER'.\n")
OrbitalContributionAnalysis = False
if (OrbitalContributionAnalysis):
    print("------- Orbital contribution analysis -------")
    print("")
    print("* Please provide the following information about PROCAR file")
    #LORBIT = int(input("Q: What Was the LORBIT value?\n"))
    LORBIT = 12
    if (LORBIT==12) or (LORBIT==2):
        LORBITPHASE = True
    print("")
    NKPTS_total = system_info[1]
    NIONS = len(iontype)
    NBANDS = system_info[2]
    jumpk = (NKPTS_total - NKPTS)
    skip = 2+(jumpk*(1+(NBANDS*(1+2+NIONS))))
    kskip_factor = NBANDS*(3+NIONS)
    
    if (LNONCOLLINEAR): 
        print("* The POSCAR is written using LNONCOLLINEAR=True. We further want to ask -")
        #SpinDecompositionAnalysis = input("o Q: Do you want to do spin decomposition analysis? Press 'y', else 'ENTER'.\n")
        kskip_factor += (NBANDS*(1+NIONS)*3)
        skip += (jumpk*(NBANDS*3*(NIONS+1)))
    if (LORBITPHASE): 
        print("* The POSCAR is written using LORBIT=2 or >=12. We further want to ask -")
        #PhaseAnalysis = input("o Q: Do you want to do phase factor analysis? Press 'y', else 'ENTER'.\n")
        kskip_factor += (NBANDS*(2+NIONS))
        skip += (jumpk*(NBANDS*(NIONS+2)))
    
    
    def read_PROCAR(filen, NIONS, kskip_factor, I, LNONCOLLINEAR,LORBITPHASE,SDA,PA):
            with open(filen, 'r') as f:
                data = [l for l in (line.strip() for line in f) if l]
            ions, ions_total,Mx,Mx_total,My,My_total,Mz,Mz_total,phase,charge = [],[],[],[],[],[],[],[],[],[]
            orbitals = data[4].split()
            kskip = I 
            while (I < len(data)):
                if (I==kskip):
                    #print(I,data[I])
                    kskip += (kskip_factor+1)
                    I += 1
                I += 2 
                ions.append([np.array(J.split()) for J in data[I:I+NIONS]])
                I += NIONS
                ions_total.append(data[I].split()[1:])
                I += 1
                if(LNONCOLLINEAR):
                    if (SDA):
                        Mx.append([np.array(J.split()) for J in data[I:I+NIONS]])
                        I += NIONS
                        Mx_total.append(data[I].split()[1:])
                        I += 1
                        My.append([np.array(J.split()) for J in data[I:I+NIONS]])
                        I += NIONS
                        My_total.append(data[I].split()[1:])
                        I += 1
                        Mz.append([np.array(J.split()) for J in data[I:I+NIONS]])
                        I += NIONS
                        Mz_total.append(data[I].split()[1:])
                        I += 1
                    else:
                        I += (3*(1+NIONS))
                if (LORBITPHASE):
                    if (PA):
                        I += 1
                        phase.append([np.array(J.split()) for J in data[I:I+NIONS]])
                        I += NIONS
                        charge.append(data[I].split()[1:])
                        I += 1
                    else:
                        I += (2+NIONS)
            return orbitals, ions, ions_total,Mx,Mx_total,My,My_total,Mz,Mz_total,phase,charge
    
    SpinDecompositionAnalysis = False; PhaseAnalysis = False
    orbital_array, ions, ions_total,Mx,Mx_total,My,My_total,Mz,Mz_total,phase,charge = read_PROCAR(procarfile, NIONS,kskip_factor,skip,LNONCOLLINEAR,LORBITPHASE,SpinDecompositionAnalysis,PhaseAnalysis)

    ions_array = np.array(ions, dtype = float)
    ions_total_array = np.array(ions_total, dtype = float)
    #normalized = input("DO you like to normalize the orbital contributions? Press 'y', else press 'ENTER'. \n")
    normalized = True
    ions_total_array = np.insert(ions_total_array,0,0,axis=1)
    if (normalized):
        ions_array = ions_array/ions_array[:,:,-1].reshape(len(ions_array),2,1)    #ions_total_array[:,-1].reshape((len(ions_total),1,1))
        ions_total_array = ions_total_array/ions_total_array[:,-1].reshape((len(ions_total),1))
    ions_array = np.array(np.vsplit(ions_array,NKPTS))
    ions_total_array = np.array(np.vsplit(ions_total_array,NKPTS))
    
    iontypemean={}
    colIni = 0
    for I in range(len(iontype)):
        colFin = colIni + iontypenumber[I]
        iontypemean[iontype[I]] = np.mean(ions_array[:,:,colIni:colFin],axis=2)
        colIni += colFin 

    update_iontype = list(iontype)+[''.join(list(iontype))]
    iontypemean[update_iontype[-1]] = ions_total_array

    print("\nION list: ",update_iontype,"\n")
    
    print("* Please provide the following details for the orbital contributions analysis.")
    sh=np.shape(y)
    scatterx = np.reshape(np.repeat(x,sh[1]),sh)
    # Masking every 2 points in between X-U,K due large number of points in this region, which
    # makes the scatter plot looks too dense.
    maskarray = np.zeros(np.shape(scatterx),dtype=bool)
    maskarray[80:121:3,:] = True
    maskarray[81:121:3,:] = True
    #element = input("1. Which element do you want to analyze? Please use above ion_list for reference.\n")
    element = "GaAs"
    colorm = ['r','m','m','m','g','k','k','k','k','k','b','c']
    oo = ['s', 'px', 'py', 'pz', 'p', 'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2' ,'d', 'tot'] 
    while (element):
        while element not in update_iontype:
            print("Error: The provided element does not exists in the POSCAR element list. Please recheck POSCAR and retype the element.")
            element = input("1. Which element do you want to analyze?\n")
    
        #orbital = input("2. Which orbital contribution [s, (px, py, pz) or p, (dxy, dyz, dz2, dxz, x2-y2) or d, tot] do you want to analyze?\n")
        orbital = "d" 
        while(orbital):
            if(LORBIT == 10) or (LORBIT == 0):
                orbital_index = orbital_array.index(orbital)
                sizearray = iontypemean[element][:,:,orbital_index]
            else:
                if (orbital=='p'):
                    sizearray = np.sum(iontypemean[element][:,:,2:5],axis=2)
                elif (orbital=='d'):
                    sizearray = np.sum(iontypemean[element][:,:,5:10], axis=2)
                else:
                    orbital_index = orbital_array.index(orbital)
                    sizearray = iontypemean[element][:,:,orbital_index]
            s = np.copy(sizearray)
            s = np.ma.masked_array(s, mask=maskarray)
            #s = sizearray[:, EMINband:EMAXband]
            #plt.ion()
            plt.scatter(scatterx,y,50*s,color=colorm[oo.index(orbital)],label=element+'-'+orbital) 
            plt.legend(markerscale=1.,fontsize=13,loc=9,ncol=2,columnspacing=1.5)
            #plt.draw()
            #plt.show()
            #orbital = input("2. Which further orbital contribution [s, (px, py, pz) or p, (dxy, dyz, dz2, dxz, x2-y2) or d, tot)] do you want to analyze? If not press 'ENTER'.\n")
            orbital = False
        #element = input("* Hi! Do you want to check for new elements? Press 'y' input new element. Else press 'Enter'.\n")
        element = False
    if (SpinDecompositionAnalysis):
        Mx_array = np.array(Mx, dtype=float)
        My_array = np.array(Mx, dtype=float)
        Mz_array = np.array(Mx, dtype=float)
        ions_total_array = np.array(ions_total, dtype=float)
        Mx_total_array = np.array(Mx_total, dtype=float)
        My_total_array = np.array(My_total, dtype=float)
        Mz_total_array = np.array(Mz_total, dtype=float)
        print("Warning: Spin decomposition analysis implementation is not comlpete yet.\n")
    if (PhaseAnalysis):
        phase_array = np.array(phase,dtype=float)
        charge_array = np.array(charge,dtype=float)
        print("Warning: Phase analysis implementation is not comlpete yet.\n")

print("\n-------- DOS analysis ----------\n")
#DOSanalysis = input("Do you want to do DOS analysis?\n ")
DOSanalysis = False
if (DOSanalysis):
    def read_DOSCAR(filen,NION):
        DosAll = np.genfromtxt(filen, dtype=float,skip_header=6,unpack=False,max_rows=NEDOS)
        pdos = {}
        readskip = 6+NEDOS
        for I in range(NION):
            pdos[str(I)] = np.genfromtxt(filen,dtype=float,skip_header=int(readskip+1), unpack=False,max_rows=NEDOS)
            readskip += NEDOS+1
        return DosAll,pdos
    
    DOS,pDOS = read_DOSCAR(doscar,len(iontype)) 
    XX=DOS[:,0]
    orb = ["","s","p","d"]
    orbitals=["","s","py","pz","px","dxy","dyz","dz2","dxz","dx2-y2","p","d"]
    orbitalsnoncollinearl=["","s","s_mx","s_my","s_mz","p","p_mx","p_my","p_mz",
            "d","d_mx","d_my","d_mz"]
    orbitalsnoncollinearlm=["","s","s_mx","s_my","s_mz","py","py_mx","py_my","py_mz",
            "pz","pz_mx","pz_my","pz_mz","px","px_mx","px_my","px_mz",
            "dxy","dxy_mx","dxy_my","dxy_mz","dyz","dyz_mx","dyz_my","dyz_mz",
            "dz2","dz2_mx","dz2_my","dz2_mz","dxz","dxz_mx","dxz_my","dxz_mz",
            "dx2-y2","dx2-y2_mx","dx2-y2_my","dx2-y2_mz","p","d"]
    print("Here is your list of orbitals:")
    print("Collinear, ISPIN=1/2, LORBIT=0/5/10 :\n",orb[1:])
    print("Collinear, ISPIN=1/2, LORBIT=1/2/>=11 :\n",orbitals[1:])
    print("Noncollinear, LORBIT=0/5/10 :\n",orbitalsnoncollinearl[1:])
    print("Noncollinear, LORBIT=1/2/>=11 :\n",orbitalsnoncollinearlm[1:])
    
    fig, ax=plt.subplots()
    ax.set_xlabel("Energy (eV)", fontsize=18)
    ax.set_ylabel("DOS",fontsize=18)
    ax.plot(XX,DOS[:,1])
    partialDOS = input("\nDo you want to do projected DOS analysis? If not press 'ENTER'.\n")
    if (partialDOS):
        if not OrbitalContributionAnalysis:
            LORBIT = int(input("Q: What was the LORBIT value?\n"))
        if not LNONCOLLINEAR:
            ISPIN = int(input("Q: What was the ISPIN value?\n"))
        ep = input("\nWhich element do you want to analyze for?\n")
        while (ep):
            e = list(iontype).index(ep)
            orbit = input("Which orbital projected DOS do you want to plot?\n")
            plt.ion()
            while(orbit):
                if (LNONCOLLINEAR):
                        if (LORBIT>=11) or (LORBIT==2) or (LORBIT==1):
                            if (orbit=="p"):
                                y=np.sum(pDOS[str(e)][:,5:17:4],axis=1)
                                ax.plot(XX,y,label=ep+"-"+"p")
                            elif (orbit=="d"):
                                y=np.sum(pDOS[str(e)][:,17::4],axis=1)
                                ax.plot(XX,y,label=ep+"-"+"d")
                            else:
                                yind = orbitalsnoncollinearlm.index(orbit)
                                ax.plot(XX, pDOS[str(e)][:,yind],label=ep+"-"+orbit)
                        else:
                                yind = orbitalsnoncollinearl.index(orbit)
                                ax.plot(XX, pDOS[str(e)][:,yind],label=ep+"-"+orbit)
                else:
                    if(ISPIN==2):
                        spin = input("which spin do you want to examine? If spin up insert 0, for spin down insert 1.\n")
                        if (LORBIT>=11) or (LORBIT==2) or (LORBIT==1):
                            if (orbit=="p"):
                                start = 3+int(spin)
                                y=np.sum(pDOS[str(e)][:,start:9:2],axis=1)
                                ax.plot(XX,y,label=ep+"-p"+"-"+spin)
                            elif (orbit=="d"):
                                start = 11+int(spin)
                                y=np.sum(pDOS[str(e)][:,start::2],axis=1)
                                ax.plot(XX,y,label=ep+"-d"+"-"+spin)
                            else:
                                obindex = orbitals.index(orbit)*2-1+int(spin)
                                ax.plot(XX,pDOS[str(e)][:,obindex],label=ep+"-"+orbit+"-"+spin)
                        else:
                            obindex = orb.index(orbit)*2-1+int(spin)
                            ax.plot(XX,pDOS[str(e)][:,obindex],label=ep+"-"+orbit+"-"+spin)
                    else:
                        if (LORBIT>=11) or (LORBIT==2) or (LORBIT==1):
                            if (orbit=="p"):
                                y=np.sum(pDOS[str(e)][:,2:5],axis=1)
                                ax.plot(XX,y,label=ep+"-"+"p")
                            elif (orbit=="d"):
                                y=np.sum(pDOS[str(e)][:,5:10],axis=1)
                                ax.plot(XX,y,label=ep+"-"+"d")
                            else:
                                obt = orbitals.index(orbit)
                                ax.plot(XX,pDOS[str(e)][:,obt],label=ep+"-"+orbit)
                        else:
                            obt = orb.index(orbit)
                            if (orbit) not in orb:
                                print("The data format is not matching. Please check.")
                                print("You have choosen LNONCOLLINEAR=FALSE, ISPIN=1,LORBIT=0/5/10")
                                exit()
                            ax.plot(XX,pDOS[str(e)][:,obt],label=ep+"-"+orbit)
                plt.legend(markerscale=2.,fontsize=12,loc=9,ncol=2,columnspacing=1.5)
                plt.draw()
                orbit = input("DO you want to check for other orbital projections for same element? If yes enter new orbital. Else press 'ENTER'\n")
            ep = input("Do you want to check for next element? If yes enter new element. Else press 'ENTER'\n")

# Save the figure and show
plt.tight_layout()
#plt.legend()
if(savefig):
	plt.savefig(file_Eg_fig,format='png',dpi=300)
else:
	plt.show()

