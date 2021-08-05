#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 15:21:18 2021

@author: bmondal
"""

import sys
import numpy as np      
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy.optimize import curve_fit

np.set_printoptions(precision=3,suppress=True)

params = {'legend.fontsize': 22,
          'figure.figsize': (8, 6),
         'axes.labelsize': 24,
         'axes.titlesize': 24,
         'xtick.labelsize': 24,
         'ytick.labelsize': 24,
         'errorbar.capsize':2}
plt.rcParams.update(params)



#%%
dirname = "/home/bmondal/MyFolder/VASP/GaAsP/mcsqs/data"
savefigfile="/home/bmondal/MyFolder/VASP/GaAsP/mcsqs/figs"

#-------------------------------------------------------------------------------
# We strated with 10 sqs structure each of which has different N-
# arrangement. Under equilibration individual 10 sqs might ended up with each having different eqm.
# lattice parameter. But remember, as they are equlibrium structures and hence, every signle of them
# must be strain free. Then with strain every structure by the same amount such that 'each' structure
# would possess x-amount of strain. This implies there might be stderr in 'lattice parameter' for each
# set of strain sqs, but there shouldn't be any stderr in 'strain'.
# If we calculate strain in the following way, 1st calculate average lattice parameter in each strain set
# and then calculate average strain from average lattice parameter; that would indirectly imply different
# sqs in a single set would possess different strain. But that's not the case in our calculation.
# On top of that, imagine someone keep the lattice constant fix for all 10 sqs structures in a set
# and then clealy, there will be no uncertainty in lattice parameter but there will be uncertainty in
# strain measurement. The disadvantage of this procedure is that there is an inherent dependency of
# strain in bandgap at each average-strain point which wouldn't be possible to quantify.
'''
#*******************************************************************************
## Wrong way::
    
with open(filename+"lattice_param.dat", "r") as f:
    data = [np.asarray(l.split(), dtype=float) for l in (line.strip() for line in f) if l]
lp = [np.split(data[i][1:],metadata[i]) for i in range(len(data))]

average_lp, stderr_lp = [], []
for Eg in lp:
    bg = []
    for Egg in Eg:
        bg.append(np.mean(Egg))
    average_lp.append(bg)

    
**calculate ave_strain and stderr from average_lp
#*******************************************************************************
## Correct way::
    # But we don't need to do this calculation because- 1. we already know the strain
    # informations when we submit the vasp run. 2. The strain calculation in the for-
    # loop wouldn't work for cases where some the strauctures in strain set didn't 
    # converge. For e.g. In our case for N conc 9.3% all the 10 equlibrium sqs converged
    # but at 5% strain only 3 of them converged. So the array size wouldn't match.
    # Also in metadata.dat we don't have the information which sqs's didn't converge.
    # So, we shouldn't compare sqs1 eqm. with sqs10 at 5% strain.
    # But anyway just for few check perposes one can use this.

average_st, stderr_st = [], []
for Eg in lp[0:3]:
    bg, sg = [], []
    eqm = Eg[0]
    for Egg in Eg:
        strain = (Egg - eqm)/eqm * 100
        bg.append(np.mean(strain))
        sg.append(np.std(strain))
    average_st.append(bg)
    stderr_st.append(sg)

print(np.array(average_st))
print(np.array(stderr_st))
'''
#-------------------------------------------------------------------------------
def createdata(filename):
    with open(filename+"bandgap.dat", "r") as f:
        data = np.array([np.asarray(l.split("|"), dtype=str) for l in (line.strip() for line in f) if l], dtype=str)
    
    with open(filename+"GetBandgapSupercell.out", "r") as f:
        metadata1 = np.array([np.asarray(l, dtype=str) for l in (line.strip() for line in f) if l.startswith("c") ], dtype=str)
        f.seek(0)
        strainmetadata = np.array([np.asarray(l.split("/")[-2], dtype=str) for l in (line.strip() for line in f) if l.endswith("N00.9") ], dtype=str)
    
    strain_raw = np.core.defchararray.replace(strainmetadata,"S","")
    strain = np.array(np.core.defchararray.replace(strain_raw,"Eqm","0"), dtype=float)
    NC = np.array(data[:,0],dtype=float)
    
    metadata2 = np.split(metadata1,len(NC))
    
    average_bandgap, stderr_bandgap, count, NearestConfig,NearestConfigBG = [], [], [], [],[]
    for (ii,i) in enumerate(data):
        bg, sg, c, Nconfig, NconfigBG = [], [], [NC[ii]], [NC[ii]], [NC[ii]]
        metadata = [metad.split("  ") for metad in metadata2[ii]] # Note: split uses 2 spaces
        for (jj,j) in enumerate(i[1:-1]):
            Egg = np.array(j.split(",")[:-1], dtype=np.float)
            c.append(len(Egg))
            if not len(Egg): Egg = np.nan
            Eggmean = np.mean(Egg)
            mpos = np.argmin(np.abs(Egg-Eggmean))
            Nconfig.append(metadata[jj][mpos])
            NconfigBG.append(Egg[mpos])           
            bg.append(Eggmean)
            sg.append(np.std(Egg))   
        average_bandgap.append(bg)
        stderr_bandgap.append(sg)
        count.append(c)
        NearestConfig.append(Nconfig)
        NearestConfigBG.append(NconfigBG)
    return NC, average_bandgap, stderr_bandgap, strain, count, NearestConfig, NearestConfigBG

def writeNearestConfig(filename, XX):
    with open(filename+"NearestConfig.txt", "w") as f:
        np.savetxt(f, XX, fmt="%s") 
        
def ReadNearestConfigAfterBandgap(filename):
    data = np.genfromtxt(filename+"NearestConfigAfterBandgap.dat",dtype=str,delimiter="|")
    NC = np.array(data[:,0],dtype=float)
    
    ddt = []
    for line in data[:,1:-1]:
        dd =np.array([l.split(" ") for l in line], dtype=float)
        ddt.append(dd)
    restdata = np.array(ddt)
    return NC, restdata

def ReadBlochWeightSummary(filename):
    data = np.genfromtxt(filename+"BlochWeightGLDX.dat", dtype=str,delimiter="|")
    
    restdata = []
    for line in data[:,1:-1]:
        dd =[h[:-1] for h in [l.split(";") for l in line]]
        ddt = []
        for line in dd:
            ddtt = [l.split(",") for l in line]
            ddt.append(ddtt)
        restdata.append(np.array(ddt, dtype=np.float))
    return restdata

def ReadBlochWeightEnergy(filename):
    data = np.genfromtxt(filename+"BlochWeightEnergyGLDX.dat", dtype=str,delimiter="|")
    
    restdata = []
    for line in data[:,1:-1]:
        dd =[h[:-1] for h in [l.split(";") for l in line]]
        ddt = []
        for line in dd:
            ddtt = [l.split(",") for l in line]
            ddt.append(ddtt)
        restdata.append(np.array(ddt, dtype=np.float))
    return restdata

def ReadSPPROCAR(filename):
    data = np.genfromtxt(filename+"s-p-PROCAR.dat", dtype=str,delimiter="|")
    farray = []
    for line in data:
        arraym = []
        for l in line[1:-1]:
            arraym.append(np.asarray(l.split(" "), dtype=float))
        farray.append(arraym)
    
    farray = np.array(farray)
    return farray    
        
def ReadHighResBandgap(filename):
    data = np.genfromtxt(filename+"HighRes-Bandgap.dat", dtype=str,delimiter="|")
    dd = []
    for line in data:
        dd.append([l.split("/") for l in line[1:]])
    return np.array(dd, dtype=float)

def ReadHighResBW(filename):
    data = np.genfromtxt(filename+"HighRes-BlochWeightGLDX.dat", dtype=str,delimiter="|")
    dd = []
    for line in data:
        dd.append([l.split("/") for l in line[1:]])
    alld = np.array(dd) 
    strainp = np.array(alld[:,:,0], dtype=float)
    ddd = []
    for line in alld[:,:,1]:
        l = [lin.split(";")[:-1] for lin in line]
        dddd = []
        for ll in l:
            dddd.append([k.split(",") for k in ll])
        ddd.append(dddd)
    return strainp, np.array(ddd, dtype=float)

def ReadHighResBWEnergy(filename):
    data = np.genfromtxt(filename+"HighRes-BlochWeightEnergyGLDX.dat", dtype=str,delimiter="|")
    dd = []
    for line in data:
        dd.append([l.split("/") for l in line[1:]])
    alld = np.array(dd) 
    ddd = []
    for line in alld[:,:,1]:
        l = [lin.split(";")[:-1] for lin in line]
        dddd = []
        for ll in l:
            dddd.append([k.split(",") for k in ll])
        ddd.append(dddd)
    return np.array(ddd, dtype=float)

#%%
#pure GaAs isotropic case
filename = "/home/bmondal/MyFolder/VASP/GaAs/GaAscompression.dat"
GaAsdata = np.loadtxt(filename)
compression_lp = GaAsdata[:,0]
compression_Eg =  GaAsdata[:,-1]
compression_strain = (compression_lp - compression_lp[0]) /compression_lp[0]*100.
comp_strain5 = np.argwhere(compression_strain >= -5).flatten()
compression_strain = compression_strain[comp_strain5]
compression_Eg = compression_Eg[comp_strain5]
pos = np.argsort(compression_strain)
compression_strain = compression_strain[pos]
compression_Eg = compression_Eg[pos]


filename = "/home/bmondal/MyFolder/VASP/GaAs/GaAsexpansion.dat"
GaAsdata = np.loadtxt(filename)
expansion_lp = GaAsdata[:,0]
expansion_Eg =  GaAsdata[:,-1]
expansion_strain = (expansion_lp - expansion_lp[0]) /expansion_lp[0]*100.
expa_strain5 = np.argwhere(expansion_strain <= 5).flatten()
expansion_strain = expansion_strain[expa_strain5]
expansion_Eg = expansion_Eg[expa_strain5]
pos = np.argsort(expansion_strain)
expansion_strain = expansion_strain[pos]
expansion_Eg = expansion_Eg[pos]

GaAs_strain = np.concatenate((compression_strain,expansion_strain))
GaAs_Eg = np.concatenate((compression_Eg,expansion_Eg))

#%% ============================================================================
filename = dirname+"/ThermExComp"

##stage1 # Strain vs band gap
stage2 = True; spprocar=False; highstrain=0 # Blochweight and orbital contribution and High strain
plotallBWs =0; plotCBsVBs = 0; Join=1; GaAsP_P37_BW3=0
GaAsP_P37_DIT1 = 1; GaAsP_P37_DIT2=0; allbandgaphighresolution=0
stage3 = 1 # High resolution

if stage3:
    stage2 = True
if GaAsP_P37_DIT1 or GaAsP_P37_DIT2:
    stage3=1

#*******************************************************************************

# Strain vs band gap
NC, average_bandgap, stderr_bandgap, strain, count, NearestConfig, NearestConfigBG = createdata(filename)
print("SQS summary:\nN%\t#sqs(S",strain,')\n')
print(count)
print(NearestConfig)
#***************************************************************************
# Saveing only compressive strain Nearest configurations
# For GaAsN/S-5/N07.4 2nd closest sqs is conf05
compre_strin_id1 = np.argwhere(strain<=0.001).flatten() + 1 #add 1 to shift by 1 index for N% in NearestConfig
print("Compressive strain sqs final order after re-ordering:\n",strain[compre_strin_id1-1])
compre_strin_id = np.insert(compre_strin_id1,0,0)
writeNearestConfig(filename, np.array(NearestConfig)[:,compre_strin_id])
ReargNearestConfigBG = np.array(NearestConfigBG)[:,compre_strin_id]
print("Nearest neighbour configurations bandgap were:\n", ReargNearestConfigBG)

#***************************************************************************
# Linear fit for isotropic tensile strain
fit = False
#%%
#*******************************************************************************
if stage2:
    # Compare the band gap after reevaluation of neares neighbour
    # Just for checking if due to reevalution how much band gap changes. This shouldn't be > 0.005 (or 5 meV)
    _, restdata = ReadNearestConfigAfterBandgap(filename)
    NconfigBG = restdata[:,:,-1]
    print("Diff. b/w NearestConfig Eg Before & After :\n", (ReargNearestConfigBG[:,1:] - NconfigBG))
    #***************************************************************************
    # Bloch weight
    BlochWights = ReadBlochWeightSummary(filename)
    #***************************************************************************
    # Energies of Gamma and D point
    BlochWightGDEnergy = ReadBlochWeightEnergy(filename)
    #***************************************************************************
    # S-P orbital contributions from PROCAR
    if spprocar:
        spprocar = ReadSPPROCAR(filename)
    #***************************************************************************
    # > 5% strain
    if highstrain:
        filename = dirname+"/HighStrain-"
        highSt_NC, highSt_average_bandgap, highSt_stderr_bandgap, highSt_strain, \
            highSt_count, highSt_NearestConfig, highSt_NearestConfigBG = createdata(filename)
        writeNearestConfig(filename, np.array(highSt_NearestConfig))
        ReargNearestConfigBG = np.array(highSt_NearestConfigBG)
        print("Nearest neighbour configurations bandgap were:\n", ReargNearestConfigBG)
        
#*******************************************************************************
if stage3:
    # High resolution data
    HR_all = ReadHighResBandgap(filename)
    HR_Eg = HR_all[:,:,1]
    HR_strain = HR_all[:,:,0]
    
    HR_strain_BW, HR_BW = ReadHighResBW(filename)
    HR_BW_Energy = ReadHighResBWEnergy(filename)

#-------------------------------------------------------------------------------
#%%
# Sorting strain for getting the plots better in line mode
strain_sort_pos = np.argsort(strain)
mystrain = strain[strain_sort_pos]
bwxx = (-1)*mystrain[10:]

titletxt = "Thermal expansion"
xlabeltxt = "Strain (%)"
#*******************************************
fig, ax = plt.subplots()
#ax.set_title(titletxt)
ax.set_xlabel(xlabeltxt)
ax.set_ylabel(r"E$_{\mathrm{g}}$ (eV)")
ax.set_xlim(-5,5)
ax.set_ylim(0.2,2.1)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
##ax.xaxis.set_ticks(strain)
ax.axvline(color='k', ls='--')
#********************************************
fig2, ax2 = plt.subplots()
#ax2.set_title(titletxt)
ax2.set_xlabel(xlabeltxt)
ax2.set_ylabel(r"Rate=E$_{\mathrm{g}}$(eV)/strain(%)")
#ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
##ax2.xaxis.set_ticks(strain)
ax2.axvline(color='k', ls='--')
#**********************************************
if stage2:
    fig3, ax3 = plt.subplots()
    #ax3.set_title(titletxt+r"[BW=($\Gamma$:L:X::$\Delta_{\mathrm{m}}$:$\Delta$E($\Delta_{\mathrm{m}}^f-\Gamma ^f)$)]")
    #ax3.set_title(r"BW=($\Gamma$:L:X::$\Delta_{\mathrm{m}}$:$\Delta$E($\Delta_{\mathrm{m}}^f-\Gamma ^f)$)")
    ax3.set_xlim(-5.0,0.0)
    ax3.set_xlabel(xlabeltxt)
    ax3.set_ylabel(r"E$_{\mathrm{g}}$ (eV)")
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    #ax3.yaxis.set_minor_formatter(FormatStrFormatter(''))
    ax3.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ##ax3.xaxis.set_ticks(strain)
    #************************************************
    if plotallBWs:
        norow=2; nocol=5
        fig4, ax4 = plt.subplots(nrows=norow,ncols=nocol,sharey=False)
        ax4[1,0].set_xlabel(xlabeltxt)
        ax4[0,0].set_ylabel("Bloch weight (%)")
        ax4[1,0].set_ylabel("Bloch weight (%)")
    #********************************************
    if spprocar:  
        fig5, ax5 = plt.subplots(nrows=norow,ncols=nocol,sharey=False)
        fig5.suptitle("From VASP projection")
        ax5[1,0].set_xlabel("Strain (%)")
        ax5[0,0].set_ylabel("Orbital projection (%)")
        ax5[1,0].set_ylabel("Orbital projection (%)")
    #********************************************   
    if plotCBsVBs:
        fig6, ax6 = plt.subplots()
        ax6.set_title(titletxt)
        ax6.set_xlabel(xlabeltxt)
        ax6.set_ylabel(r"E$_{\mathrm{CB}}$ (eV)")
        #ax6.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        #ax6.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ##ax7.xaxis.set_ticks(strain)
        #********************************************
        fig7, ax7 = plt.subplots()
        ax7.set_title(titletxt)
        ax7.set_xlabel(xlabeltxt)
        ax7.set_ylabel(r"E$_{\mathrm{VB}}$ (eV)")
        #ax7.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        #ax7.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ##ax7.xaxis.set_ticks(strain)
#********************************************
if fit:
    fig8, ax8 = plt.subplots()
    ax8.set_title(titletxt)
    ax8.set_xlabel("Tensile strain (%)")
    ax8.set_ylabel(r"E$_{\mathrm{g}}$ (eV)")
    #ax8.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    #ax8.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ##ax8.xaxis.set_ticks(strain)
#********************************************
#plot GaAs
#ax.plot(GaAs_strain, GaAs_Eg,marker="o",ls="-.",c='k') #, label="0.00")
#**********************************************

frompos = [10,10,10,10,10,10,10,10,10,10];topos = [21,21,21,21,21,21,21,21,21,21]
bk_x, bk_y = [], []
TS_x, TS_y, TS_nc = [GaAs_strain[9]], [GaAs_Eg[9]], [0.0]
SMT_x, SMT_y = [0], [6.67]
for ii in range(4,5): #range(0,len(NC)):
    
    #---------------------------------------------------------------------------
    # Plot strain vs bandgap
    yy = np.array(average_bandgap[ii])[strain_sort_pos]
    yyerr = np.array(stderr_bandgap[ii])[strain_sort_pos]
    if NC[ii]<=0.5:
        ax.plot(mystrain, yy, marker='o', label = '{:.1f}'.format(NC[ii]))
    else:
        ax.errorbar(mystrain, yy, yerr=yyerr, marker='o', label = '{:.1f}'.format(NC[ii]))
    
    if (fit):
        # Linear fit tensile strain
        xx = mystrain[frompos[ii]:topos[ii]]
        lyy = yy[frompos[ii]:topos[ii]]
        A = np.vstack([xx, np.ones(len(xx))]).T
        m, c = np.linalg.lstsq(A, lyy, rcond=None)[0]
        n_xx = np.linspace(0,-(c/m),100)
        yyy = m*np.array(n_xx) + c
        SMT_x.append(NC[ii])
        SMT_y.append(-(c/m))
        
        if NC[ii]<=0.5:
            pp = ax8.plot(xx, lyy, marker='o', label = '{:.1f}'.format(NC[ii]))
        else:
            pp = ax8.errorbar(xx, lyy, yerr=yyerr[frompos[ii]:topos[ii]], marker='o', label = '{:.1f}'.format(NC[ii]))
        ax8.plot(n_xx,yyy, color=pp[0].get_color())
        ax8.legend(loc=1, ncol=2, labelspacing=0.3, columnspacing=0.8, handlelength=0.7,handletextpad=0.5)
    
    ax.legend(loc=1, ncol=2, labelspacing=0.3, columnspacing=0.8, handlelength=0.7,handletextpad=0.5)
    
    #-----------------------------------------------------------------------
    # Plot derivative
    yy_deriv = abs(np.diff(yy))/0.5 # 0.5=strain1-strain2=Delta_strain
    xx_deriv = np.insert(mystrain,11,0.)[1:-1] #-4.5 to +4.5 strain range
    
    ax2.plot(xx_deriv, yy_deriv, marker='o', label = '{:.1f}'.format(NC[ii]))
    ax2.legend(loc=0, ncol=2, labelspacing=0.3, columnspacing=0.8, handlelength=0.7,handletextpad=0.5)
    plt.tight_layout()
    #-----------------------------------------------------------------------
    # Calculate the breake points
    break_point_pos = np.argmax(yy_deriv)  
    bk_x.append(mystrain[break_point_pos])
    bk_y.append(yy[break_point_pos])
    #-----------------------------------------------------------------------
    # Nearest neighbours
    #ax.scatter(strain,NearestConfigBG[ii][1:],marker='*',color='k')
    #---------------------------------------------------------------------------
    
    if stage2:
        # Compressive strain Bloch weight  
        BW_x = mystrain[:11]
        BW_y = yy[:11]
        if NC[ii]<=0.5:
            ax3.plot(BW_x, BW_y, marker='o', ls='-',label = '{:.1f}'.format(NC[ii]))
        else:
            ax3.errorbar(BW_x, BW_y, yerr=yyerr[:11], marker='o', ls ='-',label = '{:.1f}'.format(NC[ii])) 
            
        ax3.legend(loc=2, ncol=2, labelspacing=0.3, columnspacing=0.8, handlelength=0.7,handletextpad=0.5)
        #-----------------------------------------------------------------------------
        # Compressive strain Bloch weight 2
        J = np.round(BlochWights[ii] * 100.).astype(int) 
        JE = BlochWightGDEnergy[ii] 
        sdone=1
        ## Spinor 1
        G = J[:,0,0]; L1 = J[:,0,1]; L2 = J[:,0,2]; L3 = J[:,0,3]; L4 = J[:,0,4]
        X1 = J[:,0,5]; X2 = J[:,0,6]; X3 = J[:,0,7]
        # L = np.amax(J[:,0,1:5], axis=1)
        # X = np.amax(J[:,0,5:8], axis=1)
        L = np.round(np.mean(BlochWights[ii][:,0,1:5], axis=1)* 100.).astype(int)
        X = np.round(np.mean(BlochWights[ii][:,0,5:8], axis=1)* 100.).astype(int)
        D = J[:,0,8]
        
        JEspinor1 = np.diff(JE[:,0])
        ## Spinor 2
        # G_2 = J[:,1,0]; L1_1 = J[:,1,1]; L2_2 = J[:,1,2]; L3_2 = J[:,1,3]; L4_2 = J[:,1,4]
        # D_2 = J[:,1,8]; X1_2 = J[:,1,5]; X2_2 = J[:,1,6]; X3_2 = J[:,1,7]
        # JEspinor2 = np.diff(JE[:,1])
        for i in range(len(bwxx)):
            ## Spinor 1
            # text=str(G[i])+":"+str(L1[i])+","+str(L2[i])+","+str(L3[i])+","+str(L4[i])+":" \
            #     +str(X1[i])+','+str(X2[i])+","+str(X3[i])+"::"+str(D[i])+":"+'{:0.2f}'.format(JEspinor1[i][0])
            text=str(G[i])+":"+str(L[i])+":"+str(X[i]) #+"::"+str(D[i])+":"+'{:0.2f}'.format(JEspinor1[i][0])
            # Condition: 1. Cut off energy diff between G and D is considered 1 meV
            #            2. D bloch weight needs to be >20% 
            if (sdone and ((G[i]<=L[i] or G[i]<=X[i]) or ((JEspinor1[i][0]<1e-3) and D[i]>20))): 
                TS_x.append(bwxx[i])
                TS_y.append(BW_y[len(bwxx)-1-i])
                TS_nc.append(NC[ii])
                sdone=0
            # Text positions are shifted by 1e-3 in x-dir and 1e-5 in y-direction
            if i<2:
                changepos=-1.08
            else:
                changepos=0
            if i<5:
                changeposy=+0.008
            else:
                changeposy=0
            #ax3.text(bwxx[i]+5e-2+changepos, BW_y[len(bwxx)-1-i]-1e-2+changeposy,s="("+text+")")
            #ax3.text(bwxx[i]+1e-3, BW_y[len(bwxx)-1-i]-1e-2,s="("+text+")",rotation=15)
            #ax3.text(bwxx[i]+1e-2, BW_y[len(bwxx)-1-i],s="("+text+")",rotation=15)
        
        if GaAsP_P37_BW3:
            ax3.axvline(x=-1.25,ymax=0.75,color='r',ls='--')
            ax3.axvline(x=-2.25,ymax=0.95,color='r',ls='--')
            ax3.axvline(x=-4.8,ymax=0.6,color='r',ls='--')
            tx = [-0.8,-1.8,-3.5,-4.98]; ty=[1.6,1.6,1.6,1.6]
            ts=['$\Gamma$','L','$\Delta_{\mathrm{m}}$','X']
            for I in range(len(ts)):
                ax3.text(tx[I],ty[I],ts[I],color='r',size=24)
    
        if plotallBWs:
            ind1=(ii//nocol); ind2=(ii%nocol)
            ax4[ind1,ind2].set_title('N% = {:.1f}'.format(NC[ii]))
            ## Spinor 1
            ax4[ind1,ind2].plot(bwxx,G,"o-",label="$\Gamma$")
            ax4[ind1,ind2].plot(bwxx,L,"o-",label="L")
            ax4[ind1,ind2].plot(bwxx,D,"o-",label="$\Delta_{m}$")
            ax4[ind1,ind2].plot(bwxx,X,"o-",label="X")
            ## Spinor 2
            # ax4.plot(bwxx,G2,"o--")
            # ax4.plot(bwxx,L2,"o--")
            # ax4.plot(bwxx,D2,"o--")
            # ax4.plot(bwxx,X2,"o--")
            ax4[ind1,ind2].legend(loc=9, ncol=2, labelspacing=0.3, columnspacing=1, handlelength=0.7,handletextpad=0.5)
        #-----------------------------------------------------------------------
        if spprocar:
            # Orbital projection
            s = spprocar[ii][:,0]
            p = np.sum(spprocar[ii][:,1:4], axis=1)
            tot = spprocar[ii][:,-1]
            pPLUSd = tot - s
            
            per_s = s/tot*100.0
            per_p = p/tot*100.0
            per_pPLUSD = pPLUSd/tot*100.
            
            ax5[ind1,ind2].set_title('N% = {:.1f}'.format(NC[ii]))
            ax5[ind1,ind2].plot(bwxx, per_s, "o-", label='s')
            ax5[ind1,ind2].plot(bwxx, per_p, "o-", label='p')
            ax5[ind1,ind2].plot(bwxx, per_pPLUSD, "o-", label="p+d")
            ax5[ind1,ind2].legend(loc=9, ncol=2, labelspacing=0.3, columnspacing=1, handlelength=0.7,handletextpad=0.5)

        if highstrain:
            ax.errorbar(highSt_strain, highSt_average_bandgap[ii], yerr=highSt_stderr_bandgap[ii], marker='o', label = '{:.2f}'.format(NC[ii]))
            
        #-----------------------------------------------------------------------
        if plotCBsVBs:
            ax6.plot(bwxx, restdata[ii,:,5], marker='o', label = '{:.1f}'.format(NC[ii]))
            ax6.legend(loc=1, ncol=2, labelspacing=0.3, columnspacing=1, handlelength=0.7,handletextpad=0.5)
            #-----------------------------------------------------------------------
            ax7.plot(bwxx, restdata[ii,:,3], marker='o', label = '{:.1f}'.format(NC[ii]))
            ax7.legend(loc=1, ncol=2, labelspacing=0.3, columnspacing=1, handlelength=0.7,handletextpad=0.5)
        
#---------------------------------------------------------------
if stage2 and Join:
    # Join the breake points
    ax.plot(bk_x[1:], bk_y[1:], ls="dotted",lw=2, color='k')
    #---------------------------------------------------------------
    # Join the Transition points
    ax.plot(TS_x, TS_y, ls="dotted",lw=2, color='r')
    #ax3.plot(TS_x[1:], TS_y[1:], ls="dotted",lw=2, color='r')
#--------------------------------------------------------------
if stage3:
    # Compressive strain Bloch weight high resolution
    TS_x, TS_y = [GaAs_strain[9]], [GaAs_Eg[9]]
    for ii in range(4,5): #range(0,len(NC)):
        #ax.plot(HR_strain[ii], HR_Eg[ii], marker='*', c='r', ls='')
        #ax3.plot(HR_strain[ii], HR_Eg[ii], marker='*',ls='-')
        ax3.plot(HR_strain[ii], HR_Eg[ii], marker='*', c='r', ls='')
        
        #-------------------------------------------------------------------
        # Compressive strain Bloch weight 2
        J = np.round(HR_BW[ii] * 100.).astype(int)
        JE = HR_BW_Energy[ii]
        sdone=1
        ## Spinor 1
        G = J[:,0,0]; L1 = J[:,0,1]; L2 = J[:,0,2]; L3 = J[:,0,3]; L4 = J[:,0,4]
        X1 = J[:,0,5]; X2 = J[:,0,6]; X3 = J[:,0,7]
        L = np.amax(J[:,0,1:5], axis=1)
        X = np.amax(J[:,0,5:8], axis=1)
        #L = np.round(np.mean(HR_BW[ii][:,0,1:5], axis=1)* 100.).astype(int)
        #X = np.round(np.mean(HR_BW[ii][:,0,5:8], axis=1)* 100.).astype(int)
        D = J[:,0,8]
        JEspinor1 = np.diff(JE[:,0])
        bwxx=HR_strain_BW[ii]
        if GaAsP_P37_DIT1:
            changepos=[0,0,-1.1,0,0,0];changeposy=[-6e-3,-6e-3,-7e-3,-5e-3,0,5e-3]
        else:
            changepos = changeposy = [0]*6
        for i in range(len(bwxx)):
            ## Spinor 1
            #text=str(G[i])+":"+str(L1[i])+","+str(L2[i])+","+str(L3[i])+","+str(L4[i])+":" \
            #    +str(X1[i])+','+str(X2[i])+","+str(X3[i])+"::"+str(D[i])+":"+'{:0.2f}'.format(JEspinor1[i][0])
            text=str(G[i])+":"+str(L[i])+":"+str(X[i])+"::"+str(D[i])+":"+'{:0.2f}'.format(JEspinor1[i][0])
            # Condition: 1. Cut off energy diff between G and D is considered 1 meV
            #            2. D bloch weight needs to be >20% 
            
            # Text positions are shifted by 1e-3 in x-dir and 1e-5 in y-direction
            ax3.text(bwxx[i]+5e-2+changepos[i], HR_Eg[ii][i]+1e-5+changeposy[i],s="("+text+")") 
            #ax3.text(bwxx[i]+1e-3, HR_Eg[ii][i]+1e-5,s="("+text+")")
            
            if (sdone and ((G[i]<L[i] or G[i]<X[i]) or ((JEspinor1[i][0]<1e-3) and D[i]>20))):
                TS_x.append(bwxx[i])
                TS_y.append(HR_Eg[ii][i])
                sdone=0
    if GaAsP_P37_DIT2:
        ax3.axvline(x=-1.4,ymin=0.5,ymax=0.82,color='m',ls='--')
        ax3.text(-1.57,1.69,'DIT',color='m',fontsize=24)
    #--------------------------------------------------------------
    # Join the high resolution Transition points
    #ax.plot(TS_x, TS_y, ls="dotted",lw=2, color='g')
    ax3.plot(TS_x[1:], TS_y[1:], ls="dotted",lw=2, color='m')
    if allbandgaphighresolution:
        ax3.set_xlim(-2.03,-0.93)
        ax3.set_ylim(1.7,1.95)
        ax3.yaxis.set_ticks([1.7,1.8,1.9])
        rect = plt.Rectangle((-2.25,1.62),1.5,0.4,fill=0,lw=2,color='r')
        ax.add_patch(rect)
#-------------------------------------------------------------------
#save plots
fig.tight_layout()
#fig.savefig(savefigfile+'/GaAsP-All2.eps',format='eps',dpi=300)
#fig.savefig(savefigfile+'/GaAsP-P3.7.eps',format='eps',dpi=300)

fig3.tight_layout()
#
#fig3.savefig(savefigfile+'/GaAsP-P3.7-BW3.eps',format='eps',dpi=300)
fig3.savefig(savefigfile+'/GaAsP-P3.7-DIT1.eps',format='eps',dpi=300)
#fig3.savefig(savefigfile+'/GaAsP-highres.eps',format='eps',dpi=300)

#--------------------------------------------------------------
if stage2:
    # DIT vs %N
    plt.figure()
    #plt.title("Direct-Indirect transitin point")
    TS_x = -1*(np.array(TS_x)+0.1)
    TS_x[0] += 0.1
    plt.plot(TS_nc, TS_x, "o-",c='m')
    plt.xlabel("Phosphorus (%)")
    plt.ylabel("Compressive strain (%)")
    #plt.show()
    #--------------------------------------------------------------
    # Fitting DIT vs N%
    xfit=TS_nc; yfit= np.array(TS_x)
    popt = np.polyfit(xfit,yfit,4)
    trendpoly = np.poly1d(popt)
    ncx = np.linspace(TS_nc[0],TS_nc[-1],100)
    #popt, pcov = curve_fit(func1, k, y)
    #plt.plot(ncx, trendpoly(ncx))
    plt.tight_layout()
    #plt.savefig(savefigfile+'/GaAsP-strain-DIT-mean.eps',format='eps',dpi=300)
#--------------------------------------------------------------
if fit:
    # Tensile inset
    left, bottom, width, height = [0.55, 0.55, 0.2, 0.3]
    ax9 = fig8.add_axes([left, bottom, width, height])
    ax9.set_xlabel("Nitrogen (%)")
    ax9.set_ylabel("SMT strain (%)")
    ax9.plot(SMT_x, SMT_y, "o", ls='')

