#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 14:57:00 2021

@author: bmondal
"""
import sys
import numpy as np      
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib import animation
import glob
import os
sys.path.insert(0,'/home/bmondal/MyFolder/VASP/GaAsP/mcsqs')
import FancyFunctions
from scipy.ndimage.filters import gaussian_filter
import re
from scipy.interpolate import griddata

np.set_printoptions(precision=3,suppress=True)

params = {'legend.fontsize': 18,
          'figure.figsize': (8, 6),
         'axes.labelsize': 24,
         'axes.titlesize': 24,
         'xtick.labelsize':24,
         'ytick.labelsize': 24,
         'errorbar.capsize':2}
plt.rcParams.update(params)

#%%
def ReadBandgap(filename):
    try:
        data = np.genfromtxt(filename+"BiaxialBandgap.dat",dtype=str,delimiter="|")
        NC = np.array(data[:,0],dtype=float)
        ddt = []
        for line in data[:,1:-1]:
            dd =[l.split() for l in line]
            ddt.append(dd)
        type_=1
    except:
        with open(filename+"BiaxialBandgap.dat") as f:
            data = [line.rstrip().split("|") for line in f]
        ddt = [];NC=[]
        for line in data:
            if len(line)>1:
                dd =np.array([l.split() for l in line[1:-1] if not l.startswith('#')], dtype=float)
                ddt.append(dd)
                NC.append(line[0])
        NC = np.array(NC,dtype=float)
        type_=0
    return NC, ddt, type_
#%% Blochweight data
fuldir='/home/bmondal/MyFolder/VASP/GaAsP/mcsqs/Biaxial_new/'
dirnamesavefig='/home/bmondal/MyFolder/VASP/GaAsP/mcsqs/figs'
files = glob.glob(fuldir + 'AllBweights/**/BW_CB_spinor1.txt', recursive=True)
data = {}; S_array=[]
for J,file in enumerate(files):
    data[J] = np.genfromtxt(file, dtype=float)
    S_array.append(file.split("/")[-4:-2])
    
sort_n, sortedSample  = zip(*sorted(enumerate(S_array), key= lambda x:(float(x[1][0][1:]), float(x[1][1][1:]))))

data_sorted={}; K=-1;KK=0;testsample= sortedSample[0][0]
ind2 = []; checkfile={}
# K: index over strain
# KK: index over N conc
for J in range(len(sort_n)):
    #if (testsample in sortedSample[J][0]):
    if re.search(r'^' + testsample + r'$', sortedSample[J][0]):
        K+=1 
    else:
        testsample = sortedSample[J][0]
        ind2.append(K)
        K=0; KK+=1   
    data_sorted[KK,K]=data[sort_n[J]]
    checkfile[KK,K]=sortedSample[J]
ind2.append(K)
# data_sorted: Bloch weight data in sorted form
# sortedSample: Sorted S_array
mytext={}; DIT=[]
for I in range(KK+1): 
    sdone=1;sdone2=0
    dit=[I]
    for J in range(ind2[I]+1):
        #print(I,J)
        BWE=data_sorted[I,J][-1,3] - data_sorted[I,J][0,3] #Delta_m - Gamma
        BW = np.round(data_sorted[I,J][:,-1]*100).astype(int)
        #BW = data_sorted[I,J][:,-1]*100
        bw_L = np.amax(BW[1:216])
        bw_d = np.amax(BW[216:])
        mytext[I,J] =str(BW[0])+":"+str(bw_L)+"::"+str(bw_d)+":"+'{:0.2f}'.format(BWE)
        #mytext[I,J] =str(np.round(BW[0]-bw_L).astype(int))
        if (bw_L<BW[0] and BWE>0): 
            if sdone:
                sdone=0; sdone2=1
                if J!=0: dit.append(J) #exclude J=0 the leftmost point
        else:
            if sdone2:
                sdone2=0;dit.append(J-1)
    DIT.append(dit)        
#%% Band gap data
NC, restdata, datatype = ReadBandgap(fuldir)
if datatype: #1==array, 2==list
    restdata = np.array(restdata, dtype=float)
    strain = restdata[:,:,0]
    bandgap = restdata[:,:,-1]
    eqm_lattice = restdata[:,10,1]
else:
    strain = [np.array(dd[:,0]) for dd in restdata]
    bandgap = [np.array(dd[:,-1]) for dd in restdata]
    eqm_lattice = np.array([dd[np.argwhere(abs(dd[:,0])<0.001)[0],1][0] for dd in restdata])
        
subname = ['GaAs','GaP','Si','Test']
sublattice = np.array([[5.689],[5.475],[5.43],[5.612]]) # Substrate lattice parameters
myco = ['lightsteelblue','red','cornflowerblue','royalblue']
substrain = (sublattice - eqm_lattice)/eqm_lattice*100

BWalternative = 0
drawarrow=1; harraow=1; drawsub=0; heatmap_bound=1;DITfit=1;draw_heatmap=1
#%%
titletext = "Biaxial strain (001)"
xlabeltext = "In-plane strain (%)"
#************************************************
fig, ax = plt.subplots()
#ax.set_title(titletext)
ax.set_xlabel(xlabeltext)
ax.set_ylabel(r"E$_{\mathrm{g}}$")
#ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax.xaxis.set_ticks(strain)
ax.axvline(color='k', ls='--') 

fig1, ax1 = plt.subplots()
#ax1.set_title(titletext)
ax1.set_ylabel(r"Strain(%)")
#ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax1.xaxis.set_ticks(strain)

ax1.axhline(color='k', ls='--')  
# ax1.axvline(x=38.5,ymax=0.55,color='k', ls='--') 
# ax1.text(40,-4,'38.5% P',size=24,c='g')
# ax1.text(15,-1.2,'DIRECT E$_{\mathrm{g}}$',rotation='vertical',size=24,c='b')
# ax1.text(55,-1.6,'INDIRECT E$_{\mathrm{g}}$',rotation=270,size=24,c='b')
ax1.text(10,-1.2,'DIRECT',rotation='vertical',size=24,c='b')
ax1.text(70,-1.5,'INDIRECT',rotation=270,size=24,c='b')
#ax1.axvline(x=30,ymin=0.2,ymax=0.7,color='k', ls='-.',lw=3) 
#ax1.annotate("", xy=(25,-2), xytext=(25, -3.5),arrowprops=dict(arrowstyle="<->"))

if drawsub:
    NC_index = np.argwhere(((NC>35.5) & (NC<39))) # delete extra points in between 35-39
    SubX = np.delete(NC, NC_index)
    for J in range(len(substrain)-1):   
        Suby = np.delete(substrain[J], NC_index)
        ax1.plot(SubX,Suby,'o--',lw=2,color=myco[J],label=subname[J])  
        ax1.legend()

#%%  
left_DIT, right_DIT, fig1dit=[],[],[] #[[0,3.52]]

print("Check the below corresponds are ok. N -> BW-file")
for ii in range(0,12): #len(NC)):
    #ax.plot(strain, bandgap) 
    strain_sorted_index = np.argsort(strain[ii])
    strain_sorted = strain[ii][strain_sorted_index]
    bandgap_sorted = bandgap[ii][strain_sorted_index]
    ax.plot(strain_sorted, bandgap_sorted, marker='o', label = '{:.1f}'.format(NC[ii]))
    for jj in range(len(strain_sorted)):
        print('{:.1f}'.format(NC[ii]),'->',checkfile[ii,jj])
        ax.text(strain_sorted[jj], bandgap_sorted[jj],s="("+mytext[ii,jj]+")")
    if len(DIT[ii])>1:
        for I in range(len(DIT[ii])-1):
            st_tmp = strain_sorted[DIT[ii][I+1]]
            if st_tmp<0:
                left_DIT.append((st_tmp,bandgap_sorted[DIT[ii][I+1]]))
            else:
                right_DIT.append((st_tmp,bandgap_sorted[DIT[ii][I+1]]))
            fig1dit.append([NC[ii],st_tmp])
ax.plot(*zip(*left_DIT),'r--')
ax.plot(*zip(*right_DIT),'r--')
ax.legend(loc=1, ncol=2, labelspacing=0.3, columnspacing=0.8, handlelength=0.7,handletextpad=0.5)

fig1dit = np.array(fig1dit)
if BWalternative:
    fig1ditlater = np.copy(fig1dit)

#%%
strain_cut_off = 8.2 # Delete values greater than cut_off% strain
if heatmap_bound:
    strainarray_cutoff_ind = np.argwhere(abs(fig1dit[:,1])>strain_cut_off)
    fig1dit = np.delete(fig1dit, strainarray_cutoff_ind,axis=0)

#fig1dit_positive = fig1dit[np.where(fig1dit[:,1]>=0)]
#fig1dit_negative = fig1dit[np.where(fig1dit[:,1]<0)]
#ax1.plot(fig1dit_positive[:-2,0],fig1dit_positive[:-2,1],'o-',c='m')
#ax1.plot(fig1dit_negative[:,0],fig1dit_negative[:,1],'o-',c='m')
#ax1.plot([fig1dit_positive[-3,0],fig1dit_negative[-1,0]],[fig1dit_positive[-3,1],fig1dit_negative[-1,1]],'-',c='m')

if DITfit:
    z = np.polyfit(fig1dit[:-2,1], fig1dit[:-2,0],5)
    p = np.poly1d(z)
    xx = np.linspace(-strain_cut_off,5,100)
    _ = ax1.plot(fig1dit[:-2,0], fig1dit[:-2,1],'mo',p(xx),xx, 'm-')

if drawarrow:
        if harraow:
        # Horizontal arrows
            FancyFunctions.rainbowarrow(13,0.32 ,50,0,hl=8,width=0.3,mm_max_lim=0.4,mm_min_lim=0.6,cmap=plt.cm.prism,vertical=False)
            FancyFunctions.rainbowarrow(13,-0.62,50,0,hl=8,width=0.3,mm_max_lim=0.4,mm_min_lim=0.6,cmap=plt.cm.prism,reverse_ar=True,front=True)  
        else:
        #verical arrows
            FancyFunctions.rainbowarrow(24,0,0,4,hl=0.5,width=2,mm_max_lim=0.4,mm_min_lim=0.6,cmap=plt.cm.prism,fraction=0.58,reverse_ar=False,front=False,vertical=True)
            FancyFunctions.rainbowarrow(24,-5,0,5,hl=0.5,width=2,mm_max_lim=0.4,mm_min_lim=0.56,cmap=plt.cm.prism,fraction=0.56,reverse_ar=True,front=False,vertical=True,colorflip=True)

#%%
ax1.set_xticks([0,100])
ax1.set_xticklabels(['As','P'])
ax1.set_ylim(-6,5)
ax1.set_xlim(0,100)

#%% Scatter plot bandgap
# for ii in range(0,len(NC)):
#     xxaxis = [NC[ii]]*len(strain[ii])
#     im = ax1.scatter(xxaxis, strain[ii],cmap=plt.cm.RdYlBu_r, c=bandgap[ii],s=50) 
    
#%% Heatmap plot bandgap
if draw_heatmap:
    Xarray=[]
    for ii in range(0,len(NC)):
        xxaxis = [NC[ii]]*len(strain[ii])
        Xarray+= xxaxis
    
    Xarray = np.array(Xarray)
    strainarray = np.concatenate(strain, axis=0) #strain.flatten()
    Bandgaparray = np.concatenate(bandgap, axis=0) #bandgap.flatten()
    
    if heatmap_bound:
        strain_cut_off = 6
        strainarray_cutoff_ind = np.argwhere(abs(strainarray)>strain_cut_off)
        NC_index = np.argwhere(((Xarray>35.5) & (Xarray<39)))  # Delete heatmap data for N 36-39 data
        strainarray_cutoff_ind = np.concatenate((strainarray_cutoff_ind,NC_index))
        
        strainarray = np.delete(strainarray, strainarray_cutoff_ind)
        Bandgaparray = np.delete(Bandgaparray, strainarray_cutoff_ind)
        Xarray = np.delete(Xarray, strainarray_cutoff_ind)
        
       
        
    
    #%%
    # heatmap, xedges, yedges = np.histogram2d(Xarray, strainarray, weights=Bandgaparray) #, bins=25)
    # extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]]
    # filtered_arr=gaussian_filter(heatmap.T, sigma=2)
    # im1 = ax1.imshow(filtered_arr, extent=extent,cmap=plt.cm.RdYlBu_r, interpolation='bilinear',origin='lower',aspect='auto')
    # im1 = ax1.imshow(heatmap.T, extent=extent,cmap=plt.cm.RdYlBu_r, interpolation='bilinear',origin='lower',aspect='auto')

    extenty_b = np.amin(strainarray)
    extenty_t = np.amax(strainarray)
    extent = [0,100,extenty_b,extenty_t]
    grid_x, grid_y = np.mgrid[0:100:100j, extenty_b:extenty_t:40j]
    points = np.stack((Xarray, strainarray),axis=-1)
    grid_z0 = griddata(points,Bandgaparray, (grid_x, grid_y), method='nearest')
    filtered_arr=gaussian_filter(grid_z0.T, sigma=2)
    im1 = ax1.imshow(filtered_arr, aspect='auto',extent=extent, origin='lower',cmap=plt.cm.RdYlBu_r,interpolation='bilinear')
    #im1 = ax1.imshow(filtered_arr, aspect=8,extent=extent, origin='lower',cmap=plt.cm.RdYlBu_r,interpolation='bilinear')
    #im1 = ax1.imshow(grid_z0.T, extent=extent, origin='lower',cmap=plt.cm.RdYlBu_r, aspect='auto',interpolation='bilinear')
    
    #im1 = ax1.scatter(Xarray, strainarray,cmap=plt.cm.RdYlBu_r, c=Bandgaparray,marker='s',s=0)
    
    #cax = fig1.add_axes([0.15, 1, 0.75, 0.03])
    #fig1.colorbar(im, cax=cax, orientation='horizontal')
    
    
    cbar = fig1.colorbar(im1,format='%.1f')
    cbar.ax.set_ylabel('E$_{\mathrm{g}}$ (eV)', fontsize = 18, weight="bold")
    # cbar = fig.colorbar(im1,format='%.1f', ax=[ax1],location='top',fraction=0.05)
    # cbar.ax.set_ylabel('E$_{\mathrm{g}}$(eV)', fontsize = 18, weight="bold",rotation=0,labelpad=50,va='center_baseline')

#%%
if BWalternative:
    fig1dit_positive = fig1ditlater[np.where(fig1ditlater[:,1]>=0)]
    fig1dit_negative = fig1ditlater[np.where(fig1ditlater[:,1]<0)]
    plt.figure()
    plt.xlim(-6.1,5.1)
    plotdata = []
    for ii in range(0,12): #len(NC)):
        strain_sorted_index = np.argsort(strain[ii])
        strain_sorted = strain[ii][strain_sorted_index]
        
        contrib = []
        for jj in range(len(strain_sorted)):
            contrib.append(int(mytext[ii,jj]))
        plt.plot(strain_sorted, contrib, 'o-',label = '{:.1f}'.format(NC[ii]))
        plotdata.append([fig1dit_positive[ii,-1],contrib[np.argwhere(strain_sorted==fig1dit_positive[ii,-1])[0,0]]]) 
        plotdata.append([fig1dit_negative[ii,-1],contrib[np.argwhere(strain_sorted==fig1dit_negative[ii,-1])[0,0]]])
        #break
    
    plotdata = np.array(plotdata)
    plotdataa = plotdata[np.argsort(plotdata[:,0])]
    plt.xlabel("In-plane strain(%)")
    plt.ylabel("BW$_{\Gamma} -$BW$_{\mathrm{L}}$")
    plt.legend()  
    plt.axvline(x=0,c="k",ls='--')      
    plt.axhline(y=0,c="k",ls='--')
    plt.plot(plotdataa[:,0],plotdataa[:,1],marker='*',color='magenta',ls='-.',lw=2) 
    plt.tight_layout()

#%%


#fig1.tight_layout()
plt.tight_layout()


#plt.savefig(dirnamesavefig+'/GaAsP-Eg-Phasediag-Biaxialstrain-p3.eps', format='eps',bbox_inches = 'tight',dpi=300)

plt.show()

