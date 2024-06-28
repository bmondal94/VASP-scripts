#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  2 20:28:54 2021

@author: bmondal
"""
import sys
import numpy as np      
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy.optimize import curve_fit
from brokenaxes import brokenaxes
import matplotlib.transforms
import matplotlib.path
from matplotlib.collections import LineCollection

np.set_printoptions(precision=3,suppress=True)

params = {'legend.fontsize': 20,
          'figure.figsize': (8, 6),
         'axes.labelsize': 28,
         'axes.titlesize': 28,
         'xtick.labelsize': 28,
         'ytick.labelsize': 28,
         'errorbar.capsize':2}
plt.rcParams.update(params)

dirnamesavefig='/home/bmondal/MyFolder/VASP/GaAsP/mcsqs/figs'

#%%
def rainbowarrow(basex,basey,longx,longy,hl,width=0.1,fraction=1/2,vertical=False,reverse_ar=False,front=False):
    '''
    Parameters
    ----------
    basex : float
        x-cordinate of arrow starting point. (left-bottom)
    basey : float
        y-cordinate of arrow starting point. (left-bottom)
    longx : float
        Length of arrow along x including arrow head.
    longy : float
        Length of arrow along y including arrow head.
    hl : float
        Arrow head length.
    width : float, optional
        Width of arrow. The default is 0.1.
    fraction : float, optional
        Fraction of arrow for color. The default is 1/2.
    vertical : bool, optional
        Arrow in vertical direction. The default is False.
    reverse_ar : bool, optional
        Arrow in reverse direction. The default is False.
    front : bool, optional
        Colors in front part of the arrow. The default is False.

    Returns
    -------
    None.

    '''
    colors=np.linspace(1, 0,50)
    cmap=plt.cm.RdYlBu_r
    if reverse_ar:
        basex+=longx
        basey+=longy
        longx = -longx
        longy = -longy
        hl = -hl
        colors=np.flip(colors)
        
    #plt.arrow(basex,basey+width*0.5,longx,longy,width=width,length_includes_head=True,head_length=hl,fc='k')
    
    
    if vertical:
        colors=np.flip(colors)
        frac=longy*fraction
        rect=plt.Rectangle((basex,basey),width, longy-hl,linewidth = 0,fc = 'k')
        XY=[[basex+width*2,longy-hl+basey],[basex-width,longy-hl+basey],[basex+width*0.5,basey+longy]]
    else:
        frac=longx*fraction
        rect=plt.Rectangle((basex,basey),longx-hl,width, linewidth = 0,fc = 'k')
        XY=[[basex+longx-hl,basey+width*2],[basex+longx-hl,basey-width],[basex+longx,basey+width*0.5]]
    plt.gca().add_patch(rect)
    
    c=cmap(colors)
    if front:
        calco=c[-1]
        if vertical:
            x=np.linspace(basey+frac,basey+longy-hl,50) 
        else:
            x=np.linspace(basex+frac,basex+longx-hl,50)
    else:
        calco='k'
        if vertical:
            x=np.linspace(basey, basey+frac,50)
        else:
            x=np.linspace(basex, basex+frac,50)
    
    trian=plt.Polygon(XY,linewidth = 0,fc = calco)
    plt.gca().add_patch(trian)
    diffn = x[1]-x[0]
    for i in range(len(x)):
        if vertical:
            rect=plt.Rectangle((basex,x[i]),width, diffn,linewidth = 0,fc = c[i]) 
        else:
            rect=plt.Rectangle((x[i],basey),diffn,width, linewidth = 0,fc = c[i])
        plt.gca().add_patch(rect)
        

    
#%%
# GaAsP
NC1=[0.0, 0.4629629629629629, 0.9259259259259258, 1.8518518518518516, 2.7777777777777777, \
    3.7037037037037033, 4.62962962962963, 5.555555555555555, 7.4074074074074066, \
        9.25925925925926, 11.11111111111111] 
strain1=[-1.785, -1.5,   -1.5,   -1.5,   -1.4,   -1.4,   -1.3,   -1.3,   -1.2,   -1.1,   -1.0  ]
strain_mean1=[-1.785053963789766, -1.5, -1.5, -1.5, -1.4, -1.4, -1.4, -1.4, -1.3, -1.2, -1.1]

#GaPAs
NC2=[0.0, 0.4629629629629629, 0.9259259259259258, 1.8518518518518516, 2.7777777777777777, \
    3.7037037037037033, 4.62962962962963, 5.555555555555555, 7.4074074074074066, \
        9.25925925925926, 11.11111111111111]
strain2=[2.628440078700301, 2.6, 2.6, 2.5, 2.5, 2.5, 2.4, 2.4, 2.3, 2.2, 2.2]
strain_mean2=[2.628440078700301, 2.6, 2.6, 2.5, 2.5, 2.4, 2.4, 2.3, 2.3, 2.2, 2.1]

#Middle P conc.
## GaAsP side
NC3=[14.814814814814813, 19.90740740740741, 25.0, 30.09259259259259, 35.18518518518518, \
 39.81481481481482, 44.907407407407405, 50.0, 55.092592592592595]
strain3=[-0.9, -0.6, -0.4, -0.3, -0.1, 0.2, 0.3, 0.5, 0.7]
strain_mean3=[-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.2, 0.5, 0.7]

## GaPAs side
NC4=[85.18518518518519, 80.0925925925926, 75.0, 69.9074074074074, 64.81481481481481, 60.18518518518518]
strain4=[2.0, 1.8, 1.6, 1.4, 1.2, 1.0]
strain_mean4=[2.0, 1.8, 1.6, 1.3, 1.1, 0.9]

#%%
NC=np.hstack([np.array(NC1),np.flip(100-np.array(NC2))])
strain=strain1+strain2[-1::-1]
strain_mean=strain_mean1+strain_mean2[-1::-1]

NCall=np.hstack([np.array(NC1+NC3+NC4[-1::-1]),np.flip(100-np.array(NC2))])
strainall=strain1+strain3+strain4[-1::-1]+strain2[-1::-1]
strainall_mean=strain_mean1+strain_mean3+strain_mean4[-1::-1]+strain_mean2[-1::-1]

#%%
allpoints=0; add_rect=0; phasediag=1;drawarrow=0;drawarrow2=0;harrow=0;fit=0

plt.figure()
#plt.xlabel("Phosphorus (%)")
plt.xticks([0,100],['As','P'])
plt.ylabel('Strain (%)')
plt.xlim(-0,100)
plt.ylim(-2.01,2.85)

if phasediag:
    plt.plot(NCall,strainall,'o-',c='m') #,label='Max')
    #plt.plot(NCall,strainall_mean,'*-',c='r',label='Average')
    if fit:
        # Linear fit 
        A = np.vstack([NCall, np.ones(len(NCall))]).T
        m, c = np.linalg.lstsq(A, strainall, rcond=None)[0]
        n_xx = np.linspace(0,100,num=101)
        yyy = m*n_xx+ c
        plt.plot(n_xx,yyy,label='Linear fit: m={:.2f},c={:.2f}'.format(m,c))
        m, c = np.linalg.lstsq(A, strainall_mean, rcond=None)[0]
        yyy = m*n_xx+ c
        #plt.plot(n_xx,yyy,c='r',label='Linear fit: m={:.2f},c={:.2f}'.format(m,c))
        plt.legend(handlelength=1)
    # plt.text(25,1,'DIRECT E$_{\mathrm{g}}$',size=24,c='b',rotation=35)
    # plt.text(55,-1,'INDIRECT E$_{\mathrm{g}}$',size=24,c='b',rotation=35)
    plt.text(25,1,'DIRECT',size=24,c='b',rotation=35)
    plt.text(55,-1,'INDIRECT',size=24,c='b',rotation=35)
    
    plt.axhline(c="k",ls='--')
    #plt.axvline(x=37.0,ymax=0.42,color='g',ls='--')
    # plt.axvline(x=37.5,ymax=0.42,color='r',ls='--')
    #plt.text(13,-1.5,'37.0% P',size=24,c='g')
    # plt.text(39,-1.5,'37.5% P',size=24,c='r')
    if drawarrow:
        if harrow:
            # Horizontal arrows
            rainbowarrow(25,0.75,50,0,hl=8,fraction=0.62,width=0.1,vertical=False)
            rainbowarrow(25,0.45,50,0,hl=8,fraction=0.49,reverse_ar=True,front=True) 
        else:
            #verical arrows
            rainbowarrow(46,-0.8,0,2.5,hl=0.5,width=2,fraction=0.47,reverse_ar=False,front=True,vertical=True)
            rainbowarrow(54,-0.8,0,2.5,hl=0.5,width=2,fraction=0.41,reverse_ar=True,front=False,vertical=True)
else:
    plt.plot(NC,strain,'o-',c='m')
    #plt.plot(NC,strain_mean,'o-',c='m')
    plt.axhline(c="k",ls='--')
    if add_rect:
        cir=plt.Rectangle((35,-0.2), 6,0.4,color='r',lw='2',fill=0)
        plt.gca().add_patch(cir)
        plt.axvline(x=38,ymax=0.42,color='b',ls='--')
        plt.text(38.5,-1,'38.0% P',size=24,c='b')
    if allpoints:
        plt.plot(NC3,strain3,'*',c='purple')
        plt.plot(NC4,strain4,'*',c='purple')

#plt.legend(handlelength=1)
plt.tight_layout()
#plt.savefig(dirnamesavefig+'/strain-P3.eps', format='eps',dpi=300)
#plt.savefig(dirnamesavefig+'/GaAsP-Eg-Phasediag3.eps', format='eps',dpi=300)
#plt.savefig(dirnamesavefig+'/GaAsP-StrainEg-fit.eps', format='eps',dpi=300)
#plt.savefig(dirnamesavefig+'/GaAsP-average-mean2.eps', format='eps',dpi=300)

'''
#%%
fig = plt.figure(tight_layout=False, figsize=(16, 8.5))
bax = brokenaxes(xlims=((-0.1, 12), (88.5, 100.1)), ylims=((-1.85, -1), (2, 2.7)), \
                 hspace=.2,wspace=0.1,d=0.01,despine=False)
bax.plot(NC1,strain1,'o-',c='m')
bax.plot(np.flip(100-np.array(NC2)),strain2[-1::-1],'o-',c='m')
bax.set_ylabel('Strain (%)',labelpad=100)
bax.set_xlabel("Phosphorus (%)",labelpad=40)

plt.savefig(dirnamesavefig+'/strain-P.eps', format='eps',dpi=300)
#plt.tight_layout()
'''
#%%
draw_heatmap =1;heatmap_bound=1;HighResolutiondata=1;OnlyBandgapRegion=1
sys.path.insert(0,'/home/bmondal/MyFolder/VASP/GaAsP/mcsqs')
import FancyFunctions
if draw_heatmap:
    from scipy.interpolate import griddata
    from scipy.ndimage.filters import gaussian_filter
    sys.path.insert(0,'/home/bmondal/MyFolder/VASP/GaAsP/mcsqs')
    import ReadFunctionsAllextraPlots
    #%%
    filename = "/home/bmondal/MyFolder/VASP/GaAsP/mcsqs/data/ThermExComp"
    NC, average_bandgap1, stderr_bandgap1, strain = ReadFunctionsAllextraPlots.createdata(filename)
    HR_NC1 = np.copy(NC)
    NC1 = np.repeat(NC,len(strain))
    strain1 = np.tile(strain,len(NC))
    average_bandgap1 = np.concatenate(np.array(average_bandgap1))
    filename = "/home/bmondal/MyFolder/VASP/GaPAs/mcsqs/data/ThermExComp"
    NC, average_bandgap4, stderr_bandgap4, strain = ReadFunctionsAllextraPlots.createdata(filename)
    HR_NC4 = np.copy(NC)
    NC4 = np.repeat(NC,len(strain))
    strain4 = np.tile(strain,len(NC))
    average_bandgap4 = np.concatenate(np.array(average_bandgap4))
    
    filename ="/home/bmondal/MyFolder/VASP/GaAsP/mcsqs/data/HighConc-ThermExComp"
    NC, restdata = ReadFunctionsAllextraPlots.ReadNearestConfigAfterBandgap(filename)
    HR_NC2 = np.copy(NC)
    NconfigBG2 = restdata[:,:,-1]
    latticep=restdata[:,:,0]
    eqml = latticep[:,0].reshape((len(latticep),1))
    strain = (latticep - eqml)/eqml*100.
    strain[-3]=strain[-2]=strain[-1]=np.array([-1,-0.5,0.0,0.5,1,1.5])
    NC2 = np.repeat(NC,np.shape(strain)[1])
    strain2 = np.concatenate(strain)
    average_bandgap2 = np.concatenate(NconfigBG2)
    
    filename ="/home/bmondal/MyFolder/VASP/GaPAs/mcsqs/data/HighConc-ThermExComp"
    NC, restdata = ReadFunctionsAllextraPlots.ReadNearestConfigAfterBandgap(filename)
    HR_NC3 = np.copy(NC)
    NconfigBG3 = restdata[:,:,-1]
    latticep=restdata[:,:,0]
    eqml = latticep[:,0].reshape((len(latticep),1))
    strain = (latticep - eqml)/eqml*100.
    NC3 = np.repeat(NC,np.shape(strain)[1])
    strain3 = np.concatenate(strain)
    average_bandgap3 = np.concatenate(NconfigBG3)
    
    
    
    Xarray = np.concatenate((NC1,NC2,NC3,NC4))
    strainarray =np.concatenate((strain1,strain2,strain3,strain4))
    Bandgaparray = np.concatenate((average_bandgap1,average_bandgap2,average_bandgap3,average_bandgap4))
    
    if HighResolutiondata:
        #High Resolution part
        filename = "/home/bmondal/MyFolder/VASP/GaAsP/mcsqs/data/ThermExComp"
        HR_all = ReadFunctionsAllextraPlots.ReadHighResBandgap(filename)
        HR_Eg = HR_all[:,:,1]
        HR_strain = HR_all[:,:,0]
        HR_NC1 = np.repeat(HR_NC1, np.shape(HR_strain)[1])
        HR_Eg1 = np.concatenate(HR_Eg)
        HR_strain1 = np.concatenate(HR_strain)
        
        filename = "/home/bmondal/MyFolder/VASP/GaPAs/mcsqs/data/ThermExComp"
        HR_all = ReadFunctionsAllextraPlots.ReadHighResBandgap(filename)
        HR_Eg = HR_all[:,:,1]
        HR_strain = HR_all[:,:,0]
        HR_NC4 = np.repeat(HR_NC4, np.shape(HR_strain)[1])
        HR_Eg4 = np.concatenate(HR_Eg)
        HR_strain4 = np.concatenate(HR_strain)
        
        filename ="/home/bmondal/MyFolder/VASP/GaAsP/mcsqs/data/HighConc-ThermExComp"
        Nconc, HR_all = ReadFunctionsAllextraPlots.ReadHighResBandgapHighConc(filename)
        HR_Eg = HR_all[:,:,1]
        HR_strain = HR_all[:,:,0]
        HR_NC2 = np.repeat(HR_NC2, np.shape(HR_strain)[1])
        HR_Eg2 = np.concatenate(HR_Eg)
        HR_strain2 = np.concatenate(HR_strain)
        
        filename ="/home/bmondal/MyFolder/VASP/GaPAs/mcsqs/data/HighConc-ThermExComp"
        Nconc, HR_all = ReadFunctionsAllextraPlots.ReadHighResBandgapHighConc(filename)
        HR_Eg = HR_all[:,:,1]
        HR_strain = HR_all[:,:,0]
        HR_NC3 = np.repeat(HR_NC3, np.shape(HR_strain)[1])
        HR_Eg3 = np.concatenate(HR_Eg)
        HR_strain3 = np.concatenate(HR_strain)
        
        Xarray = np.hstack((Xarray,HR_NC1,HR_NC2,HR_NC3,HR_NC4))
        strainarray =np.hstack((strainarray,HR_strain1,HR_strain2,HR_strain3,HR_strain4))
        Bandgaparray = np.hstack((Bandgaparray,HR_Eg1,HR_Eg2,HR_Eg3,HR_Eg4))
        
    if OnlyBandgapRegion:
        filename = "/home/bmondal/MyFolder/VASP/GaAsP/mcsqs/data/ThermExComp"
        NC, strain, Bandg = ReadFunctionsAllextraPlots.ReadOnlyBandgap(filename)
        
        Xarray = np.hstack((Xarray,NC))
        strainarray =np.hstack((strainarray,strain))
        Bandgaparray = np.hstack((Bandgaparray,Bandg))
        
    
    #%%
    if heatmap_bound:
        strain_cut_off = 3
        strainarray_cutoff_ind = np.argwhere(abs(strainarray)>strain_cut_off)
        strainarray = np.delete(strainarray, strainarray_cutoff_ind)
        Bandgaparray = np.delete(Bandgaparray, strainarray_cutoff_ind)
        Xarray = np.delete(Xarray, strainarray_cutoff_ind)

    extenty_b = np.amin(strainarray)
    extenty_t = np.amax(strainarray)
    extent = [0,100,extenty_b,extenty_t]
    grid_x, grid_y = np.mgrid[0:100:100j, extenty_b:extenty_t:40j]
    points = np.stack((Xarray, strainarray),axis=-1)
    grid_z0 = griddata(points,Bandgaparray, (grid_x, grid_y), method='nearest')
    filtered_arr=gaussian_filter(grid_z0.T, sigma=2)
    im1 = plt.imshow(filtered_arr, aspect='auto',extent=extent, origin='lower',cmap=plt.cm.RdYlBu_r,interpolation='bilinear')
    #im1 = ax1.imshow(grid_z0.T, extent=extent, origin='lower',cmap=plt.cm.RdYlBu_r, aspect='auto',interpolation='bilinear')    
    #im1 = plt.scatter(Xarray, strainarray,cmap=plt.cm.RdYlBu_r, c=Bandgaparray,marker='o',s=40)   
    cbar = plt.colorbar(im1,format='%.1f')
    cbar.ax.set_ylabel('E$_{\mathrm{g}}$ (eV)', fontsize = 18, weight="bold")
    #cbar = plt.colorbar(im1,format='%.1f', ax=[ax1],location='top',fraction=0.05)
    #cbar.ax.set_ylabel('E$_{\mathrm{g}}$(eV)', fontsize = 18, weight="bold",rotation=0,labelpad=40,va='center_baseline')
if drawarrow2:
    # Horizontal arrows
    FancyFunctions.rainbowarrow(13,0.1,50,0,hl=8,width=0.1,fraction=0.51,mm_max_lim=0.4,mm_min_lim=0.6,cmap=plt.cm.prism,vertical=False)
    FancyFunctions.rainbowarrow(13,-0.2,50,0,hl=8,width=0.1,fraction=0.58,mm_max_lim=0.4,mm_min_lim=0.6,cmap=plt.cm.prism,reverse_ar=True,front=True)  
    #verical arrows
    FancyFunctions.rainbowarrow(33,-1.3,0,2.5,hl=0.5,width=2,mm_max_lim=0.4,mm_min_lim=0.6,cmap=plt.cm.prism,fraction=0.46,reverse_ar=False,front=True,vertical=True)
    FancyFunctions.rainbowarrow(38,-1.3,0,2.5,hl=0.5,width=2,mm_max_lim=0.4,mm_min_lim=0.56,cmap=plt.cm.prism,fraction=0.41,reverse_ar=True,front=False,vertical=True,colorflip=True)
#%%
plt.tight_layout()
#plt.savefig(dirnamesavefig+'/GaAsP-Eg-Phasediag-ThermExComp.eps', format='eps',dpi=300)
#plt.savefig(dirnamesavefig+'/GaAsP-Eg-Phasediag4.eps', format='eps',dpi=300)
