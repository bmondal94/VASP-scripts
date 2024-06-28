#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 11:05:02 2021

@author: bmondal
"""

import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(precision=5,suppress=True)

params = {'legend.fontsize': 18,
          'figure.figsize': (8, 6),
         'axes.labelsize': 18,
         'axes.titlesize': 18,
         'xtick.labelsize': 18,
         'ytick.labelsize': 18,
         'errorbar.capsize':2}
plt.rcParams.update(params)

#%%

dirname = '/home/bmondal/clusterf/project_DFT_substrate_effect/band_diagram/Ref_level_GaAs/GaAs_110/'

## ------- Read the 1st part of ion information -------------------------------
data_sf = np.loadtxt(dirname+'LOCPOT', skiprows=1, max_rows=1)
data_lv = np.loadtxt(dirname+'LOCPOT', skiprows=2, max_rows=3)*data_sf
data_lv_mag = np.linalg.norm(data_lv, axis=1)

## ------- Read LOCPOT --------------------------------------------------------
data_natom = np.sum(np.loadtxt(dirname+'LOCPOT', skiprows=6, max_rows=1))
imgnor = int(8+data_natom+1) # 8: 1st 8 is heading, data_natom: no. lines with atom coordinates, 1: 1 vacant line at end
data_t = np.loadtxt(dirname+'LOCPOT', skiprows=imgnor, max_rows=1, dtype=int)
data_no_t = np.prod(data_t)
data_no = int(data_no_t//5)
data_locpot = np.genfromtxt(dirname+'LOCPOT', skip_header=imgnor+1, max_rows=data_no)
Extra = int(data_no_t%5)
if Extra:
    data_locpot = np.hstack((data_locpot.flatten(), np.genfromtxt(dirname+'LOCPOT', \
                                skip_header=imgnor+1+len(data_locpot), max_rows=Extra)))

## ---- Reshape LOCPOT: -------------------------------------------------------
## WRITE(IU,FORM) (((C(NX,NY,NZ),NX=1,NGXF),NY=1,NGYF),NZ=1,NGZF)
data_locpot = np.reshape(data_locpot, tuple((data_t[2],data_t[1],data_t[0])))

## --------- dx, dy and dz in mesh --------------------------------------------
dxyz = data_lv_mag/ data_t

## ----------- Averaing along Vacuum ------------------------------------------
ave = input("Which axis do you want to average along (x/X, y/Y, z/Z)?\n")
if ave.upper() == 'Z':
    X_lim, nn = data_lv_mag[2], data_t[2]
    New_data = data_locpot.reshape((data_t[2],data_t[0] * data_t[1])) \
        * dxyz[0] * dxyz[1]
    new_data_ave = np.sum(New_data, axis=1)/(data_lv_mag[0]*data_lv_mag[1])

elif ave.upper() == 'Y':
    X_lim, nn = data_lv_mag[1], data_t[1]
    New_data = data_locpot.transpose((1, 0, 2)).reshape((data_t[2],data_t[1] * data_t[0])) \
        * dxyz[1] * dxyz[0]
    new_data_ave = np.sum(New_data, axis=1)/(data_lv_mag[1]*data_lv_mag[0])

else:
    X_lim, nn = data_lv_mag[0], data_t[0]
    New_data = data_locpot.transpose((0,2,1)).transpose((1, 0, 2)).reshape((data_t[0],data_t[1] * data_t[2])) \
        * dxyz[1] * dxyz[2]
    new_data_ave = np.sum(New_data, axis=1)/(data_lv_mag[1]*data_lv_mag[2])


print(f"Average axis = {ave}")
XX = np.linspace(0, X_lim, int(nn), endpoint=False)    
vl = np.amax(new_data_ave)
print(f'Vacuum level (eV) = {vl:.3f}')
 
#%% --------------- Plotting --------------------------------------------------
fig, ax = plt.subplots()
ax.set_xlabel(r"Distance ($\AA$)")
ax.set_ylabel("Potential (eV)")
ax.plot(XX, new_data_ave, '-')

plt.show()    

#%% ------------- These are data generated using VASPkit ----------------------
# data = np.loadtxt(dirname+'PLANAR_AVERAGE.dat')
# ax.plot(data[:,0], data[:,1], '-')




