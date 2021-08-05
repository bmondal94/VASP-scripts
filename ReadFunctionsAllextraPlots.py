#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 14:18:28 2021

@author: bmondal
"""
import numpy as np 
def createdata(filename):
    with open(filename+"bandgap.dat", "r") as f:
        data = np.array([np.asarray(l.split("|"), dtype=str) for l in (line.strip() for line in f) if l], dtype=str)
    
    with open(filename+"GetBandgapSupercell.out", "r") as f:
        f.seek(0)
        strainmetadata = np.array([np.asarray(l.split("/")[-2], dtype=str) for l in (line.strip() for line in f) if l.endswith("N00.9") ], dtype=str)
    
    strain_raw = np.core.defchararray.replace(strainmetadata,"S","")
    strain = np.array(np.core.defchararray.replace(strain_raw,"Eqm","0"), dtype=float)
    NC = np.array(data[:,0],dtype=float)
    
    average_bandgap, stderr_bandgap = [], []
    for (ii,i) in enumerate(data):
        bg, sg= [], []
        for (jj,j) in enumerate(i[1:-1]):
            Egg = np.array(j.split(",")[:-1], dtype=np.float)
            if not len(Egg): Egg = np.nan
            Eggmean = np.mean(Egg)
            bg.append(Eggmean)
            sg.append(np.std(Egg))   
        average_bandgap.append(bg)
        stderr_bandgap.append(sg)

    return NC, average_bandgap, stderr_bandgap, strain

def ReadNearestConfigAfterBandgap(filename):
    data = np.genfromtxt(filename+"NearestConfigAfterBandgap.dat",dtype=str,delimiter="|")
    NC = np.array(data[:,0],dtype=float)
    
    ddt = []
    for line in data[:,1:-1]:
        dd =np.array([l.split(" ") for l in line], dtype=float)
        ddt.append(dd)
    restdata = np.array(ddt)
    return NC, restdata

def ReadHighResBandgap(filename):
    data = np.genfromtxt(filename+"HighRes-Bandgap.dat", dtype=str,delimiter="|")
    dd = []
    for line in data:
        dd.append([l.split("/") for l in line[1:]])
    return np.array(dd, dtype=float)

def ReadHighResBandgapHighConc(filename):
    data = np.genfromtxt(filename+"HighRes-Bandgap.dat", dtype=str,delimiter="|")
    dd = []
    for line in data:
        dd.append([l.split("/") for l in line[1:]])
    Nconc = np.asarray(data[:,0],dtype=float)
    return Nconc, np.array(dd, dtype=float)

def ReadOnlyBandgap(filename):
    with open(filename+"OnlyBandgap.dat") as f:
        data = [line.rstrip().split("|") for line in f]
        ddt = [];NC=[];strain=[]
        for line in data:
            if len(line)>1:
                dd =np.array([l.split() for l in line[1:-1] if not l.startswith('#')], dtype=float)
                ddt.append(dd[:,1])
                strain.append(dd[:,0])
                #print(line[0])
                NC+=[line[0]]*len(dd[:,0])
        NC = np.array(NC, dtype=float)
    return NC, np.concatenate(strain), np.concatenate(ddt)