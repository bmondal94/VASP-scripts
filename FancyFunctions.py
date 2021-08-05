#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 22:15:19 2021

@author: bmondal
"""
from scipy.spatial import ConvexHull
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import mpl_toolkits.mplot3d as a3
from mpl_toolkits.mplot3d import Axes3D

#%%
params = {'legend.fontsize': 18,
          'figure.figsize': (8, 6),
         'axes.labelsize': 24,
         'axes.titlesize': 24,
         'xtick.labelsize':24,
         'ytick.labelsize': 24,
         'errorbar.capsize':2}
plt.rcParams.update(params)
#%%
def Box3d2dBox():
    x = np.array([ 0.16257299, -0.370805  , -1.09232295,  1.62570095,
                  -1.62570095,  1.09232295,  0.370805  , -0.16257299])
    y = np.array([-1.71022499, -0.81153202, -0.52910602, -0.36958599,
                   0.369587  ,  0.52910602,  0.81153202,  1.71022499])
    z = np.array([ 0.22068501, -1.48456001,  1.23566902,  0.469576  ,
                  -0.469576  , -1.23566902,  1.48456001, -0.22068501])
    
    verts = np.c_[x,y,z]
    hull = ConvexHull(verts)
    simplices = hull.simplices
    
    org_triangles = [verts[s] for s in simplices]
    
    class Faces():
        def __init__(self,tri, sig_dig=12, method="convexhull"):
            self.method=method
            self.tri = np.around(np.array(tri), sig_dig)
            self.grpinx = list(range(len(tri)))
            norms = np.around([self.norm(s) for s in self.tri], sig_dig)
            _, self.inv = np.unique(norms,return_inverse=True, axis=0)
    
        def norm(self,sq):
            cr = np.cross(sq[2]-sq[0],sq[1]-sq[0])
            return np.abs(cr/np.linalg.norm(cr))
    
        def isneighbor(self, tr1,tr2):
            a = np.concatenate((tr1,tr2), axis=0)
            return len(a) == len(np.unique(a, axis=0))+2
    
        def order(self, v):
            if len(v) <= 3:
                return v
            v = np.unique(v, axis=0)
            n = self.norm(v[:3])
            y = np.cross(n,v[1]-v[0])
            y = y/np.linalg.norm(y)
            c = np.dot(v, np.c_[v[1]-v[0],y])
            if self.method == "convexhull":
                h = ConvexHull(c)
                return v[h.vertices]
            else:
                mean = np.mean(c,axis=0)
                d = c-mean
                s = np.arctan2(d[:,0], d[:,1])
                return v[np.argsort(s)]
    
        def simplify(self):
            for i, tri1 in enumerate(self.tri):
                for j,tri2 in enumerate(self.tri):
                    if j > i: 
                        if self.isneighbor(tri1,tri2) and \
                           self.inv[i]==self.inv[j]:
                            self.grpinx[j] = self.grpinx[i]
            groups = []
            for i in np.unique(self.grpinx):
                u = self.tri[self.grpinx == i]
                u = np.concatenate([d for d in u])
                u = self.order(u)
                groups.append(u)
            return groups
    
        def order_along_axis(self,faces,axis):
            midpoints = np.array([f.mean(axis=0) for f in faces])
            s = np.dot(np.array(axis),midpoints.T)
            return np.argsort(s)
    
        def remove_last_n(self, faces, order, n=1):
            return np.array(faces)[order][::-1][n:][::-1]
    
    
    
    
    f = Faces(org_triangles, sig_dig=4)
    g = f.simplify()
    order = f.order_along_axis(g, [0,1,0])
    g = f.remove_last_n(g, order, 3)
    
    # Reduce dimension, ommit y axis:
    g2D = g[:,:,[0,2]]
    
    
    fig = plt.figure()
    ax = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122)
    
    
    colors = np.random.rand(len(g),3)
    
    pc = a3.art3d.Poly3DCollection(g,  facecolors=colors, 
                                       edgecolor="k", alpha=0.9)
    ax.add_collection3d(pc)
    
    pc2 = PolyCollection(g2D,  facecolors=colors, 
                                       edgecolor="k", alpha=0.9)
    ax2.add_collection(pc2)
    ax2.autoscale()
    ax2.set_aspect("equal")
    
    
    ax.set_xlim([-1.5,2])
    ax.set_ylim([-1.5,2])
    ax.set_zlim([-1.5,2])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.show()

def rainbowarrow(basex,basey,longx,longy,hl,width=0.1,fraction=1/2,mm_min_lim=0,mm_max_lim=1,cmap=plt.cm.RdYlBu_r,vertical=False,reverse_ar=False,front=False,colorflip=False):
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
    mm_min_lim : float, optional
        Color scale minimum between 0-1. The default is 0.
    mm_max_lim : float, optional
        Color scale maxima between 0-1. The default is 1.
    cmap: color map, optional
        The default is plt.cm.RdYlBu_r
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
    colors=np.linspace(mm_max_lim, mm_min_lim,50)
    
    if reverse_ar:
        basex+=longx
        basey+=longy
        longx = -longx
        longy = -longy
        hl = -hl
        colors=np.flip(colors)
    if colorflip:
        colors=np.flip(colors)
    #plt.arrow(basex,basey+width*0.5,longx,longy,width=width,length_includes_head=True,head_length=hl,fc='k')
    
    
    if vertical:
        colors=np.flip(colors)
        frac=longy*fraction
        rect=plt.Rectangle((basex,basey),width, longy-hl,linewidth = 0,fc = 'k', antialiased=True)
        XY=[[basex+width*2,longy-hl+basey],[basex-width,longy-hl+basey],[basex+width*0.5,basey+longy]]
    else:
        frac=longx*fraction
        rect=plt.Rectangle((basex,basey),longx-hl,width, linewidth = 0,fc = 'k', antialiased=True)
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
    
    trian=plt.Polygon(XY,linewidth = 0,fc = calco, antialiased=True)
    plt.gca().add_patch(trian)
    diffn = x[1]-x[0]
    for i in range(len(x)):
        if vertical:
            rect=plt.Rectangle((basex,x[i]),width, diffn,linewidth = 0,fc = c[i], antialiased=True) 
        else:
            rect=plt.Rectangle((x[i],basey),diffn,width, linewidth = 0,fc = c[i], antialiased=True)
        plt.gca().add_patch(rect)
        
if __name__ == '__main__':
    #Box3d2dBox()
    rainbowarrow
