# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 16:58:37 2016

@author: swoboj
"""

import os, glob,inspect,getopt,sys
import shutil
import scipy as sp
import scipy.io as sio
from RadarDataSim.IonoContainer import IonoContainer

def readmattsdata(filename,datadir,outdir,keepspec=[1,2,3,7],angle=20.5):
    d2r=sp.pi/180.
    angr=d2r*angle
    lsp=7
    # Read in Data
    
    inst = sio.loadmat(os.path.join(datadir,filename))
    xg=inst['xg'][0,0]
    x1v = xg['xp']# added to avoid gratting lobes.
    x3v = xg['zp']
    
    [x3mat,x1mat] = sp.meshgrid(x3v,x1v);
    
    x2 = x1v*sp.cos(angr)#x
    x3 = x1v*sp.sin(angr)#y
    x1 = x3v#z
    
    E = x1mat*sp.cos(angr)#x
    N = x1mat*sp.sin(angr)#y
    U = x3mat
    lxs=x3mat.size
    
    Time_Vector = sp.column_stack([inst['t'],inst['t']+15])
    ns = sp.transpose(inst['ns'],[2,1,3])
    print('Loaded densities...');
    
    ns= sp.reshape(ns,[lxs,lsp])
    Ts = sp.transpose(inst['Ts'],[2,1,3])
    print('Loaded temperatures...')
    
    Ts=sp.reshape(Ts,[lxs,lsp])
    vs = sp.transpose(inst['vsx1'],[2,1,3]);
    print('Loaded parallel velocities...\n');
    
    # derive velocity from ExB
    Ex,Ez=sp.gradient(-1*inst['Phi'],xg['x'],xg['z'])
    
    vx1=-1*Ex/xg['Bmag']
    # from looking at the data it seems that the velocity is off by a factor of
    # 10.
    vx1=vx1.flatten()/10.
    vs=sp.reshape(vs,[lxs,lsp])
    vs=sp.sum(ns[:,:lsp-1]*vs[:,:lsp-1],2)/ns[:,lsp]
    v_e= vx1*sp.cos(angr)
    v_n = vx1*sp.sin(angr)
    v_u = vs;
    #change units of velocity to km/s
    Velocity = sp.reshape(sp.concatenate([v_e,v_n,v_u],1),[lxs,1,3]);
    # reduce the number of spcecies
#    if islogical(keepspec)
#        keepspec(end)=true;
#        keepspecnum = sum(keepspec);
#    elseif ~any(keepspec==numspec)
#        keepspec = [keepspec(:),numspec];
#        keepspecnum = length(keepspec);
#    else
#        keepspecnum = length(keepspec);
#    end
    nsout = sp.reshape(ns[:,keepspec],[lxs,1,keepspec])
    Tsout = sp.reshape(Ts[:,keepspec],[lxs,1,keepspec])
    
    
    
    # Put everything in to ionocontainer structure
    Cart_Coords = sp.column_stack([E.flatten(),N.flatten(),U.flatten()])*1e-3
    Sphere_Coords = cart2sphere(Cart_Coords)
    x = x2*1e-3;
    y = x3*1e-3;
    z = x1*1e-3;
    Param_List = sp.concatenate((sp.expand_dims(nsout,nsout.ndim),sp.expand_dims(Tsout,Tsout.ndim)),-1);
    Species = sp.array(['O+','NO+','N2+','O2+','N+', 'H+','e-'])
    Species = Species[keepspec]
    fname=os.path.splitext(filename)[0]
    ionoout=IonoContainer(Cart_Coords,Time_Vector,Param_List,Species,Velocity)
    
    ionoout.saveh5()
        
def cart2sphere(coordlist):
    r2d = 180.0/sp.pi
    d2r = sp.pi/180.0
    
    X_vec = coordlist[:,0]
    Y_vec = coordlist[:,1]
    Z_vec = coordlist[:,2]
    R_vec = sp.sqrt(X_vec**2+Y_vec**2+Z_vec**2)
    Az_vec = sp.arctan2(Y_vec,X_vec)*r2d
    El_vec = sp.arcsin(Z_vec/R_vec)*r2d
    sp_coords = sp.array([R_vec,Az_vec,El_vec]).transpose()
    return sp_coords