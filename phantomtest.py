#!/usr/bin/env python
"""
Created on Wed Oct  5 13:27:40 2016

@author: John Swoboda
"""

import scipy as sp
import os,glob
from RadarDataSim.IonoContainer import IonoContainer, makeionocombined
from RadarDataSim.operators import RadarSpaceTimeOperator
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from PlaneProcPlot import plotacf
from PlaneProcUtil import makesimpledata
import pdb
def phantest(inputfile,outdir,configfile,timevec):
    
    
    ne_red=1e-10
    nemin,nemax=[0.*ne_red,3e11*ne_red]
    xlim = [-200.,360.]
    xticks = [-150.,0.,150.,300.]
    defmap = 'viridis'# color map
    fs=18# fontsize
    fscb=12
    figsize = (10,7)
    ylim = [100.,500.]
    
#    iono1=makesimpledata(inputfile,timevec= None,begx=0.,begz=300.,vx=500.)
    iono1=makeionocombined(os.path.split(inputfile)[0])
    iono2=iono1.deepcopy()
    iono2.Param_List[:,1:]=0.
#    sum1=sp.sum(iono1.Param_List,axis=1)
#    iono1.Param_List[:,:,:]=0.
#    iono1.Param_List[:,0]=sum1
    RSTO1 = RadarSpaceTimeOperator(iono1,configfile,timevec,mattype='matrix')
    RSTO2 = RadarSpaceTimeOperator(iono2,configfile,timevec,mattype='sim')
    iono1new=RSTO1.mult_iono(iono1)
    iono2new=RSTO2.mult_iono(iono2)
    
    rngrdr =iono1new.Sphere_Coords[:,0].astype('float32')
    sign1 = sp.sign(iono1new.Sphere_Coords[:,1])
    el = iono1new.Sphere_Coords[:,2].astype('float32')
    elvec,elinv = sp.unique(el,return_inverse=True)
    nbeams = len(elvec)
    nrg = len(rngrdr)/nbeams
    
    Rngrdrmat = sp.reshape(rngrdr,(nrg,nbeams))
    Signmat = sp.reshape(sign1,(nrg,nbeams))
    Elmat = sp.reshape(el,(nrg,nbeams))

    Xmat = Rngrdrmat*Signmat*sp.cos(Elmat*sp.pi/180.)
    Zmat = Rngrdrmat*sp.sin(Elmat*sp.pi/180.)
    Ne1 = iono1new.Param_List[:,:,0].reshape(nrg,nbeams,1).real*ne_red
    Nemat1 = Ne1[:,:,0]
    Ne2 = iono2new.Param_List[:,:,0].reshape(nrg,nbeams,1).real*ne_red
    Nemat2 = Ne2[:,:,0]    
    fig ,axmat= plt.subplots(nrows=1,ncols=2,facecolor='w',figsize=figsize,sharey=True)

    avec=axmat.flatten()

    plt.sca(avec[0])
    avec[0].set_xlabel('X Plane in km',fontsize=fs)
    avec[0].set_ylabel('Alt in km',fontsize=fs)
    pc1 = avec[0].pcolor(Xmat,Zmat,Nemat1,cmap = defmap,vmin=nemin,vmax=nemax)
    plt.tick_params(labelsize=16)
    plt.xticks(xticks)
    avec[0].set_xlim(xlim)
    avec[0].set_ylim(ylim)
    avec[0].set_title('Operator No motion',fontsize=fs)
    tick_locator = ticker.MaxNLocator(nbins=5)
        
    cb1 = plt.colorbar(pc1, ax=avec[0])
    cb1.ax.xaxis.set_label_position('top')
    cb1.ax.set_xlabel('')
    cb1.locator = tick_locator
    cb1.update_ticks()
    
    plt.sca(avec[1])
    avec[1].set_xlabel('X Plane in km',fontsize=fs)
    avec[1].set_ylabel('Alt in km',fontsize=fs)
    pc1 = avec[1].pcolor(Xmat,Zmat,Nemat2,cmap = defmap,vmin=nemin,vmax=nemax)
    plt.tick_params(labelsize=fs)
    plt.xticks(xticks)
    avec[1].set_xlim(xlim)
    avec[1].set_ylim(ylim)
    avec[1].set_title('Operator for Motion',fontsize=fs)
    tick_locator = ticker.MaxNLocator(nbins=5)
        
    cb2 = plt.colorbar(pc1, ax=avec[1])
    cb2.ax.xaxis.set_label_position('top')
    cb2.ax.set_xlabel('')
    cb2.locator = tick_locator
    cb2.update_ticks()
    plt.tight_layout()
    return fig,avec,Nemat1,Nemat2
    
if __name__== '__main__':
    f1='~/DATA/PerryPlane/MattsDataInv/Origparams/18000.h5'
    f1=os.path.expanduser(f1)
    outdir=os.path.expanduser('~/DATA/PerryPlane/SimpleData')
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    timevec= sp.column_stack((sp.linspace(0,60.,4,False),sp.linspace(0,60.,4,False)+15.))
    configfile='~/DATA/PerryPlane/planeproc2red.ini'
    configfile=os.path.expanduser(configfile)
    fig,axvec=phantest(f1,outdir,configfile,timevec)