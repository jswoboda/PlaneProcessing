#!/usr/bin/env python
"""
Created on Wed Dec 30 16:20:38 2015

@author: John Swoboda
"""

import os, glob
import scipy as sp
import matplotlib
import pdb
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from RadarDataSim.IonoContainer import IonoContainer, makeionocombined
from RadarDataSim.utilFunctions import readconfigfile
from GeoData.GeoData import GeoData
from GeoData.utilityfuncs import readIono
from GeoData.plotting import insertinfo


#%% For sorting
def ke(item):
    if item[0].isdigit():
        return int(item.partition(' ')[0])
    else:
        return float('inf')
#%% Plot input data
def plotinputdata(testdir,imgdir):
    """This will plot all of the input data with each time step as a pcolor images of
    electron density and ion tempreture.
    Inputs
    testdir - The directory with the input data in h5 files formated for the ionocontainer structure.
    imgdir - The directory that holds the images."""

    if os.path.exists(imgdir):
        imgfiles = glob.glob(os.path.join(imgdir,'*.png'))
        for imgf in imgfiles:
            os.remove(imgf)
    else:
        os.mkdir(imgdir)

    filelist = glob.glob(os.path.join(testdir,'Origparams','*.h5'))
    numlist = [os.path.splitext(os.path.split(x)[-1])[0] for x in filelist]
    numdict = {numlist[i]:filelist[i] for i in range(len(filelist))}
    slist = sorted(numlist,key=ke)
    imcount = 0
    filetemplate = 'inputdata'
    dsetname = os.path.split(os.path.dirname(testdir))[-1]
    print "Plotting input data for "+dsetname

    if 'perryplane' in testdir.lower():
        xlim = [-100.,450.]
        xticks = [0.,150.,300.]
        allparams = True
        ncols=3
        figsize = (15,7)
    else:
        xlim = [0.,400.]
        xticks = [0.,150.,300.]
        allparams = False
        ncols=1
        figsize = (5,7)
    ylim = [100.,500.]
    f1 =  True
    for inumn, inum in enumerate(slist):
        print "{0} Input for {1} of {2}".format(dsetname,inumn,len(slist))
        ifile = numdict[inum]
        iono = IonoContainer.readh5(ifile)
        Iono1 = GeoData(readIono,[iono])
        nt = Iono1.times.shape[0]
        rng = sp.sqrt(Iono1.dataloc[:,0]**2+Iono1.dataloc[:,1]**2)*sp.sign(Iono1.dataloc[:,1])
        z = Iono1.dataloc[:,2]
        rngvec = sp.unique(rng)
        zvec = sp.unique(z)
        rngmat = rng.reshape(len(zvec),len(rngvec))
        zmat = z.reshape(len(zvec),len(rngvec))
        Ne = Iono1.data['Ne'].reshape(len(zvec),len(rngvec),nt)
        Ti = Iono1.data['Ti'].reshape(len(zvec),len(rngvec),nt)
        Te = Iono1.data['Te'].reshape(len(zvec),len(rngvec),nt)
        for itimen,itime in enumerate(Iono1.times):
            if f1:
                f1=False
                t0 = itime[0]
            fig ,axmat= plt.subplots(nrows=1,ncols=ncols,facecolor='w',figsize=figsize ,sharey=True)
            if allparams:
                avec = axmat.flatten()
            else:
                avec =[axmat]

            plt.sca(avec[0])
            plt.xticks(xticks)
            plt.tick_params(labelsize=16)
            avec[0].set_xlabel('X Plane in km',fontsize=18)
            avec[0].set_ylabel('Alt in km',fontsize=18)
            pc1 = avec[0].pcolor(rngmat,zmat,Ne[:,:,itimen],cmap = 'plasma',vmin=3e10,vmax=3e11)
            avec[0].set_xlim(xlim)
            avec[0].set_ylim(ylim)
            avec[0].set_title('Electron Density',fontsize=18)
            
           #pc1.set_norm(colors.LogNorm(vmin=5e8,vmax=5e12))
            cb1 = plt.colorbar(pc1, ax=avec[0],format='%.0e')
            cb1.ax.set_xlabel(r'm$^{-3}$')
            if allparams:
                plt.sca(avec[1])
                plt.xticks(xticks)
                plt.tick_params(labelsize=16)
                avec[1].set_xlabel('X Plane in km',fontsize=18)
                pc2 = avec[1].pcolor(rngmat,zmat,Te[:,:,itimen],cmap = 'plasma',vmin=0,vmax=5.e3)
                avec[1].set_xlim(xlim)
                avec[1].set_ylim(ylim)
                avec[1].set_title('Electron Temperature',fontsize=18)

                cb2 = plt.colorbar(pc2, ax=avec[1],format='%.0d')
                cb2.ax.set_xlabel(r'$^{\circ}$K')
                plt.sca(avec[2])
                plt.xticks(xticks)
                plt.tick_params(labelsize=16)
                avec[2].set_xlabel('X Plane in km',fontsize=18)
                pc3 = avec[2].pcolor(rngmat,zmat,Ti[:,:,itimen],cmap = 'plasma',vmin=0,vmax=5.e3)
                avec[2].set_xlim(xlim)
                avec[2].set_ylim(ylim)
                avec[2].set_title('Ion Temperature',fontsize=18)

                cb3 = plt.colorbar(pc3, ax=avec[2],format='%.0d')
                cb3.ax.set_xlabel(r'$^{\circ}$K')

            plt.tight_layout()
#            plt.subplots_adjust(top=0.85)
           # spti = fig.suptitle('Parameters at {0} seconds'.format(int(itime[0]-t0)),fontsize=24)
            fname= '{0:0>3}_'.format(imcount)+filetemplate+'.png'
            plt.savefig(os.path.join(imgdir,fname),dpi=300)
            for ax in avec:
                for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                    label.set_fontsize(20)
            imcount=imcount+1
            plt.close(fig)
#%%Plot output data
def plotoutput(testdir,imgdir,config):
    """This will plot all of the fitted data with each time step as a pcolor images of
    electron density and electron density from power mesurements.
    Inputs
    testdir - The directory with the input data in h5 files formated for the ionocontainer structure.
    imgdir - The directory that holds the images."""
    (sensdict,simparams)=readconfigfile(config)
    tvec = simparams['Timevec']
    if os.path.exists(imgdir):
        imgfiles = glob.glob(os.path.join(imgdir,'*.png'))
        for imgf in imgfiles:
            os.remove(imgf)
    else:
        os.mkdir(imgdir)
    filename = os.path.join(testdir,'Fitted','fitteddata.h5')
    iono = IonoContainer.readh5(filename)
    Iono1 = GeoData(readIono,[iono])
    rngrdr =Iono1.dataloc[:,0]
    sign1 = sp.sign(Iono1.dataloc[:,1])
    el = Iono1.dataloc[:,2]
    elvec,elinv = sp.unique(el,return_inverse=True)
    nbeams = len(elvec)
    nrg = len(rngrdr)/nbeams
    nt = Iono1.times.shape[0]
    Rngrdrmat = sp.reshape(rngrdr,(nrg,nbeams))
    Signmat = sp.reshape(sign1,(nrg,nbeams))
    Elmat = sp.reshape(el,(nrg,nbeams))

    Xmat = Rngrdrmat*Signmat*sp.cos(Elmat*sp.pi/180.)
    Zmat = Rngrdrmat*sp.sin(Elmat*sp.pi/180.)
    Ne = Iono1.data['Ne'].reshape(nrg,nbeams,nt)
    Te = Iono1.data['Te'].reshape(nrg,nbeams,nt)
    Ti = Iono1.data['Ti'].reshape(nrg,nbeams,nt)


    imcount=0
    filetemplate = 'fitteddata'

    dsetname = os.path.split(os.path.dirname(testdir))[-1]
    print "Plotting Output data for "+dsetname

    if 'perryplane' in testdir.lower():
        xlim = [-100.,450.]
        xticks = [0.,150.,300.]
        allparams=True
        ncols=3
        figsize = (15,7)
    else:
        xlim = [0.,400.]
        xticks = [0.,150.,300]
        allparams = False
        ncols=1
        figsize = (5,7)
    ylim = [100.,500]
    for itimen,itime in enumerate(Iono1.times):
        print "{0} Output for {1} of {2}".format(dsetname,itimen,len(Iono1.times))
        
        
        Nemat = Ne[:,:,itimen]
        Timat = Ti[:,:,itimen]
        Temat = Te[:,:,itimen]
        
        fig ,axmat= plt.subplots(nrows=1,ncols=ncols,facecolor='w',figsize=figsize,sharey=True)
        if allparams:
            avec=axmat.flatten()
        else:
            avec=[axmat]

        plt.sca(avec[0])
        avec[0].set_xlabel('X Plane in km',fontsize=18)
        avec[0].set_ylabel('Alt in km',fontsize=18)
        pc1 = avec[0].pcolor(Xmat,Zmat,Nemat,cmap = 'plasma',vmin=3e10,vmax=3e11)
        plt.tick_params(labelsize=16)
        plt.xticks(xticks)
        avec[0].set_xlim(xlim)
        avec[0].set_ylim(ylim)
        avec[0].set_title('Electron Density',fontsize=18)
        cb1 = plt.colorbar(pc1, ax=avec[0],format='%.0e')
        cb1.ax.set_xlabel(r'm$^{-3}$',fontsize=14)
        if allparams:

            
            plt.sca(avec[1])
            plt.tick_params(labelsize=16)
            plt.xticks(xticks)
            avec[1].set_xlabel('X Plane in km',fontsize=18)
            pc2 = avec[1].pcolor(Xmat,Zmat,Temat,cmap = 'plasma',vmin=0,vmax=5.e3)
            avec[1].set_xlim(xlim)
            avec[1].set_ylim(ylim)
            avec[1].set_title('Electron Temperature',fontsize=18)

            cb2 = plt.colorbar(pc2, ax=avec[1],format='%.0d')
            cb2.ax.set_xlabel(r'$^{\circ}$K',fontsize=14)
            plt.sca(avec[2])
            plt.xticks(xticks)
            plt.tick_params(labelsize=16)
            avec[2].set_xlabel('X Plane in km',fontsize=18)
            pc3 = avec[2].pcolor(Xmat,Zmat,Timat,cmap = 'plasma',vmin=0,vmax=5.e3)
            avec[2].set_xlim(xlim)
            avec[2].set_ylim(ylim)
            avec[2].set_title('Ion Temperature',fontsize=18)
            #        for ax in avec:
            #            for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            #                label.set_fontsize(20)

            cb3 = plt.colorbar(pc3, ax=avec[2],format='%.0d')
            cb3.ax.set_xlabel(r'$^{\circ}$K',fontsize=14)	

       # spti = fig.suptitle('Parameters at {0} seconds'.format(int(tvec[itimen])),fontsize=24)

#            ims.append([pc1,pc2])
        plt.tight_layout()

        fname= '{0:0>3}_'.format(imcount)+filetemplate+'.png'
        plt.savefig(os.path.join(imgdir,fname),dpi=300)
        imcount=imcount+1
        plt.close(fig)
def plotlines(inputlist,fitiono,alt,times,paramlist=['Ne','Te','Ti']):

    inputiono = makeionocombined(inputlist)
    Iono1 = GeoData(readIono,[inputiono])
    fitiono = IonoContainer.readh5(fitiono)
    fitGeo = GeoData(readIono,[fitiono])

    paramlist = ['Ne','Te','Ti']
    (x,y,z) = inputiono.Cart_Coords.transpose()
    r = sp.sqrt(x**2+y**2)*sp.sign(y)
    incoords = sp.column_stack((r,z,sp.ones_like(z)))
    ru,zu =[ sp.unique(r),sp.unique(z)]
    Rmat,Zmat = sp.meshgrid(ru,zu)
    zinput = sp.argmin(sp.absolute(zu-alt))
    
    (xf,yf,zf) = fitiono.Cart_Coords.transpose()
    rf = sp.sqrt(xf**2+yf**2)*sp.sign(yf)
    outcoords = sp.column_stack((rf,zf,sp.ones_like(zf)))
    rfu,zfu =[ sp.unique(rf),sp.unique(zf)]
    Rfmat,Zfmat = sp.meshgrid(rfu,zfu)
    zoutput = sp.argmin(sp.absolute(zfu-alt))
    
    fitGeo.interpolate(incoords,Iono1.coordnames,method='linear',fill_value=sp.nan,twodinterp = True,oldcoords=outcoords)

    uz = sp.unique(z)
    ur = sp.unique(r)
    (rmat,zmat) = sp.meshgrid(ur,uz)
    
    inputdata = {}
    for iparm in paramlist:
        if iparm =='Nepow':
            iparm='Ne'
        curdata = Iono1.data[iparm][:,times[0]]
        
        inputdata[iparm] = sp.reshape(curdata,Rmat.shape)[zinput]

    outputdata = {}
    for iparm in paramlist:
        curdata = fitGeo.data[iparm][:, times[1]]
        outputdata[iparm] = sp.reshape(curdata,Rmat.shape)[zinput]
    
    fig, axvec = plt.subplots(len(paramlist),1,sharey=False,figsize=(12,5*len(paramlist)))
    
    for ipn,iparam in enumerate(paramlist):
        ax = axvec[ipn]
        ax2 = ax.twinx()
        p1, = ax.plot(ru,inputdata[iparam],'b-',label='In',linewidth=3)
        p2, = ax.plot(ru,outputdata[iparam],'b--',label='Out',linewidth=3)
        ax.set_title(iparam)
        ax.set_xlabel('X Plane in km')
        gi = sp.gradient(inputdata[iparam])/sp.gradient(ru)
        go = sp.gradient(outputdata[iparam])/sp.gradient(ru)
        p3, = ax2.plot(ru,gi,'g-',label='Grad In',linewidth=3)
        p4, = ax2.plot(ru,go,'g--',label='Grad Out',linewidth=3)
        ax.yaxis.label.set_color(p1.get_color())
        ax2.yaxis.label.set_color(p3.get_color())
        lines = [p1,p2,p3,p4]
        ax.legend(lines,[l.get_label() for l in lines])
    
    plt.tight_layout()    
    return(fig)
