#!/usr/bin/env python
"""
Created on Wed Dec 30 16:20:38 2015

@author: John Swoboda
"""

import os, glob
import scipy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from RadarDataSim.IonoContainer import IonoContainer
from GeoData.GeoData import GeoData
from GeoData.utilityfuncs import readIono


#%% For sorting
def ke(item):
    if item[0].isdigit():
        return int(item.partition(' ')[0])
    else:
        return float('inf')
#%% Plot input data
def plotoutdata(testdir,imgdir):
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
    filetemplate = 'outputdata'
    dsetname = os.path.split(os.path.dirname(testdir))[-1]
    print "Plotting input data for "+dsetname

    for inumn, inum in enumerate(slist):
        print "{0} Input for {1} of {2}".format(dsetname,inumn,len(slist))
        ifile = numdict[inum]
        iono = IonoContainer.readh5(ifile)
        Iono1 = GeoData(readIono,[iono])
        nt = Iono1.times.shape[0]
        rng = sp.sqrt(Iono1.dataloc[:,0]**2+Iono1.dataloc[:,1]**2)
        z = Iono1.dataloc[:,2]
        rngvec = sp.unique(rng)
        zvec = sp.unique(z)
        rngmat = rng.reshape(len(zvec),len(rngvec))
        zmat = z.reshape(len(zvec),len(rngvec))
        Ne = Iono1.data['Ne'].reshape(len(zvec),len(rngvec),nt)
#        Ti = Iono1.data['Ti'].reshape(len(zvec),len(rngvec),nt)

        for itimen,itime in enumerate(Iono1.times):
            fig = plt.figure(facecolor='w',figsize=(14, 8))
            ax1=fig.add_subplot(1,1,1)


            ax1.set_title('Ne')

            ax1.set_xlabel('Range in km')
            ax1.set_ylabel('Alt in km')

            pc1 = ax1.pcolor(rngmat,zmat,Ne[:,:,itimen],cmap = 'jet',vmin=5e10,vmax=2e11)
            ax1.set_xlim([rngmat.min(),rngmat.max()])
            spti = fig.suptitle('Parameters at {0} seconds'.format(int(itime[0])))

            cb1 = plt.colorbar(pc1, ax=ax1,format='%.0e')
#            ax2=fig.add_subplot(1,2,2)
#            ax2.set_title('Ti')
#            ax2.set_xlabel('Range in km')
#            ax2.set_ylabel('Alt in km')
#            pc2 = ax2.pcolor(rngmat,zmat,Ti[:,:,itimen],cmap = 'jet',vmin=3000,vmax=6000)

#            ax2.set_ylim([zmat.min(),zmat.max()])

#            cb2 = plt.colorbar(pc2, ax=ax2,format='%.0e')


            fname= '{0:0>3}_'.format(imcount)+filetemplate+'.png'
            plt.savefig(os.path.join(imgdir,fname))
            imcount=imcount+1
            plt.close(fig)
#%%Plot output data
def plotoutput(testdir,imgdir):
    """This will plot all of the fitted data with each time step as a pcolor images of
    electron density and electron density from power mesurements.
    Inputs
    testdir - The directory with the input data in h5 files formated for the ionocontainer structure.
    imgdir - The directory that holds the images."""
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
    el = Iono1.dataloc[:,2]
    rngrdrvec,rndinv = sp.unique(rngrdr,return_inverse=True)
    elvec,elinv = sp.unique(el,return_inverse=True)
    nt = Iono1.times.shape[0]
    Rngrdrmat = sp.reshape(rngrdr,(len(rngrdrvec),len(elvec)))
    Elmat = sp.reshape(el,(len(rngrdrvec),len(elvec)))


    Xmat = Rngrdrmat*sp.cos(Elmat*sp.pi/180.)
    Zmat = Rngrdrmat*sp.sin(Elmat*sp.pi/180.)
    Ne = Iono1.data['Ne'].reshape(len(rngrdrvec),len(elvec),nt)
    Ti = Iono1.data['Nepow'].reshape(len(rngrdrvec),len(elvec),nt)


    imcount=0
    filetemplate = 'fitteddata'

    dsetname = os.path.split(os.path.dirname(testdir))[-1]
    print "Plotting Output data for "+dsetname
    for itimen,itime in enumerate(Iono1.times):
        print "{0} Output for {1} of {2}".format(dsetname,itimen,len(Iono1.times))
        fig = plt.figure(facecolor='w',figsize=(14, 8))
        ax1=fig.add_subplot(1,2,1)
        ax2=fig.add_subplot(1,2,2)

        ax1.set_title('Ne')
        ax2.set_title('Ne From Power')
        ax1.set_xlabel('Range in km')
        ax1.set_ylabel('Alt in km')
        ax2.set_xlabel('Range in km')
        ax2.set_ylabel('Alt in km')
        Nemat = Ne[:,:,itimen]
        Timat = Ti[:,:,itimen]
        pc1 = ax1.pcolor(Xmat,Zmat,Nemat,cmap = 'jet',vmin=5e10,vmax=3e11)

        pc2 = ax2.pcolor(Xmat,Zmat,Timat,cmap = 'jet',vmin=5e10,vmax=3e11)
        ax1.set_xlim([Xmat.min(),Xmat.max()])
        ax2.set_ylim([Zmat.min(),Zmat.max()])
        spti = fig.suptitle('Parameters at {0} seconds'.format(int(itime[0])))
#            if imcount==0:
        cb1 = plt.colorbar(pc1, ax=ax1,format='%.0e')
        cb2 = plt.colorbar(pc2, ax=ax2,format='%.0e')
#            ims.append([pc1,pc2])

        fname= '{0:0>3}_'.format(imcount)+filetemplate+'.png'
        plt.savefig(os.path.join(imgdir,fname))
        imcount=imcount+1
        plt.close(fig)
