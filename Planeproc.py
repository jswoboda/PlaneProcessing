#!/usr/bin/env python
"""
This is will create the plane based data and

@author: John Swoboda
"""

import os, glob,inspect,getopt,sys
import shutil
import pdb
import scipy as sp
from RadarDataSim.utilFunctions import makeconfigfile, readconfigfile
from RadarDataSim.IonoContainer import IonoContainer, MakeTestIonoclass
import RadarDataSim.runsim as runsim
from RadarDataSim.analysisplots import analysisdump
import matplotlib.pyplot as plt
from GeoData.GeoData import GeoData
from  GeoData.utilityfuncs import readIono



def makeline(testdir,meanaz,linewidth=1):
    """This will create """
    if not os.path.isdir(testdir):
        os.mkdir(testdir)
    datadir = os.path.join(testdir,'Origparams')
    if not os.path.isdir(datadir):
        os.mkdir(datadir)
#    (sensdict,simparams) = readconfigfile(configfile)
#    azangles = [iang[0] for iang in simparams['angles']]
#    meanaz = sp.mean(azangles)

    rngart = 10
    rngend = 700
    nr = sp.floor((rngend-rngart)/2)
    nz = 350
    nt = 31
    d2r = sp.pi/180.
    rvec = sp.linspace(rngart,rngend,nr+1)
    zvec = sp.linspace(0,2.*nz,nz+1)
    Rmat,Zmat = sp.meshgrid(rvec,zvec)
    Xmat = Rmat*sp.sin(d2r*meanaz)
    Ymat = Rmat*sp.cos(d2r*meanaz)
    coords = sp.column_stack((Xmat.flatten(),Ymat.flatten(),Zmat.flatten()))
    timevec = sp.linspace(0,900,nt)
    rng_vel = -0.5
    xvel = rng_vel*sp.sin(d2r*meanaz)
    yvel = rng_vel*sp.cos(d2r*meanaz)

    parammult = sp.ones_like(Rmat).astype('float64')
    paramadd =sp.zeros_like(Rmat).astype('float64')
    multval = 2.
    start = 450.
   
    # start file
    it=0.
    Icont1 = MakeTestIonoclass(testv=False,testtemp=False,N_0=1e11,z_0=250.0,H_0=50.0,
                               coords=coords,times=sp.array([[it,it+30.]]))
    strtfile = os.path.join(os.path.split(testdir)[0],'startdata.h5')
    if not os.path.isfile(strtfile):
        Icont1.saveh5(strtfile)
    if sp.mod(linewidth,2):
        linewidth+=1
    fwhm = sp.floor(linewidth/2.)
    for itn,it in enumerate(timevec):
        Icont1 = MakeTestIonoclass(testv=False,testtemp=False,N_0=1e11,z_0=250.0,H_0=50.0,
                                   coords=coords,times=sp.array([[it,it+30.]]))
        midloc = sp.argmin(sp.absolute(rvec-(start+rng_vel*it)))

        iloc = sp.arange(midloc-fwhm,midloc+1+fwhm).astype(sp.int64)
        parammult[:,iloc] = multval
        Paramflt = parammult.flatten()
        Icont1.Param_List[:,0,0,0] = Icont1.Param_List[:,0,0,0]*Paramflt#ion density enhancement
        Icont1.Param_List[:,0,1,0] = Icont1.Param_List[:,0,1,0]*Paramflt# electron density enhancement
        Icont1.Param_List[:,0,0,1] = Icont1.Param_List[:,0,0,1]*Paramflt# ion temp enhancement
        paramadd[:,iloc] = xvel
        Icont1.Velocity[:,0,0] = paramadd.flatten()
        paramadd[:,iloc] = yvel
        Icont1.Velocity[:,0,1] = paramadd.flatten()
        paramadd[:,iloc]=0.
        parammult[:,iloc] = 1.
        Icont1.saveh5(os.path.join(datadir,'{0} planeiono.h5'.format(int(it)) ))

def makealldata(basedir,meanaz):
    w_list = sp.arange(1,20,2)
    basestr = 'exp_width_'
    fsuffix = '{0:0'+str(int(sp.ceil(sp.log10(w_list.max()))))+'d}'
    for iwid in w_list:
        dirname = basestr+fsuffix.format(iwid)
        fulldir = os.path.join(basedir,dirname)
        print('Making Data for {0}'.format(fulldir))
        makeline(fulldir,meanaz,linewidth=iwid)
        plotoutdata(dirname,os.path.join(dirname,'Inputimages'))
#%% For sorting
def ke(item):
    if item[0].isdigit():
        return int(item.partition(' ')[0])
    else:
        return float('inf')
#%%
def plotoutdata(testdir,imgdir):


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


    for inum in slist:
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
        Ti = Iono1.data['Ti'].reshape(len(zvec),len(rngvec),nt)

        for itimen,itime in enumerate(Iono1.times):
            fig = plt.figure(facecolor='w',figsize=(14, 8))
            ax1=fig.add_subplot(1,2,1)
            ax2=fig.add_subplot(1,2,2)

            ax1.set_title('Ne')
            ax2.set_title('Ti')
            ax1.set_xlabel('Range in km')
            ax1.set_ylabel('Alt in km')
            ax2.set_xlabel('Range in km')
            ax2.set_ylabel('Alt in km')
            pc1 = ax1.pcolor(rngmat,zmat,Ne[:,:,itimen],cmap = 'jet',vmin=5e10,vmax=2e11)

            pc2 = ax2.pcolor(rngmat,zmat,Ti[:,:,itimen],cmap = 'jet',vmin=3000,vmax=6000)
            ax1.set_xlim([rngmat.min(),rngmat.max()])
            ax2.set_ylim([zmat.min(),zmat.max()])
            spti = fig.suptitle('Parameters at {0} seconds'.format(int(itime[0])))

            cb1 = plt.colorbar(pc1, ax=ax1,format='%.0e')
            cb2 = plt.colorbar(pc2, ax=ax2,format='%.0e')


            fname= '{0:0>3}_'.format(imcount)+filetemplate+'.png'
            plt.savefig(os.path.join(imgdir,fname))
            imcount=imcount+1
            plt.close(fig)

def plotoutput(testdir,imgdir):
    """ Plot fitted data"""
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
    Ti = Iono1.data['Ti'].reshape(len(rngrdrvec),len(elvec),nt)


    imcount=0
    filetemplate = 'fitteddata'

    for itimen,itime in enumerate(Iono1.times):
        fig = plt.figure(facecolor='w',figsize=(14, 8))
        ax1=fig.add_subplot(1,2,1)
        ax2=fig.add_subplot(1,2,2)

        ax1.set_title('Ne')
        ax2.set_title('Ti')
        ax1.set_xlabel('Range in km')
        ax1.set_ylabel('Alt in km')
        ax2.set_xlabel('Range in km')
        ax2.set_ylabel('Alt in km')
        Nemat = Ne[:,:,itimen]
        Timat = Ti[:,:,itimen]
        pc1 = ax1.pcolor(Xmat,Zmat,Nemat,cmap = 'jet',vmin=5e10,vmax=2e11)

        pc2 = ax2.pcolor(Xmat,Zmat,Timat,cmap = 'jet',vmin=1000,vmax=4000)
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

def runradarsims(testpath,funcnamelist=['spectrums','radardata','fitting'],configfile = 'planeproc2.ini',remakealldata=False):
    curpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

    origparamsdir = os.path.join(testpath,'Origparams')
    if not os.path.exists(testpath):
        os.mkdir(testpath)
        print "Making a path for testdata at "+testpath
    if not os.path.exists(origparamsdir):
        os.mkdir(origparamsdir)
        print "Making a path for testdata at "+origparamsdir
        makeline(configfile,testpath)
    # clear everything out
    folderlist =[]
    if 'spectrums' in funcnamelist:
        folderlist.append('Spectrums')
    if 'radardata' in funcnamelist and remakealldata:
        folderlist.append('Radardata')
        folderlist.append('ACF')
    if 'fitting' in funcnamelist:
        folderlist.append('Fitted')
#    folderlist = ['Spectrums','Radardata','ACF','Fitted']
    for ifl in folderlist:
        flist = glob.glob(os.path.join(testpath,ifl,'*.h5'))
        for ifile in flist:
            os.remove(ifile)



    runsim.main(funcnamelist,testpath,configfile,remakealldata)
    try:
        analysisdump(testpath,configfile,'Plane Example')
    except:
        print "Analysis dump failed somewhere"
    plotoutput(testpath,os.path.join(testpath,'fittedimages'))
#    dboxpath =os.path.expanduser(os.path.join('~','Dropbox'))
#    dboxsave=os.path.join(dboxpath,'Planeexample')
#    if os.path.exists(dboxpath):
#        if os.path.exists(dboxsave):
#            shutil.rmtree(dboxsave)
#
#        shutil.copytree(curpath,dboxsave)

if __name__== '__main__':
    argv = sys.argv[1:]

    outstr = 'runsim.py -f <function: origdata, spectrums, radardata, fitting or all> -i <basedir List or all> -c <config> -r <type y to remake data>'
    curpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    try:
        opts, args = getopt.gnu_getopt(argv,"hf:i:c:r:")
    except getopt.GetoptError:
        print(outstr)
        sys.exit(2)

    remakealldata = False
    basedir = os.path.join(curpath,'exp_width_01')
    for opt, arg in opts:
        if opt == '-h':
            print(outstr)
            sys.exit()
        elif opt in ("-i", "--ifile"):
            basedir = arg

        elif opt in ("-c", "--cfile"):
            outdirexist = True
            configfile = arg
        elif opt in ("-f", "--func"):
            funcname = arg

        elif opt in ('-r', "--re"):
            if arg.lower() == 'y':
                remakealldata = True

    funcnamelist= funcname.split()

    if 'origdata' in funcnamelist:
        funcnamelist.remove('origdata')
        makedirs = True
        makealldata(curpath,20.5)

    if 'all' in funcnamelist:

        funcnamelist=['spectrums','radardata','fitting']



    if basedir.lower() == 'all':
        basedirlist = glob.glob('exp_width_*')
    else:
        basedirlist = basedir.split()

    if len(funcnamelist)>0:
        for ibase in basedirlist:
            runradarsims(ibase,funcnamelist,configfile,remakealldata)
