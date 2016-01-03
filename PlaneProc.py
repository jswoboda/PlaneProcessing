#!/usr/bin/env python
"""
This is will create the plane based data sets for the ISR errors paper

@author: John Swoboda
"""

import os, glob,inspect,getopt,sys
import shutil
import pdb
import scipy as sp
import matplotlib
matplotlib.use('Agg')
from RadarDataSim.IonoContainer import MakeTestIonoclass
import RadarDataSim.runsim as runsim
from RadarDataSim.analysisplots import analysisdump
from RadarDataSim.utilFunctions import readconfigfile
from PlaneProcPlot import plotoutdata,plotoutput


def makeline(testdir,meanaz,linewidth=1,multval = 5.,start = 450.,rng_vel = -0.5):
    """This will create one of the plane datasets given the desired azmuth location
    and line width in samples. The line will be an enhancement in electron density
    Inputs:
    testdir - directory that will house the simulation data.
    meanaz - the Azimuth angle the data will lie on in degrees.
    linewidth - The width of the enhancement in samples
    multval - The enhancement's multiplicaive value compared to the background.
    start - The start location of the enhancement. In km from radar.
    rng_vel - The velocity of the enhancement - is toward the radar in km/s
    """
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

    xvel = rng_vel*sp.sin(d2r*meanaz)
    yvel = rng_vel*sp.cos(d2r*meanaz)

    parammult = sp.ones_like(Rmat).astype('float64')
    paramadd =sp.zeros_like(Rmat).astype('float64')


    # start file
    it=0.
    Icont1 = MakeTestIonoclass(testv=False,testtemp=False,N_0=1e11,z_0=250.0,H_0=50.0,
                               coords=coords,times=sp.array([[it,it+30.]]))
    strtstr = 'startdata{0}.h5'.format(str(meanaz).replace('.','_',1))
    strtfile = os.path.join(os.path.split(testdir)[0],strtstr)
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
        paramadd[:,iloc] = xvel
        Icont1.Velocity[:,0,0] = paramadd.flatten()
        paramadd[:,iloc] = yvel
        Icont1.Velocity[:,0,1] = paramadd.flatten()
        paramadd[:,iloc]=0.
        parammult[:,iloc] = 1.
        Icont1.saveh5(os.path.join(datadir,'{0} planeiono.h5'.format(int(it)) ))

def makealldata(basedir,meanaz,multval = 5.):
    """This will make data sets of different enhancement widths and one stationary enhancement.
    Inputs:
    basedir - directory that will house all the different simulation data.
    meanaz - the Azimuth angle the data will lie on in degrees.

    multval - The enhancement's multiplicaive value compared to the background.
    """

    w_list = sp.arange(1,20,2)
    basestr = 'exp_width_'
    fsuffix = '{0:0'+str(int(sp.ceil(sp.log10(w_list.max()))))+'d}'
    for iwid in w_list:
        dirname = basestr+fsuffix.format(iwid)
        fulldir = os.path.join(basedir,dirname)
        origparamsdir = os.path.join(fulldir,'Origparams')
        if not os.path.exists(fulldir):
            os.mkdir(fulldir)
            print "Making a path for testdata at "+dirname
        if not os.path.exists(origparamsdir):
            os.mkdir(origparamsdir)
            print "Making a path for testdata at "+origparamsdir
        print('Making Data for {0}'.format(fulldir))
        makeline(fulldir,meanaz,linewidth=iwid,multval=multval,)
        plotoutdata(dirname,os.path.join(fulldir,'Inputimages'))

    # for stationary data

    iwid = 1
    dirname = basestr+'stat_'+fsuffix.format(iwid)

    fulldir = os.path.join(basedir,dirname)
    origparamsdir = os.path.join(fulldir,'Origparams')
    if not os.path.exists(fulldir):
        os.mkdir(fulldir)
        print "Making a path for testdata at "+dirname
    if not os.path.exists(origparamsdir):
        os.mkdir(origparamsdir)
        print "Making a path for testdata at "+origparamsdir
    print('Making Data for {0}'.format(fulldir))
    makeline(fulldir,meanaz,linewidth=1,multval=multval,start = 150.,rng_vel = -0.0)

#%% For sorting
def ke(item):
    if item[0].isdigit():
        return int(item.partition(' ')[0])
    else:
        return float('inf')
#%%

def runradarsims(testpath,funcnamelist=['spectrums','radardata','fitting'],configfile = 'planeproc2.ini',remakealldata=False):
    """ This will run the radar simulations for all the selected data sets"""
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

def save2dropbox(testpath,imgonly=True):
    endpath = testpath.split(os.path.sep)[-1]
    imgpaths = ['Inputimages','fittedimages']
    dboxpath =os.path.expanduser(os.path.join('~','Dropbox'))
    dboxsave=os.path.join(dboxpath,'Planeexample',endpath)

    if os.path.exists(dboxpath):
        if os.path.exists(dboxsave):
            shutil.rmtree(dboxsave)
        for i in imgpaths:
            shutil.copytree(os.path.join(testpath,i),os.path.join(dboxsave,i))

if __name__== '__main__':
    argv = sys.argv[1:]

    outstr = 'Planeproc.py -f <function: origdata, spectrums, radardata, fitting or all> -i <basedir List or all> -p <testpath> -c <config> -r <type y to remake data>'
    curpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    try:
        opts, args = getopt.gnu_getopt(argv,"hf:i:c:r:p:")
    except getopt.GetoptError:
        print(outstr)
        sys.exit(2)

    remakealldata = False
    basedir = os.path.join(curpath,'exp_width_01')
    funcnamelist=[]
    for opt, arg in opts:
        if opt == '-h':
            print(outstr)
            sys.exit()
        elif opt in ('-p','--path'):
            curpath=os.path.expanduser(arg)
        elif opt in ("-i", "--ifile"):
            basedir = arg

        elif opt in ("-c", "--cfile"):
            outdirexist = True
            configfile = arg
        elif opt in ("-f", "--func"):
            funcnamelist.append(arg)

        elif opt in ('-r', "--re"):
            if arg.lower() == 'y':
                remakealldata = True


    if 'origdata' in funcnamelist:
        funcnamelist.remove('origdata')
        makedirs = True
        (sensdict,simparams) = readconfigfile(configfile)
        azangles = [iang[0] for iang in simparams['angles']]
        meanaz = sp.mean(azangles)
        makealldata(curpath,meanaz)

    if 'all' in funcnamelist:

        funcnamelist=['spectrums','radardata','fitting','plotting']


    if basedir.lower() == 'all':
        basedirlist = glob.glob(os.path.join(curpath,'exp_width_*'))
    else:
        basedirlist = basedir.split()
    plotbool = False
    if 'plotting' in funcnamelist:
        plotbool=True
        funcnamelist.remove('plotting')
    if len(funcnamelist)>0:
        for ibase in basedirlist:
            runradarsims(ibase,funcnamelist,configfile,remakealldata)
            #save2dropbox(ibase)
            if plotbool:
                plotoutdata(ibase,os.path.join(ibase,'Inputimages'))
                plotoutput(ibase,os.path.join(ibase,'fittedimages'))
    else:
        if plotbool:
            for ibase in basedirlist:
                plotoutdata(ibase,os.path.join(ibase,'Inputimages'))
                plotoutput(ibase,os.path.join(ibase,'fittedimages'))
