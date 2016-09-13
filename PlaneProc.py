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
from RadarDataSim.analysisplots import analysisdump,plotbeamparametersv2
from RadarDataSim.utilFunctions import readconfigfile
from PlaneProcPlot import plotinputdata,plotoutput,ploterrors
from RadarDataSim.IonoContainer import IonoContainer


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
    Xmat = Rmat*sp.cos(d2r*meanaz)
    Ymat = Rmat*sp.sin(d2r*meanaz)
    coords = sp.column_stack((Xmat.flatten(),Ymat.flatten(),Zmat.flatten()))
    if rng_vel==0.:
        timevec = sp.array([0.])
    else:
        timevec = sp.linspace(0.,900.,nt)

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
def translatemat(datadir,outputdir):
    flist=glob.glob(os.path.join(datadir,'V*.mat'))
    f2list=[os.path.splitext(os.path.split(i)[-1])[0].split('_')[-1] for i in flist]
    for i,ifile in enumerate(flist):
        print 'Opening {}'.format(os.path.split(ifile)[-1])
        curiono=IonoContainer.readmat(ifile)
        curfile=os.path.join(outputdir,f2list[i]+'.h5')
        curiono.saveh5(curfile)

def speedup(datadir,outputdir,s=5.):
    ext='.h5'
    ionocontlist = glob.glob(os.path.join(datadir,'*'+ext))
    (sortlist,outime,fileslist,timebeg,timelist_s)=IonoContainer.gettimes(ionocontlist)
    
    mintime=timebeg.min()
    
    for i in ionocontlist:
        curiono=IonoContainer.readh5(i)
        curiono.Time_Vector=(curiono.Time_Vector.astype(float)-mintime)/s
        curiono.Velocity=curiono.Velocity*s
        fname='{0:05d}.h5'.format(int(curiono.Time_Vector[0,0]))
        curfile=os.path.join(outputdir,fname)
        curiono.saveh5(curfile)
#%% 
def runradarsims(testpath,funcnamelist=['spectrums','radardata','fitting'],configfile = 'planeproc2.ini',remakealldata=False,fittimes=None,invtype=''):
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
#    folderlist =[]
#    if 'spectrums' in funcnamelist:
#        folderlist.append('Spectrums')
#    if 'radardata' in funcnamelist and remakealldata:
#        folderlist.append('Radardata')
#        folderlist.append('ACF')
#    if 'fitting' in funcnamelist:
#        folderlist.append('Fitted')
##    folderlist = ['Spectrums','Radardata','ACF','Fitted']
#    for ifl in folderlist:
#        flist = glob.glob(os.path.join(testpath,ifl,'*.h5'))
#        for ifile in flist:
#            os.remove(ifile)
    
    ismat = False
    for ifc in funcnamelist:
        if 'mat' in ifc:
            ismat=True
            break
    runsim.main(funcnamelist,testpath,configfile,remakealldata,fittimes,invtype=invtype)
    try:
        if ismat:
            plotdir = os.path.join(testpath,'AnalysisPlots')
            f_templ = os.path.join(plotdir,'paramsmat')
            plotbeamparametersv2([0.],configfile,testpath,fitdir = 'FittedMat',params=['Ne','Ti','Te'],filetemplate=f_templ,
                             suptitle = 'With Mat',werrors=False,nelog=False)
        else:
            analysisdump(testpath,configfile,'Plane Example')
    except:
        print "Analysis dump failed somewhere"

def save2dropbox(testpath,imgonly=True):
    endpath = testpath.split(os.path.sep)[-1]
    imgpaths = ['Inputimages','fittedimages','fittederrorimages']
    dboxpath =os.path.expanduser(os.path.join('~','Dropbox'))
    dboxsave=os.path.join(dboxpath,'Planeexample',endpath)

    if os.path.exists(dboxpath):
        if os.path.exists(dboxsave):
            shutil.rmtree(dboxsave)
        for i in imgpaths:
            shutil.copytree(os.path.join(testpath,i),os.path.join(dboxsave,i))
#%% Fixplanes
def fixspecs(basedirlist):
    for ibase in basedirlist:
        filelist = glob.glob(os.path.join(ibase,'Origparams','*.h5'))
        numlist = [os.path.splitext(os.path.split(x)[-1])[0] for x in filelist]
        numdict = {numlist[i]:filelist[i] for i in range(len(filelist))}
        slist = sorted(numlist,key=ke)

        origlist = [numdict[i] for i in slist]

        filelist = glob.glob(os.path.join(ibase,'Spectrums','*.h5'))
        numlist = [os.path.splitext(os.path.split(x)[-1])[0] for x in filelist]
        numdict = {numlist[i]:filelist[i] for i in range(len(filelist))}
        slist = sorted(numlist,key=ke)

        speclist = [numdict[i] for i in slist]
        for (iorig,ispec) in zip(origlist,speclist):
            origiono=IonoContainer.readh5(iorig)
            speciono=IonoContainer.readh5(ispec)
            speciono.Cart_Coords=origiono.Cart_Coords
            speciono.Sphere_Coords=origiono.Sphere_Coords
            os.remove(ispec)
            speciono.saveh5(ispec)
if __name__== '__main__':

    from argparse import ArgumentParser
    descr = '''
             This script will perform the basic run est for ISR sim.
            '''
    curpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    p = ArgumentParser(description=descr)
    p.add_argument('-t','--times',help='Times, as locations in the output time vector array that will be fit.',nargs='+',default=[])
    p.add_argument('-i','--idir',help='Base directory',default='all')
    p.add_argument('-c','--config',help='Config file for simlation',default = 'planeproc2_stat.ini')
    p.add_argument('-r ','--re',help='Remake data True or False.',type=bool,default=False)
    p.add_argument("-p", "--path",help='Number of pulses.',default=curpath)
    p.add_argument("-l", "--linewid",help='Line Width in number of Samples.',default=1)
    p.add_argument('-m', "--mult", help="Multiplication of enhancement.", default=5.)
    p.add_argument('-w', "--wtimes", help="Put times at top of plots.",default='n')
    p.add_argument('-f','--funclist',help='Functions to be uses',nargs='+',default=['spectrums','radardata','fitting'])#action='append',dest='collection',default=['spectrums','radardata','fitting','analysis'])
    
    args = p.parse_args()
    basedir = args.idir
    curpath = args.path
    configfile = args.config
    remakealldata = args.re
    funcnamelist = args.funclist
    fittimes = args.times
    lw =float( args.linewid)
    mult = float(args.mult)
    wtimes=args.wtimes.lower()=='y'

    if len(fittimes)==0:
        fittimes=None
    else:
        fittimes = [int(i) for i in fittimes]

    configlist = ['planeproc2.ini','planeproc2_stat.ini','dishplaneproc.ini','dishplaneproc_stat.ini']
   
    if basedir.lower() == 'all':
        basedirlist = glob.glob(os.path.join(curpath,'exp_width_*'))
    else:
        basedirlist = basedir.split()


    if 'origdata' in funcnamelist:
        funcnamelist.remove('origdata')
        makedirs = True
        (sensdict,simparams) = readconfigfile(configfile)
        azangles = [iang[0] for iang in simparams['angles']]
        meanaz = sp.mean(azangles)
        for ibase in basedirlist:
            makeline(ibase,meanaz,linewidth=lw,multval = mult)

    if 'all' in funcnamelist:

        funcnamelist=['spectrums','radardata','fitting','plotting']



    plotboolin = False
    plotboolout= False
    ploterror = False
    if 'plotting' in funcnamelist:
        plotboolin=True
        plotboolout=True
        ploterror=True
        funcnamelist.remove('plotting')
    if 'plottingin' in funcnamelist:
        plotboolin=True
        funcnamelist.remove('plottingin')
    if 'plottingout' in funcnamelist:
        plotboolout=True
        funcnamelist.remove('plottingout')
    if 'plottingerror' in funcnamelist:
        ploterror=True
        funcnamelist.remove('plottingerror')
    for ibase in basedirlist:
        if len(funcnamelist)>0:
            runradarsims(ibase,funcnamelist,configfile,remakealldata,fittimes)
            #save2dropbox(ibase)
        if plotboolin:
            plotinputdata(ibase,os.path.join(ibase,'Inputimages'),wtimes)
        if plotboolout:
            plotoutput(ibase,os.path.join(ibase,'fittedimages'),configfile,wtimes)
        if ploterror:
            ploterrors(ibase,os.path.join(ibase,'fittederroronlyimages'),configfile,wtimes)

