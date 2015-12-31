#!/usr/bin/env python
"""
Created on Wed Dec 30 16:11:56 2015

@author: John Swoboda
"""

import os, glob,inspect,getopt,sys
import shutil
import pdb
import scipy as sp
import matplotlib
matplotlib.use('Agg')
from RadarDataSim.IonoContainer import IonoContainer, MakeTestIonoclass
import RadarDataSim.runsim as runsim
from RadarDataSim.analysisplots import analysisdump
from RadarDataSim.utilFunctions import readconfigfile
from RadarDataSim.operators import RadarSpaceTimeOperator
import matplotlib.pyplot as plt
import scipy.fftpack as scfft


from PlaneProcPlot import plotoutdata,plotoutput



def runradarsims(testpath,funcnamelist=['spectrums','radardata','fitting'],configfile = 'planeproc2.ini',remakealldata=False):
    """ This will run the radar simulations for all the selected data sets"""
    origparamsdir = os.path.join(testpath,'Origparams')
    if not os.path.exists(testpath):
        os.mkdir(testpath)
        print "Making a path for testdata at "+testpath
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

def matrixcalc(testpath,funcnamelist,configfile):

    (sensdict,simparams) = readconfigfile(configfile)
    ts = simparams['t_s']
    plens = simparams['pulselength']
    plen = sp.floor(plens/ts)
    paramvec = sp.arange(plen)*ts
    # clear everything out
    folderlist =[]
    if 'matrix' in funcnamelist:
        folderlist.append('Matrix_ACF')
        folderlist.append('ACFOrig')
        folderlist.append('ACFMat')
    if 'fittingmat' in funcnamelist:
        folderlist.append('Matrix_Fitted')
    if 'matrix' in funcnamelist or 'fittingmat' in funcnamelist:
        folderlist.append('Matrix')


    for ifl in folderlist:
        fulldir = os.path.join(testpath,ifl)
        if os.path.isdir(fulldir):
            flist = glob.glob(os.path.join(fulldir,'*.h5'))
            for ifile in flist:
                os.remove(ifile)
        else:
            os.mkdir(fulldir)


    if 'matrix' in funcnamelist:
        specpath = os.path.join(testpath,'Spectrums')
        acfinpath = os.path.join(testpath,'ACFOrig')
        acfoutpath = os.path.join(testpath,'ACFMat')
        specflist = glob.glob(os.path(specpath,'*.h5'))

        acflist = []
        timelist = []
        for ifile in specflist:
            fname = os.path.split(ifile)[-1]
            newfname = os.path.join(acfinpath,fname)
            acflist.append(newfname)
            curiono = IonoContainer.readh5(ifile)
            timelist.append(curiono.Time_Vector)
            outdata = curiono.Param_List
            outdata = scfft.ifftshift(outdata,axes=-1)
            outdata = scfft.ifft(outdata,axes=-1)
            outdata = outdata[:,:,:plen]
            curiono.Param_List=outdata
            curiono.Param_Names=paramvec
            curiono.saveh5(newfname)
        timelist = sp.vstack(timelist)
        tidx = sp.argsort(timelist[:,0])
        timelist = timelist[tidx]
        acflist=acflist[tidx]
        RSTO = RadarSpaceTimeOperator(curiono,configfile,timelist)
        outiono = RSTO.mult_iono(acflist)
        outiono.saveh5(os.path.join(acfoutpath,'00 acfmat.h5'))

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
            curpath=arg
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




    if 'matrix' in funcnamelist:
        funcnamelist.remove('origdata')
        makedirs = True
        (sensdict,simparams) = readconfigfile(configfile)
        azangles = [iang[0] for iang in simparams['angles']]
        meanaz = sp.mean(azangles)

    if 'all' in funcnamelist:

        funcnamelist=['spectrums','radardata','fitting']


    if basedir.lower() == 'all':
        basedirlist = glob.glob(os.path.join(curpath,'exp_width_*'))
    else:
        basedirlist = basedir.split()

    if len(funcnamelist)>0:
        for ibase in basedirlist:
            runradarsims(ibase,funcnamelist,configfile,remakealldata)
            #save2dropbox(ibase)
