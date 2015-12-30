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
import matplotlib.pyplot as plt

from PlaneProcPlot import plotoutdata,plotoutput



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
    plotoutput(testpath,os.path.join(testpath,'fittedimages'))

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


    if 'origdata' in funcnamelist:
        funcnamelist.remove('origdata')
        makedirs = True
        (sensdict,simparams) = readconfigfile(configfile)
        azangles = [iang[0] for iang in simparams['angles']]
        meanaz = sp.mean(azangles)
        makealldata(curpath,meanaz)

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
