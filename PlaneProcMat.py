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

from PlaneProc import makeline, runradarsims
from PlaneProcPlot import plotinputdata,plotoutput,plotoutputerrors







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
    p.add_argument('-f','--funclist',help='Functions to be uses',nargs='+',default=['spectrums','applymat','fittingmat'])#action='append',dest='collection',default=['spectrums','radardata','fitting','analysis'])
        
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
    if 'plotting' in funcnamelist:
        plotboolin=True
        plotboolout=True
        funcnamelist.remove('plotting')
    if 'plottingin' in funcnamelist:
        plotboolin=True
        funcnamelist.remove('plottingin')
    if 'plottingout' in funcnamelist:
        plotboolout=True
        funcnamelist.remove('plottingout')
    for ibase in basedirlist:
        if len(funcnamelist)>0:
            runradarsims(ibase,funcnamelist,configfile,remakealldata,fittimes)
            #save2dropbox(ibase)
        if plotboolin:
            plotinputdata(ibase,os.path.join(ibase,'Inputimages'),wtimes)
        if plotboolout:
            plotoutput(ibase,os.path.join(ibase,'fittedimages'),configfile,wtimes,fitpath='FittedMat')
            plotoutputerrors(ibase,os.path.join(ibase,'fittederrorimages'),configfile,wtimes,fitpath='FittedMat')
            #save2dropbox(ibase)
