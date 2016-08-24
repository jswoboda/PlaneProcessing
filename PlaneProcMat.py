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
import cvxpy as cvx
from PlaneProc import makeline, runradarsims
from PlaneProcPlot import plotinputdata,plotoutput,plotoutputerrors,ploterrors




def invertRSTO(RSTO,Iono):
    
    
    nlout,ntout,np=Iono.Param_List.shape
    nlin=len(RSTO.Cart_Coords_In)
    time_out=RSTO.Time_Out
    time_in=RSTO.Time_In
    overlaps = RSTO.overlaps
    rin,azin,elin=RSTO.Sphere_Coords_In.transpose()
    
    anglist=RSTO.simparams['angles']
    ang_vec=sp.array([[i[0],i[1]] for i in anglist])
    
    # trim out cruft
    allrng= RSTO.simparams['Rangegatesfinal']
    dR=allrng[1]-allrng[0]
    npdir=sp.ceil(int(np)/2.)
    minangpos=ang_vec[sp.logical_and(ang_vec[:,0]<180.,ang_vec[:,0]>=0),1].min()
    rngbounds=[allrng[0]-npdir*dR,allrng[-1]+npdir*dR]
    rng_log=sp.logical_and(rin>rngbounds[0],rin<rngbounds[1])
    elbounds=elin>minangpos-2
    keeplog=sp.logical_and(rng_log,elbounds)
    keeplist=sp.where(keeplog)[0]
    nlin_red=len(keeplist)
    #make constriants
    alpha=1/(1e1)**2
    new_params=sp.zeros((nlin,len(time_in),np),dtype=Iono.Param_List.dtype)
    for itimen, itime in enumerate(time_out):
        print('Making Outtime {0:d} of {1:d}'.format(itimen+1,len(time_out)))
        allovers=overlaps[itimen]
        curintimes=[i[0] for i in allovers]
        for it_in_n,it in enumerate(curintimes):
            print('\t Making Intime {0:d} of {1:d}'.format(it_in_n+1,len(curintimes)))
            A=RSTO.RSTMat[itimen*nlout:(itimen+1)*nlout,it*nlin:(it+1)*nlin]
            Acvx=cvx.Constant(A[:,keeplist])
            for ip in range(1):
                print('\t\t Making Lag {0:d} of {1:d}'.format(ip+1,np))
                b=Iono.Param_List[:,itimen,ip]
                xr=cvx.Variable(nlin_red)
                xi=cvx.Variable(nlin_red)
                constr=alpha*cvx.norm(xr,2)
                consti=alpha*cvx.norm(xi,2)
                br=b.real
                bi=b.imag
                if ip==0:
                    objective=cvx.Minimize(cvx.norm(Acvx*xr-br,2)+constr)
                    constraints= [xr>=0]
                    prob=cvx.Problem(objective)
                    result=prob.solve(verbose=True,solver=cvx.ECOS)
                    new_params[keeplog,it,ip]=xr.value.flatten()
                else:
                    objective=cvx.Minimize(cvx.norm(Acvx*xr-br,2)+constr)
                    prob=cvx.Problem(objective)
                    result=prob.solve(verbose=True,solver=cvx.ECOS)
                    
                    objective=cvx.Minimize(cvx.norm(Acvx*xi-bi,2)+consti)
                    prob=cvx.Problem(objective)
                    result=prob.solve(verbose=True,solver=cvx.ECOS)
                    xcomp=xr.value.flatten()+1j*xi.value.flatten()
                    new_params[keeplog,it,ip]=xcomp
    ionoout=IonoContainer(coordlist=RSTO.Cart_Coords_In,paramlist=new_params,times = time_in,sensor_loc = sp.zeros(3),ver =0,coordvecs =
        ['x','y','z'],paramnames=Iono.Param_Names)
    return ionoout
    
    
def diffmat(dims,order = 'C'):
    """ This function will return a tuple of difference matricies for data from an 
        Nd array that has been rasterized. The order parameter determines whether 
        the array was rasterized in a C style (python) of FORTRAN style (MATLAB).
    """
    dims=[int(i) for i in dims]
    xdim = dims[0]
    ydim = dims[1]
    dims[0]=ydim
    dims[1]=xdim

    if order.lower() == 'c':
        dims = dims[::-1]

    outD = []
    for idimn, idim in enumerate(dims):
        if idim==0:
            outD.append(sp.array([]))
            continue
        e = sp.ones(idim)
        dthing = sp.vstack((-e,e))
        D = sp.sparse.spdiags(dthing,[0,1],idim-1,idim).toarray()
        D = sp.vstack((D,D[-1]))
        if idim>0:
            E = sp.sparse.eye(sp.prod(dims[:idimn]))
            D = sp.sparse.kron(D,E)

        if idimn<len(dims)-1:
            E = sp.sparse.eye(sp.prod(dims[idimn+1:]))
            D = sp.sparse.kron(E,D)

        outD.append(D)
    if order.lower() == 'c':
        outD=outD[::-1]
    Dy=outD[0]
    Dx = outD[1]
    outD[0]=Dx
    outD[1]=Dy

    return tuple(outD)
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

        funcnamelist=['spectrums','applymat','fittingmat','plotting']



    plotboolin = False
    plotboolout= False
    ploterror=False
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
            plotoutput(ibase,os.path.join(ibase,'fittedimagesmat'),configfile,wtimes,fitpath='FittedMat')
        if ploterror:
            plotoutputerrors(ibase,os.path.join(ibase,'fittedmaterrorimages'),configfile,wtimes,fitpath='FittedMat')
            ploterrors(ibase,os.path.join(ibase,'fittederroronlyimages'),configfile,wtimes,fitpath='FittedMat')
            #save2dropbox(ibase)
