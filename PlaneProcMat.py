#!/usr/bin/env python
"""
Created on Wed Dec 30 16:11:56 2015

@author: John Swoboda
"""

import os, glob,inspect,getopt,sys
import shutil
import pdb
import scipy as sp
import numbers
import matplotlib
import pickle
matplotlib.use('Agg')
from RadarDataSim.IonoContainer import IonoContainer, MakeTestIonoclass,makeionocombined
import RadarDataSim.runsim as runsim
from RadarDataSim.radarData import makeCovmat
from RadarDataSim.analysisplots import analysisdump
from RadarDataSim.utilFunctions import readconfigfile,spect2acf
from RadarDataSim.operators import RadarSpaceTimeOperator
import matplotlib.pyplot as plt
import scipy.fftpack as scfft
import cvxpy as cvx
from PlaneProc import makeline, runradarsims
from PlaneProcPlot import plotinputdata,plotoutput,ploterrors,plotalphaerror,plotLcurve




def invertRSTO(RSTO,Iono,alpha_list=1e-2,invtype='tik',rbounds=[100,200]):
    """ """
    
    nlout,ntout,np=Iono.Param_List.shape
    nlin=len(RSTO.Cart_Coords_In)
    time_out=RSTO.Time_Out
    time_in=RSTO.Time_In
    overlaps = RSTO.overlaps
    xin,yin,zin=RSTO.Cart_Coords_In.transpose()
    z_u=sp.unique(zin)
    rplane=sp.sqrt(xin**2+yin**2)*sp.sign(xin)
    r_u=sp.unique(rplane)
    n_z=z_u.size
    n_r=r_u.size
    dims= [n_r,n_z]
    
    rin,azin,elin=RSTO.Sphere_Coords_In.transpose()
    
    anglist=RSTO.simparams['angles']
    ang_vec=sp.array([[i[0],i[1]] for i in anglist])
    
    # trim out cruft
    
    zmin,zmax=[150,500]
    rpmin,rpmax=rbounds#[-50,100]#[100,200]
    altlog= sp.logical_and(zin>zmin,zin<zmax)
    rplog=sp.logical_and(rplane>rpmin,rplane<rpmax)
    allrng= RSTO.simparams['Rangegatesfinal']
    dR=allrng[1]-allrng[0]
    npdir=sp.ceil(int(np)/2.)
    posang_log1= sp.logical_and(ang_vec[:,0]<=180.,ang_vec[:,0]>=0)
    negang_log1 = sp.logical_or(ang_vec[:,0]>180.,ang_vec[:,0]<0)
    azin_pos = sp.logical_and(azin<=180.,azin>=0)
    azin_neg = sp.logical_or(azin>180.,azin<0)
    minangpos=0
    minangneg=0
    
    
    if sp.any(posang_log1):
        minangpos=ang_vec[posang_log1,1].min()
    if sp.any(negang_log1):
        minangneg=ang_vec[negang_log1,1].min()
    
    rngbounds=[allrng[0]-npdir*dR,allrng[-1]+npdir*dR]
    rng_log=sp.logical_and(rin>rngbounds[0],rin<rngbounds[1])
    elbounds_pos=sp.logical_and(azin_pos,elin>minangpos)
    elbounds_neg=sp.logical_and(azin_neg,elin>minangneg)
    
    elbounds=sp.logical_or(elbounds_pos,elbounds_neg)
    keeplog=sp.logical_and(sp.logical_and(rng_log,elbounds),sp.logical_and(altlog,rplog))
    keeplist=sp.where(keeplog)[0]
    nlin_red=len(keeplist)
    # set up derivative matrix
    dx,dy=diffmat(dims)
    dx_red=dx[keeplist][:,keeplist]
    dy_red=dy[keeplist][:,keeplist]
    # need the sparse vstack to make srue things stay sparse
    D=sp.sparse.vstack((dx_red,dy_red))
    # New parameter matrix
    new_params=sp.zeros((nlin,len(time_out),np),dtype=Iono.Param_List.dtype)
    if isinstance(alpha_list,numbers.Number):
        alpha_list=[alpha_list]*np
    ave_datadif=sp.zeros((len(time_out),np))
    ave_data_const = sp.zeros_like(ave_datadif)
    q=1e10
    for itimen, itime in enumerate(time_out):
        print('Making Outtime {0:d} of {1:d}'.format(itimen+1,len(time_out)))
        #allovers=overlaps[itimen]
        #curintimes=[i[0] for i in allovers]
        #for it_in_n,it in enumerate(curintimes):
        #print('\t Making Intime {0:d} of {1:d}'.format(it_in_n+1,len(curintimes)))
        #A=RSTO.RSTMat[itimen*nlout:(itimen+1)*nlout,it*nlin:(it+1)*nlin]
        A=RSTO.RSTMat[itimen*nlout:(itimen+1)*nlout,itimen*nlin:(itimen+1)*nlin]
        Acvx=cvx.Constant(A[:,keeplist])
        for ip in range(np):
            alpha=alpha_list[ip]
            print('\t\t Making Lag {0:d} of {1:d}'.format(ip+1,np))
            datain=Iono.Param_List[:,itimen,ip]/q
            xr=cvx.Variable(nlin_red)
            xi=cvx.Variable(nlin_red)
            if invtype.lower()=='tik':
                constr=alpha*cvx.norm(xr,2)
                consti=alpha*cvx.norm(xi,2)
            elif invtype.lower()=='tikd':
                constr=alpha*cvx.norm(D*xr,2)
                consti=alpha*cvx.norm(D*xi,2)
            elif invtype.lower()=='tv':
                constr=alpha*cvx.norm(D*xr,1)
                consti=alpha*cvx.norm(D*xi,1)
            br=datain.real
            bi=datain.imag
            if ip==0:
                objective=cvx.Minimize(cvx.norm(Acvx*xr-br,2)+constr)
                constraints= [xr>=0]
                prob=cvx.Problem(objective)
                result=prob.solve(verbose=True,solver=cvx.SCS,use_indirect=True)
#                    new_params[keeplog,it,ip]=xr.value.flatten()
                xcomp=sp.array(xr.value).flatten()*q
            else:
                objective=cvx.Minimize(cvx.norm(Acvx*xr-br,2)+constr)
                prob=cvx.Problem(objective)
                result=prob.solve(verbose=True,solver=cvx.SCS,use_indirect=True)
                
                objective=cvx.Minimize(cvx.norm(Acvx*xi-bi,2)+consti)
                prob=cvx.Problem(objective)
                result=prob.solve(verbose=True,solver=cvx.SCS,use_indirect=True)
                xcomp=sp.array(xr.value + 1j*xi.value).flatten()*q
#                    new_params[keeplog,it,ip]=xcomp
            new_params[keeplog,itimen,ip]=xcomp
            ave_datadif[itimen,ip]=sp.sqrt(sp.nansum(sp.absolute(A[:,keeplist].dot(xcomp)-datain*q)**2))
            if invtype.lower()=='tik':
                sumconst=sp.sqrt(sp.nansum(sp.power(sp.absolute(xcomp),2)))
            elif invtype.lower()=='tikd':
                dx=D.dot(xcomp)
                sumconst=sp.sqrt(sp.nansum(sp.power(sp.absolute(dx),2)))
            elif invtype.lower()=='tv':
                dx=D.dot(xcomp)
                sumconst=sp.nansum(sp.absolute(dx))
            ave_data_const[itimen,ip]=sumconst
            # set up nans                    
            new_params[sp.logical_not(keeplog),itimen]=sp.nan
    datadif=sp.nanmean(ave_datadif,axis=0)
    constval=sp.nanmean(ave_data_const,axis=0)
    ionoout=IonoContainer(coordlist=RSTO.Cart_Coords_In,paramlist=new_params,times = time_out,sensor_loc = sp.zeros(3),ver =0,coordvecs =
        ['x','y','z'],paramnames=Iono.Param_Names)
        
    return (ionoout,datadif,constval)
    
def runinversion(basedir,configfile,acfdir='ACF',invtype='tik',alpha=1e-2):
    """ """
    costdir = os.path.join(basedir,'Cost')
    
    pname=os.path.join(costdir,'cost{0}-{1}.pickle'.format(acfdir,invtype))
    pickleFile = open(pname, 'rb')
    alpha_arr=pickle.load(pickleFile)[-1]
    pickleFile.close()
    
    ionoinfname=os.path.join(basedir,acfdir,'00lags.h5')
    ionoin=IonoContainer.readh5(ionoinfname)
    
    dirio = ('Spectrums','Mat','ACFMat')
    inputdir = os.path.join(basedir,dirio[0])
    
    dirlist = glob.glob(os.path.join(inputdir,'*.h5'))
    (listorder,timevector,filenumbering,timebeg,time_s) = IonoContainer.gettimes(dirlist)
    Ionolist = [dirlist[ikey] for ikey in listorder]
    if acfdir.lower()=='acf':
        
        acfloc='ACFInv'
    elif acfdir.lower()=='acfmat':
        mattype='matrix'
        acfloc='ACFMatInv'
    mattype='sim'
    RSTO = RadarSpaceTimeOperator(Ionolist,configfile,timevector,mattype=mattype)  
    if 'perryplane' in basedir.lower():
        rbounds=[-500,500]
    else:
        rbounds=[0,500]
    
    ionoout=invertRSTO(RSTO,ionoin,alpha_list=alpha_arr,invtype=invtype,rbounds=rbounds)[0]
    outfile=os.path.join(basedir,acfloc,'00lags{0}.h5'.format(invtype))
    ionoout.saveh5(outfile)
    if acfdir=='ACF':
        lagsDatasum=ionoout.Param_List
        # !!! This is done to speed up development 
        lagsNoisesum=sp.zeros_like(lagsDatasum)
        Nlags=lagsDatasum.shape[-1]
        pulses_s=RSTO.simparams['Tint']/RSTO.simparams['IPP']
        Ctt=makeCovmat(lagsDatasum,lagsNoisesum,pulses_s,Nlags)
        outfile=os.path.join(basedir,acfloc,'00sigs{0}.h5'.format(invtype))
        ionoout.Param_List=Ctt
        ionoout.Param_Names=sp.repeat(ionoout.Param_Names[:,sp.newaxis],Nlags,axis=1)
        ionoout.saveh5(outfile)
        
def mkalphalist(pnamefile):
    
    pickleFile = open(pnamefile, 'rb')
    dictlist = pickle.load(pickleFile)
    alpha_list,errorlist,datadif,constdif,errorlaglist=dictlist[:5]
    
    pickleFile.close()
    os.remove(pnamefile)

    errorlagarr=sp.array(errorlaglist)
    alphar=sp.array(alpha_list)
    errlocs=sp.argmin(errorlagarr,axis=0)
    alout=alphar[errlocs]
    pickleFile = open(pnamefile, 'wb')
    pickle.dump([alpha_list,errorlist,datadif,constdif,errorlaglist,alout],pickleFile)
    pickleFile.close()
    
def parametersweep(basedir,configfile,acfdir='ACF',invtype='tik'):
    """ This function will run the inversion numerious times with different constraint
        parameters. This will create a directory called cost and place.
        Input
        basedir - The directory that holds all of the data for the simulator.
        configfile - The ini file for the simulation.
        acfdir - The directory within basedir that hold the acfs to be inverted.
        invtype - The inversion method that will be tested. Can be tik, tikd, and tv.
        """

    alpha_sweep=sp.logspace(-3.5,3.5,25)
    costdir = os.path.join(basedir,'Cost')
    ionoinfname=os.path.join(basedir,acfdir,'00lags.h5')
    ionoin=IonoContainer.readh5(ionoinfname)
    
    dirio = ('Spectrums','Mat','ACFMat')
    inputdir = os.path.join(basedir,dirio[0])
    
    dirlist = glob.glob(os.path.join(inputdir,'*.h5'))
    (listorder,timevector,filenumbering,timebeg,time_s) = IonoContainer.gettimes(dirlist)
    Ionolist = [dirlist[ikey] for ikey in listorder]
    
    RSTO = RadarSpaceTimeOperator(Ionolist,configfile,timevector,mattype='Sim')
    
    
    ionospec=makeionocombined(dirlist)
    tau,acfin=spect2acf(ionospec.Param_Names,ionospec.Param_List)
    nloc,ntimes=acfin.shape[:2]
    
    # get the original acf
    ambmat=RSTO.simparams['amb_dict']['WttMatrix']
    np=ambmat.shape[0]
    acfin_amb=sp.zeros((nloc,ntimes,np),dtype=acfin.dtype)
        
    for iloc,locarr in enumerate(acfin):
        for itime,acfarr in enumerate(locarr):
            acfin_amb[iloc,itime]=sp.dot(ambmat,acfarr)
    acfin_amb=acfin_amb[:,0]
    if not os.path.isdir(costdir):
        os.mkdir(costdir)
    # pickle file stuff 
    pname=os.path.join(costdir,'cost{0}-{1}.pickle'.format(acfdir,invtype))

    alpha_list=[]
    errorlist=[]
    errorlaglist=[]
    datadiflist=[]
    constlist=[]
    
    alpha_list_new=alpha_sweep.tolist()
    for i in alpha_list:
        if i in alpha_list_new:
            alpha_list_new.remove(i)
    
    for i in alpha_list_new:
        ionoout,datadif,constdif=invertRSTO(RSTO,ionoin,alpha_list=i,invtype=invtype)
        
        datadiflist.append(datadif)
        constlist.append(constdif)
        acfout=ionoout.Param_List[:,0]
        alpha_list.append(i)
        outdata=sp.power(sp.absolute(acfout-acfin_amb)/sp.absolute(acfin_amb),2)
        aveerror=sp.nanmean(outdata,axis=0)
        errorlaglist.append(aveerror)
        errorlist.append(sp.nansum(aveerror))
        
    pickleFile = open(pname, 'wb')
    pickle.dump([alpha_list,errorlist,datadiflist,constlist,errorlaglist],pickleFile)
    pickleFile.close()
    mkalphalist(pname)
    alphaarr=sp.array(alpha_list)
    errorarr=sp.array(errorlist)
    errorlagarr=sp.array(errorlaglist)
    datadif=sp.array(datadiflist)
    constdif=sp.array(constlist)
    fig,axlist,axmain=plotalphaerror(alphaarr,errorarr,errorlagarr)
    fig.savefig(os.path.join(costdir,'lcurve{0}-{1}.png'.format(acfdir,invtype)))
    
    fig,axlist=plotLcurve(alphaarr,datadif,constdif)
    fig.savefig(os.path.join(costdir,'lcurve{0}-{1}.png'.format(acfdir,invtype)))
    
def diffmat(dims,order = 'C'):
    """ This function will return a tuple of difference matricies for data from an 
        Nd array that has been rasterized. The order parameter determines whether 
        the array was rasterized in a C style (python) of FORTRAN style (MATLAB).
        Inputs:
            dims- A list of the size of the x,y,z.. dimensions.
            order- Specifies the vectorization of the matrix
        Outputs:
            dx,dy,dy... - The finite difference operators for a vectorized array.
                If these are to be stacked together as one big operator then
                sp.sparse.vstack should be used.
    """
    # flip the dimensions around
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

        outD.append(sp.sparse.csc_matrix(D))
    if order.lower() == 'c':
        outD=outD[::-1]
    Dy=outD[0]
    Dx = outD[1]
    outD[0]=Dx
    outD[1]=Dy
    return tuple(outD)
    
def cgmat(A,x,b,M=None,max_it=100,tol=1e-8):
    """ This function will performa conjuguate gradient search to find the inverse of
        an operator A, given a starting point x, and data b.
    """
    if M is None:
        M= sp.diag(A)
    bnrm2 = sp.linalg.norm(b)
    r=b-A.dot(x)
    rho=sp.zeros(max_it)
    for i in range(max_it):
        z=sp.linalg.solve(M,r)
        rho[i] = sp.dot(r,z)
        if i==0:
            p=z
        else:
            beta=rho/rho[i-1]
            p=z+beta*p

        q=A.dot(p)
        alpha=rho/sp.dot(p,q)
        x = x+alpha*p
        r = r-alpha*q
        error = sp.linalg.norm( r ) / bnrm2
        if error <tol:
            return (x,error,i,False)

    return (x,error,max_it,True)

    
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
    p.add_argument('-k', "--ktype",help='The type of constraint can be tik tikd tv.',default='')
    p.add_argument('-a', "--acftype",help='The ACF directory that will have the inversions applied to it',default='ACFMat')
    p.add_argument('-g', "--gamma",help='The Parameter gamma for the constraints',default=1e-2)
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
    acffolder=args.acftype
    invtype=args.ktype
    gamma=float(args.gamma)
    
        
    if len(fittimes)==0:
        fittimes=None
    else:
        fittimes = [int(i) for i in fittimes]

    configlist = ['planeproc2.ini','planeproc2_stat.ini','dishplaneproc.ini','dishplaneproc_stat.ini']
   
    if basedir.lower() == 'all':
        basedirlist = glob.glob(os.path.join(curpath,'exp_width_*'))
    else:
        basedirlist = basedir.split()

    if 'paramsweep' in funcnamelist:
        
        parametersweep(basedir,configfile,acfdir=acffolder,invtype=invtype)
        funcnamelist.remove('paramsweep')
    if 'invertdata' in funcnamelist:
        runinversion(basedir,configfile,acfdir=acffolder,invtype=invtype,alpha=gamma)
        funcnamelist.remove('invertdata')
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
    plotmat=False
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
    if 'plottingmat' in funcnamelist:
        plotmat=True
        funcnamelist.remove('plottingmat')
    for ibase in basedirlist:
        if len(funcnamelist)>0:
            runradarsims(ibase,funcnamelist,configfile,remakealldata,fittimes,invtype=invtype)
            #save2dropbox(ibase)
        if plotboolin:
            plotinputdata(ibase,os.path.join(ibase,'Inputimages'),wtimes)
        if plotboolout:
            #
            plotoutput(ibase,os.path.join(ibase,'fittedimages{}'.format(invtype)),configfile,wtimes,fitpath='FittedInv',fitfile='fitteddata{0}.h5'.format(invtype))
        if plotmat:
            plotoutput(ibase,os.path.join(ibase,'fittedimagesmat'),configfile,wtimes,fitpath='FittedMat')
        if ploterror:
            ploterrors(ibase,os.path.join(ibase,'fittederroronlyimages'),configfile,wtimes,fitpath='FittedMat')
            #save2dropbox(ibase)
