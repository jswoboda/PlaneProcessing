#!/usr/bin/env python
"""
Created on Wed Dec 30 16:20:38 2015

@author: John Swoboda
"""
import os,glob,shutil
import numpy as np
import numpy.ma as ma
import scipy as sp
import scipy.io as sio
import pdb
from RadarDataSim.IonoContainer import IonoContainer,MakeTestIonoclass
from RadarDataSim.utilFunctions import Chapmanfunc

def makesimpledata(inputfile,timevec= None,begx=0.,begz=300.,vx=500.):
    if timevec is None:
        timevec= sp.column_stack((sp.linspace(0,60.,4,False),sp.linspace(0,60.,4,False)+15.))
    d2r=sp.pi/180.
    Iono=IonoContainer.readh5(inputfile)
    x,y,z=Iono.Cart_Coords.transpose()
    H_0=50.
    z_0=begz
    N_0=1e11
    basefunc=Chapmanfunc(z,H_0,z_0,N_0)
    angr= d2r*Iono.Sphere_Coords[:,1].max()
    r=sp.sqrt(x**2+y**2)*sp.sign(x)
    nx=r.size
    np=1
    enh=2.
    sigx=30.
    sigz=30.
    outdata=sp.zeros((nx,len(timevec),np))
    vel=sp.zeros((nx,len(timevec),3))
    for it,t in enumerate(timevec): 
        xc=begx+vx*t[0]*1e-3
        xwin=sp.absolute(r-xc)<sigx
        xwin=xwin.astype(x.dtype)
#        zwin=sp.absolute(z-begz)<sigz
#        zwin=zwin.astype(z.dtype)
#        x2=sp.exp(2-1/(1-(r-xc)**2/sigx**2)-1/(1-(z-begz)**2/sigz**2))
        x2=sp.exp(1-1/(1-(r-xc)**2/sigx**2))
        xar=ma.array(x2,mask=sp.logical_not(xwin),fill_value=0.0)
        outdata[:,it,0]=(enh*xar.filled() +1.)*basefunc
        vr=vx*xwin*sp.ones_like(xwin)
        vel[:,it,0]=vr*sp.sin(angr)
        vel[:,it,1]=vr*sp.cos(angr)
    ionoout=IonoContainer(coordlist=Iono.Cart_Coords,paramlist=outdata,times = timevec,sensor_loc = sp.zeros(3),ver =0,coordvecs =
        ['x','y','z'],paramnames=sp.array(['0']),velocity=vel)
    return ionoout
    
        
def changefilenames(folder,exten,inttime,filetemplate,folder2=None):
    if folder2 is None:
        folder2=folder
    if not os.path.isdir(folder2):
        os.mkdir(folder2)
    
    files = glob.glob(os.path.join(folder,'*.'+exten))
    fnames = [os.path.split(i)[-1] for i in files]
    fnums = [int(i.split('_')[0]) for i in fnames]
    times = [i*inttime for i in fnums]
    fnumsort = [i for (v, i) in sorted((v, i) for (i, v) in enumerate(fnums))]
    maxtime = max(times)
    fmtstr ='{0:0>'+ str(int(np.floor(np.log10(maxtime))+1))+'}'
    newnames = [fmtstr.format(i)+filetemplate+'.'+exten for i in times]

    for inum in fnumsort:
        curfile = files[inum]
        newfilename = os.path.join(folder2,newnames[inum])
        shutil.copy2(curfile,newfilename)

#%% Translation and speed up
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
#%% Converting files
def convertMattsfiles(filename,datadir,outdir,keepspec=[0,1,2,6],angle=15.,offset=0.):
    """ 
    This function will convert a set of files from Matt Zettegrens simulator to formated h5 files that 
    RadarDataSim can read. 
    Inputs
        filename - This is a .mat file that will be converted.
        datadir - This is the directory location of the file.
        outdir - The directory that the data will be stored.
        angle - The angle in degrees that the thing will lie on.
        keepspec - A numpy array of numbers.
        offset - Offset of the x-axis.
        outdir - The directory this is getting saved.
    """
    d2r=sp.pi/180.
    angr=d2r*angle
    lsp=7
    # Read in Data
    
    inst = sio.loadmat(os.path.join(datadir,filename))
    xg=inst['xg'][0,0]
    x1v = xg['xp']+offset# added to avoid gratting lobes.
    x3v = xg['zp']
    
    [x1mat,x3mat] = sp.meshgrid(x1v,x3v);
    

    print('Reading {0}'.format(filename))
    E = x1mat*sp.sin(angr)#x
    N = x1mat*sp.cos(angr)#y
    U = x3mat
    lxs=x3mat.size
    Time_Vector = sp.column_stack([inst['t'],inst['t']+15])
    ns =inst['ns']
    print('\tLoaded densities...');
    
    ns= sp.reshape(ns,[lxs,lsp])
    Ts =inst['Ts']
    print('\tLoaded temperatures...')
    
    Ts=sp.reshape(Ts,[lxs,lsp])
    vs = inst['vsx1']
    print('\tLoaded parallel velocities...\n');
    
    # derive velocity from ExB
    Ez,Ex=sp.gradient(-1*inst['Phi'])
    dx=sp.diff(xg['x'].flatten())[0]
    dEdx=Ex/dx
    vx1=-1*dEdx/xg['Bmag']
    # from looking at the data it seems that the velocity is off by a factor of
    # 10.
    vx1=vx1.flatten()/10.
    vs=sp.reshape(vs,[lxs,lsp])
    vs=sp.sum(ns[:,:(lsp-1)]*vs[:,:(lsp-1)],1)/ns[:,lsp-1]
    v_e= vx1*sp.sin(angr)
    v_n = vx1*sp.cos(angr)
    v_u = vs
    #change units of velocity to km/s
    Velocity = sp.reshape(sp.column_stack([v_e,v_n,v_u]),[lxs,1,3])
    # reduce the number of spcecies
#    if islogical(keepspec)
#        keepspec(end)=true;
#        keepspecnum = sum(keepspec);
#    elseif ~any(keepspec==numspec)
#        keepspec = [keepspec(:),numspec];
#        keepspecnum = length(keepspec);
#    else
#        keepspecnum = length(keepspec);
#    end
    nsout = sp.reshape(ns[:,keepspec],[lxs,1,len(keepspec)])
    Tsout = sp.reshape(Ts[:,keepspec],[lxs,1,len(keepspec)])
    
    
    
    # Put everything in to ionocontainer structure
    Cart_Coords = sp.column_stack([E.flatten(),N.flatten(),U.flatten()])*1e-3

    Param_List = sp.concatenate((sp.expand_dims(nsout,nsout.ndim),sp.expand_dims(Tsout,Tsout.ndim)),-1);
    Species = sp.array(['O+','NO+','N2+','O2+','N+', 'H+','e-'])
    Species = Species[keepspec]
    Species= [i for i in Species]
    fpart=os.path.splitext(filename)[0].split('_')[-1]
    fout=os.path.join(outdir,fpart+'.h5')
    ionoout=IonoContainer(Cart_Coords,Param_List,Time_Vector,ver=0,species=Species,velocity=Velocity)
    ionoout.Sphere_Coords[:,1]=ionoout.Sphere_Coords[:,1].astype('float32')
    ionoout.Sphere_Coords[:,1]=ionoout.Sphere_Coords[:,1].astype('float64')
    ionoout.saveh5(fout)
def convertMattsfiles_dir(datadir,outdir,keepspec=[0,1,2,6],angle=15.,offset=0.):
    fnames=glob.glob(os.path.join(datadir,'2*.mat'))
    filenames=[os.path.split(i)[-1] for i in fnames]
    for i in filenames:
        convertMattsfiles(i,datadir,outdir,keepspec=keepspec,angle=angle,offset=offset)
    
if __name__== '__main__':

    from argparse import ArgumentParser
    descr = '''
             This script will perform the basic run est for ISR sim.
            '''
    p = ArgumentParser(description=descr)
    p.add_argument('-f','--fold',help='Folder',default='~/Dropbox/PerryPlane/fittedimages240seconds/')
    p.add_argument('-p','--path2',help='Output Folder',default='~/Dropbox/PerryPlane/fittedimages240seconds2/')
    p.add_argument('-i','--inttime',help='Integration Time',type=int,default = 240)
    
    args = p.parse_args()
    folder = os.path.expanduser(args.fold)
    folder2 = os.path.expanduser(args.path2)
    inttime = args.inttime
    filetemplate = "_{0}_int".format(inttime)
    exten= 'png'
    changefilenames(folder,exten,inttime,filetemplate,folder2)
