#!/usr/bin/env python
"""
Created on Wed Dec 30 16:20:38 2015

@author: John Swoboda
"""
import os,glob,shutil
import numpy as np
import scipy as sp
import scipy.io as sio
from RadarDataSim.IonoContainer import IonoContainer


def makesimpledata():
    
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


def convertMattsfiles(filename,datadir,outdir,keepspec=[0,1,2,6],angle=20.5,offset=0.):
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
    

    
    E = x1mat*sp.sin(angr)#x
    N = x1mat*sp.cos(angr)#y
    U = x3mat
    lxs=x3mat.size
    
    Time_Vector = sp.column_stack([inst['t'],inst['t']+15])
    ns =inst['ns']
    print('Loaded densities...');
    
    ns= sp.reshape(ns,[lxs,lsp])
    Ts =inst['Ts']
    print('Loaded temperatures...')
    
    Ts=sp.reshape(Ts,[lxs,lsp])
    vs = inst['vsx1']
    print('Loaded parallel velocities...\n');
    
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
    fpart=os.path.splitext(filename)[0]
    fout=os.path.join(outdir,fpart+'.h5')
    ionoout=IonoContainer(Cart_Coords,Param_List,Time_Vector,ver=0,species=Species,velocity=Velocity)
    
    ionoout.saveh5(fout)
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
