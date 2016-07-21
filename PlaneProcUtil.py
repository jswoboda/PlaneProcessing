#!/usr/bin/env python
"""
Created on Wed Dec 30 16:20:38 2015

@author: John Swoboda
"""
import os,glob,shutil
import numpy as np
import scipy as sp
import scipy.io as io
from RadarDataSim.IonoContainer import IonoContainer

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


def convertMattsfiles(fname_list,angle,keepspec=sp.array([0,1,2,6]),offset=0.,outdir=os.getcwd()):
    """ 
    This function will convert a set of files from Matt Zettegrens simulator to formated h5 files that 
    RadarDataSim can read. 
    Inputs
        fname_list - This is a list of .mat files that will be converted.
        angle - The angle in degrees that the thing will lie on.
        keepspec - A numpy array of numbers.
        offset - Offset of the x-axis.
        outdir - The directory this is getting saved.
    """
    angle_r = angle*sp.pi/180.
    if isinstance(fname_list,str):
        fname_list=[fname_list]

    for ifile in fname_list:
        dirname,fname = os.path.split(ifile)
        origfname,ext = os.path.splitext(fname)
        filetime= origfname.split('_')[-1]
        outfile = os.path.join(outdir,filetime+'.h5')

        matdict = io.loadmat(ifile)
        
        x1v=matdict['xg']['xp'][0][0]
        x3v=matdict['xg']['zp'][0][0]

        x1v=x1v+offset*1e3
        x1mat,x3mat = sp.meshgrid(x1v,x3v)


        E = x1mat.flatten()*sp.cos(angle_r)
        N = x1mat.flatten()*sp.sin(angle_r)
        U = x3mat.flatten()
        cart_coords = sp.column_stack((E,N,U))*1e-3
        lxs=x3mat.size
        
        t1 = matdict['t'][0][0]
        Time_Vector = sp.array([[t1,t1+15]])
        
        ns= matdict['ns']
        lsp = ns.shape[-1]
        ns = ns.reshape(lxs,lsp)
        # set up the tempretures
        Ts = matdict['Ts']
        Ts = Ts.reshape(lxs,lsp)
        # Set up the velocities
        vs= matdict['vsx1']
        vs = vs.reshape(lxs,lsp)
        vion=sp.sum(vs[:,:lxs-1]*ns[:,:lxs-1],1)/ns[:,-1]

        # reduce the number of species
        Ts = Ts[:,keepspec]
        Ts = Ts[:,sp.newaxis,:]
        ns = ns[:,keepspec]
        ns = ns[:,sp.newaxis,:]

        Param_List=sp.concatenate((ns[:,:,:,sp.newaxis],Ts[:,:,:,sp.newaxis]),3)
        v_e = vion*sp.cos(angle_r)
        v_n = vion*sp.sin(angle_r)
        v_u = sp.zeros_like(v_n)
        Velocity = sp.column_stack((v_e,v_n,v_u))
        Velocity = Velocity[:,sp.newaxis,:]
        Species = ['O+','NO+','N2+','O2+','N+', 'H+','e-']
        Species = [Species[i] for i in keepspec]
        Iono1 = IonoContainer(cart_coords,Param_List,times=Time_Vector,ver=0,species=Species,velocity=Velocity)
        
        Iono1.saveh5(outfile)
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
