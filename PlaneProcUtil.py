#!/usr/bin/env python
"""
Created on Wed Dec 30 16:20:38 2015

@author: John Swoboda
"""
import os,glob,shutil
import numpy as np


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