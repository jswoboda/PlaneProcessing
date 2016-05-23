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
    if ~os.path.isdir(folder2):
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

        
