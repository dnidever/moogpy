#!/usr/bin/env python

import os
import time
import numpy as np

def sort_mooglinelist,file1,file2,overwrite=overwrite 
     
    #+ 
    # 
    # SORT_MOOGLINELIST 
    # 
    # This program reads in a MOOG-formatted linelist and 
    # sorts it by species. 
    # 
    # INPUTS: 
    #  file1       Filename of input MOOG linelist 
    #  file2       Filename of output sorted MOOG linelist to create 
    #                The default is "sortlines.txt" 
    #  /overwrite  Overwrite the output file if it exists already. 
    # 
    # OUTPUTS: 
    #  A species-sorted MOOG linelist is created 
    # 
    # USAGE: 
    #  IDL>sort_mooglinelist,'lines1.txt','lines2.txt' 
    # 
    # By D.Nidever  Sep.2009 
    #- 
     
    # Not enough inputs 
    if n_params() == 0: 
        print('Syntax - sort_mooglinelist,file1,file2,overwrite=overwrite' 
        return 
     
    if os.path.exists(file1) == 0: 
        print(file1,' NOT FOUND' 
        return 
     
    print('Reading file >>',file1,'<<' 
     
    readline,file1,lines 
    arr = strsplitter(lines,' ',/extract) 
    lines0 = lines 
     
    # First line is header 
    header = 0 
    if valid_num(arr[0,0]) == 0: 
        print('Header line' 
        header = 1 
        head = lines[0] 
        arr = arr[:,1::] 
        lines = lines[1::] 
 
     
     
     
    # Sort by lines 
    species = float(arr[1,:]) 
    si = np.argsort(species) 
    ui = np.uniq(species,np.argsort(species)) 
    ui2 = np.uniq(int(species),np.argsort(int(species))) 
    print('N lines = ',str(len(lines),2) 
    print('N elements = ',str(len(ui2),2) 
    print('N species = ',str(len(ui),2) 
     
    # Number of lines per species 
    print('' 
    print('Species Nlines' 
    uspecies = species[ui] 
    for i in range(len(ui)): 
        gd , = np.where(species == uspecies[i],ngd) 
        print(uspecies[i],ngd,format='(F6.1,I5)' 
 
    print('' 
     
    # Making output list 
    lines2 = lines[si] 
    species2 = species[si] 
     
    # Put iron files firt 
    gd , = np.where(int(species2) == 26,ngd) 
    if ngd > 0: 
        ironlines = lines2[gd] 
        otherlines = lines2 
        REMOVE,gd,otherlines 
        lines2 = [ironlines,otherlines] 
     
     
     
    # Write output file 
    if len(file2) == 0: 
        outfile='sortlines.txt' 
    else: 
        outfile=file2 
    if header == 1 : 
        lines2=[head,lines2] 
    if os.path.exists(outfile) == 1 and not keyword_set(overwrite): 
        print(outfile,' ALREADY EXISTS' 
        return 
    print('Writing sorted file to >>',outfile,'<<' 
    WRITELINE,outfile,lines2 
     
    #stop 
     
 
