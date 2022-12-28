#!/usr/bin/env python

import os
import time
import numpy as np

def moog_loadsum,sumfile,str,element_str,atmos,solar=abund_solar,silent=silent,                 verbose=verbose,stp=stp 
     
    #+ 
    # 
    # MOOG_LOADSUM 
    # 
    # This loads the MOOG "sumout" file 
    # 
    # INPUTS: 
    #  sumfile    The MOOG "sumout" filename 
    #  /silent    No output to the screen. 
    #  /verbose   Output all of the MOOG "sumout" information. 
    #  /stp       Stop at the end of the program. 
    # 
    # OUTPUTS: 
    #  str        IDL structure of the summary file.  For each line there 
    #               are: WAVE, SPECIES, EP, LOGGF, EW, LOGRW, ABUND, 
    #               DELAVG, ELEMENT. 
    #  element_str  A structure for each unique element: ELEMENT, ATOMNUM, 
    #                 NLINES, ABUND, M_H, M_FE. 
    #  atmos      IDL structure with the atmosphere parameters (Teff,logg,microtrub,metal) 
    #               used to run MOOG. 
    #  =solar     The MOOG solar abundance structure. 
    # 
    # USAGE: 
    #  IDL>moog_loadsum,'sumout.txt',str,element_str,atmos 
    # 
    # By D.Nidever  Sep 2009 
    #- 
     
    nsumfile = len(sumfile) 
     
    # Not enough inputs 
    if nsumfile == 0: 
        print('Syntax - moog_loadsum,sumfile,str,element_str,atmos,solar=solar,silent=silent,' 
        print('                      verbose=verbose,stp=stp' 
        return 
     
    # Check the sumfile file 
    if os.path.exists(sumfile) == 0: 
        print(sumfile,' NOT FOUND' 
        return 
     
     
    # Read the file 
    READLINE,sumfile,lines 
     
    # Loop through the lines 
    nlines = len(lines) 
    species = '' 
    undefine,str,str1 
    for i in range(nlines): 
         
        iline = lines[i] 
         
        # Beginning of new species 
        if strmid(iline,0,17) == 'Abundance Results': 
            species = strmid(iline,29,10) 
            species = strcompress(species,/remove_all) 
 
         
        arr = strsplit(iline,' ',/extract) 
        isnum = valid_num(arr[0]) 
         
        # Header line 
        if str(arr[0],2) == 'Teff': 
            equalind1 = strpos(iline,'=',0) 
            teff = float(strmid(iline,equalind1+1,6)) 
            equalind2 = strpos(iline,'=',equalind1+1) 
            logg = float(strmid(iline,equalind2+1,6)) 
            equalind3 = strpos(iline,'=',equalind2+1) 
            microturb = float(strmid(iline,equalind3+1,6)) 
            equalind4 = strpos(iline,'=',equalind3+1) 
            metal = float(strmid(iline,equalind4+1,6)) 
            atmos = {teff:teff,logg:logg,microturb:microturb,metal:metal} 
         
        if isnum == 1: 
            ncol = len(arr) 
            #fieldnames = ['WAVE','EP','LOGGF','EW','LOGRW','ABUND','DELAVG'] 
            #str = arr2str(reform(arr,ncol,1),fieldnames=fieldnames,/silent) 
            arr2 = float(arr) 
            str1 = {WAVE:arr2[0],SPECIES:species,EP:arr2[1],LOGGF:arr2[2],EW:arr2[3],LOGRW:arr2[4],           ABUND:arr2[5],DELAVG:arr2[6]} 
            PUSH,str,str1 
 
         
 
     
    nstr = len(str) 
    if not keyword_set(silent) : 
        print('N lines = ',str(nstr,2) 
     
    # Verbose output 
    if keyword_set(verbose): 
        print('MOOG OUTPUT' 
        print('' 
        printline,lines 
        print('' 
     
    # Atmospheric parameters 
    if not keyword_set(silent): 
        print('Atmospheric Parameters:' 
        print('  Teff=',stringize(atmos.teff,ndec=0),' logg=',stringize(atmos.logg,ndec=2),        ' vt=',stringize(atmos.microturb,ndec=2),' [M/H]=',stringize(atmos.metal,ndec=2) 
     
    # Get element names 
    two = strmid(str.species,0,2) 
    first = strmid(str.species,0,1) 
    second = strmid(str.species,1,1) 
    element = two 
    bd , = np.where(second == 'I',nbd) 
    if nbd > 0 : 
        element[bd] = first[bd] 
     
    ui = np.uniq(element,np.argsort(element)) 
    uelement = element[ui] 
    nel = len(ui) 
     
    # Add element names 
    add_tag,str,'ELEMENT','',str 
    str.element = element 
     
    if not keyword_set(silent) : 
        print('Nelements = ',str(nel,2) 
     
     
    # MOOG abundance and element information 
    # 
    # MOOG gives the abundances in log eps(X) = log( N(X)/N(H) ) + 12.0 
    #  It also uses the Anders and Grevesse (1989) solar abundances 
    #  and a meteoritic Fe abundances of 7.52. 
    # 
    # metal deficiencies: [M/H] = log(N(X)/N(H)) - log(N(X)/N(H))_solar 
    #  or [X/H] = log eps(x) - log eps(x)_solar 
    #  e.g. [Fe/H] = log eps(Fe) - 7.52 
     
    # From MOOG src/Batom/f file 
    elnames = ['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',           'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',           'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',           'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',           'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',           'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',           'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',           'Lu','Hf','Ta','Wl','Re','Os','Ir','Pt','Au','Hg',           'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',           'Pa','U ','Np','Pu','Am'] 
    elnames = str(elnames,2) 
    #c 
    #c  xabu = the set of current solar (when available) or meteorite 
    #c  abundances, scaled to log(h) = 12.00 .  The data are from Anders 
    #c  and Grevesse (1989, Geochim.Cosmichim.Acta, v53, p197) and the solar 
    #c  values are adopted except for a) uncertain solar data, or b) Li, Be, 
    #c  and B, for which the meteoritic values are adopted. 
    #c  I was told to use the new Fe value of 7.52 as adopted in Sneden 
    #c  et al. 1992 AJ 102 2001. 
    #      data xabu/ 
    # 
    elabund = [12.00,10.99, 3.31, 1.42, 2.88, 8.56, 8.05, 8.93, 4.56, 8.09,           6.33, 7.58, 6.47, 7.55, 5.45, 7.21, 5.5 , 6.56, 5.12, 6.36,           3.10, 4.99, 4.00, 5.67, 5.39, 7.52, 4.92, 6.25, 4.21, 4.60,           2.88, 3.41, 2.37, 3.35, 2.63, 3.23, 2.60, 2.90, 2.24, 2.60,           1.42, 1.92, 0.00, 1.84, 1.12, 1.69, 1.24, 1.86, 0.82, 2.0,           1.04, 2.24, 1.51, 2.23, 1.12, 2.13, 1.22, 1.55, 0.71, 1.50,           0.00, 1.00, 0.51, 1.12, 0.33, 1.1 , 0.50, 0.93, 0.13, 1.08,           0.12, 0.88, 0.13, 0.68, 0.27, 1.45, 1.35, 1.8 , 0.83, 1.09,           0.82, 1.85, 0.71, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.12,           0.00, 0.00, 0.00, 0.00, 0.00] 
     
    abund_solar = replicate({atomnum:0L,element:'',abund:0.0},len(elnames)) 
    abund_solar.atomnum = lindgen(len(elnames))+1 
    abund_solar.element = elnames 
    abund_solar.abund = elabund 
     
     
    # Add ATOMNUM to structure 
    ADD_TAG,str,'atomnum',-1L,str 
    for i in range(nstr): 
        ind , = np.where(abund_solar.element == str[i].element,nind) 
        if nind > 0 : 
            str[i].atomnum=abund_solar[ind[0]].atomnum 
 
     
     
    # Getting Fe abundance 
    ind , = np.where(str.element == 'Fe',nind) 
    if nind > 0: 
        abundfe = np.median(str[ind].abund) 
    else: 
        print('NO Fe abundance - using solar' 
        abundfe = 7.52 
    Fe_H = abundfe - 7.52# [Fe/H] 
     
    # Start the "element" structure 
    element_str = replicate({atomnum:0L,element:'',nlines:0L,abund:0.0,M_H:0.0,M_Fe:0.0},nel) 
     
    # Loop through the elements 
    for i in range(nel): 
        ind , = np.where(str.element == uelement[i],nind) 
        str1 = str[ind] 
        if nind > 1: 
            medabund = np.median(str1.abund) 
        else: 
            medabund=str1[0].abund 
         
        # Getting solar abundance 
        indsol , = np.where(abund_solar.element == uelement[i],nindsol) 
        atomnum = abund_solar[indsol[0]].atomnum 
        abund_sol = abund_solar[indsol[0]].abund 
         
        # Metal deficiency 
        #  metal deficiencies: [M/H] = log(N(X)/N(H)) - log(N(X)/N(H))_solar 
        #    or [X/H] = log eps(X) - log eps(X)_solar 
        #    e.g. [Fe/H] = log eps(Fe) - 7.52 
        M_H = medabund - abund_sol 
         
        # Metal deficiency versus Fe 
        # [X/H] = [X/Fe] + [Fe/H] and therefore 
        # [X/Fe] = [X/H] - [Fe/H] 
        M_Fe = M_H - Fe_H 
         
        # Put it in the structure 
        element_str[i].atomnum = atomnum 
        element_str[i].element = uelement[i] 
        element_str[i].nlines = nind 
        element_str[i].abund = medabund 
        element_str[i].M_H = M_H 
        element_str[i].M_Fe = M_Fe 
         
        #form='(A-4,I-4,I-4,3F8.2)' 
        #print,format=form,uelement[i],atomnum,nind,medabund,M_H,M_Fe 
        #form='(A-4,A-14,A-15)' 
        #print,format=form,uelement[i],'Abund = '+strtrim(string(medabund,format='(F9.2)'),2),  ;     'Nlines = '+strtrim(nind,2) 
 
     
    # Sort the element array by atomnum 
    si = np.argsort(element_str.atomnum) 
    element_str = element_str[si] 
     
     
    # Print out the Element information 
    if not keyword_set(silent): 
        print('-------------------------------------' 
        print('Elem A.Num Nline Abund  [M/H] [M/Fe]' 
        print('-------------------------------------' 
        for i in range(nel): 
            istr = element_str[i] 
            form='(A3,I5,I5,F8.2,2F7.2)' 
            print(format=form,istr.element,istr.atomnum,istr.nlines,istr.abund,istr.M_H,istr.M_Fe 
 
        print('-------------------------------------' 
     
     
    # Delete the MOOG sumout file 
    #if keyword_set(delete) then begin 
    #  if not keyword_set(silent) then print,'Deleting MOOG "sumout" file' 
    #  FILE_DELETE,sumfile,/allow 
    #endif 
     
     
    if keyword_set(stp) : 
        import pdb; pdb.set_trace() 
     
 
