#!/usr/bin/env python

import os
import time
import numpy as np

def moog_driver,temp,logg,metal,alpha,outfile 
     
    # This runs moog, multiple times on small chunks 
    # if necessary 
     
    ntemp = len(temp) 
    nlogg = len(logg) 
    nmetal = len(metal) 
    nalpha = len(alpha) 
     
    # Not enough inputs 
    if (ntemp == 0 or nlogg == 0 or nmetal == 0 or nalpha == 0): 
        print('moog_driver,temp,logg,metal,alpha,outfile' 
        return 
     
    #template = 'test.inp5' 
    linefile = 'vald_linelist.dat' 
    #model = 'atmos' 
    #temp = 5050. 
    #logg = 4.5 
    #metal = -1.5 
    #alpha = 0.3 
     
    print('Teff = ',str(temp,2) 
    print('Logg = ',str(logg,2) 
    print('[Fe/H] = ',str(metal,2) 
    print('[alpha/Fe] = ',str(alpha,2) 
     
    # Prepared the MODEL atmosphere file. 
    modeldir = '/net/halo/dln5q/kurucz/' 
    # Get the RIGHT alpha value 
    alpha_atmos = alpha 
    if alpha < 0.0 : 
        alpha_atmos = 0.0 
    if alpha > 0.4 : 
        alpha_atmos = 0.4 
    if metal < -2.5 : 
        alpha_atmos = 0.0 
    if alpha_atmos != alpha : 
        print('Using [alpha/Fe]=',str(alpha_atmos,2),' for the Kurucz model' 
    KURUCZ_MODEL,temp,logg,metal,alpha_atmos,filename,/moog 
     
    # THESE ARE WHAT KURUCZ CONSIDERED "ALPHA" ELEMENTS WHEN HE MADE 
    # THE MODEL ATMOSPHERES!!!!! 
    # For model atmospheres only the alpha elements are: O, Ne, Mg, Si, S, 
    # Ar, Ca, and Ti.  From Kirby. 
    # O-8, Ne-10, Mg-12, Si-14, S-16, Ar-18, Ca-20, Ti-22 
     
     
    # Copy the MOOG model atmosphere file to this directory 
    CD,current=curdir 
    FILE_COPY,modeldir+filename+'.moog',curdir,/allow,/over 
     
    # Prepare the end of the MOOG model atmosphere file 
    #NATOMS        3    0.10 
    #       6.0      7.70       8.0      8.30       7.0      7.50 
    #NMOL         21 
    #     606.0     106.0     607.0     608.0     107.0     108.0     112.0 
    #707.0 
    #     708.0     808.0      12.1   60808.0   10108.0     101.0   60606.0 
    #839.0 
    #     840.0     822.0      22.1      40.1      39.1 
    strmetal = str(string(metal,format='(F10.2)'),2) 
    if float(alpha) != 0.0: 
        natoms='6' 
    else: 
        natoms='0' 
    addon = 'NATOMS        '+natoms+'    '+strmetal 
    # ADD IN THE ALPHA ABUNDANCES 
    # For MOOG the "alpha" elements are: Mg, Si, S, Ar, Ca and Ti. 
    # Mg-12, Si-14, S-16, Ar-18, Ca-20, Ti-22 
    # according to Kirby 
     
    # This needs to be in log(epsilon) format: log(eps) = log(N(X)/N(H)) + 12.0 
    # [A/B] = log(n(A)/n(B))- log(n(A)/n(B))_solar 
    #  where n() is number density 
    # 12 + log[n(Fe)/n(H)]_sol = 7.52  ; solar [Fe/H] 
     
    # CONVERT TO LOG EPSILON FORMAT 
    # log eps(X) = log(N(X)/N(H)) + 12.0 
    # [X/H] = log(n(X)/n(H)) - log(n(X)/n(H))_solar 
    # log(n(X)/n(H)) = [X/H] + log(n(X)/n(H))_solar 
    # log eps(X) = [X/H] + log(n(X)/n(H))_solar + 12.0 
    # log eps(X) = [X/Fe] + [Fe/H] + (log(n(X)/n(H))_solar + 12.0) 
    #INPUT ABUNDANCES: (log10 number densities, log H=12) 
    #      Default solar abundances: Anders and Grevesse 1989 
    # This is the (log(n(X)/n(H))_solar + 12.0) value 
    # from Batom.f 
    # Mg(12)= 7.58 
    # Si(14)= 7.55 
    # S (16)= 7.21 
    # Ar(18)= 6.56 
    # Ca(20)= 6.36 
    # Ti(22)= 4.99 
    logeps_mg = alpha + metal + 7.58 
    logeps_si = alpha + metal + 7.55 
    logeps_s  = alpha + metal + 7.21 
    logeps_ar = alpha + metal + 6.56 
    logeps_ca = alpha + metal + 6.36 
    logeps_ti = alpha + metal + 4.99 
    PUSH,addon,'      12.0      '+str(string(logeps_mg,format='(F8.2)'),2)# Mg 
    PUSH,addon,'      14.0      '+str(string(logeps_si,format='(F8.2)'),2)# Si 
    PUSH,addon,'      16.0      '+str(string(logeps_s,format='(F8.2)'),2)# S 
    PUSH,addon,'      18.0      '+str(string(logeps_ar,format='(F8.2)'),2)# Ar 
    PUSH,addon,'      20.0      '+str(string(logeps_ca,format='(F8.2)'),2)# Ca 
    PUSH,addon,'      22.0      '+str(string(logeps_ti,format='(F8.2)'),2)# Ti 
    # NO MOLECULES FOR NOW!!! 
    WRITELINE,filename+'.moog',addon,/append 
    model = filename+'.moog' 
     
    # Mg, Si, S, Ar, Ca and Ti. 
    # Mg-12, Si-14, S-16, Ar-18, Ca-20, Ti-22 
    # [alpha/H] = [Fe/H] + [alpha/Fe] 
    # 
    # 12 + log[n(Fe)/n(H)]_sol = 7.52  ; solar [Fe/H] 
    # 
    # Microturbulence equation: vt = 2.700 - 0.509*logg 
     
     
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 
    # RUN MOOG 
     
    # Wavelength parameters 
    w0 = 3900.0#3805.0 
    w1 = 7000.0#7000.0 
    dw = 0.02 
    fwhm = 1.0# Gaussian broadening 
     
    step = 100.0#100.0  ; chunk size in A 
     
    # Read in the linelist 
    READLINE,linefile,linelist,count=nlinelist,comment='#' 
    dum = strsplitter(linelist,' ',/extract) 
    lwave = float(reform(dum[0,:])) 
    wavemin = min(lwave) 
    wavemax = max(lwave) 
     
    # Read in the parameter file 
    #READLINE,template,params 
     
    t0 = time.time() 
     
    # Loop through the chunks 
    nsteps = ceil((w1-w0)/step) 
    for i in range(nsteps): 
         
        # What wavelength range for this chunk 
        wstart = i*step+w0 
        wend = (wstart+step) < w1 
         
        # Add 2A overhang on both ends 
        wstart1 = (wstart-2.0) > wavemin 
        wend1 = (wend+2.0) < wavemax 
        #wstart1 = wstart 
        #wend1 = wend 
         
        print(i+1,wstart,wend,format='(I4,2F10.3)' 
         
        # Make temporary linelist file 
        templist = MKTEMP('line') 
        gd , = np.where(lwave >= wstart1 and lwave <= wend1,ngd) 
        tlinelist = linelist[gd] 
        WRITELINE,templist,tlinelist 
         
         
        # This is what the MOOGsynth input file looks like 
        #terminal       'xterm' 
        #model_in       'atmos' 
        #lines_in       'vald.dat' 
        #atmosphere    1 
        #molecules     2 
        #lines         1 
        #flux/int      0 
        #damping       0 
        #abundances    1   1 
        #        3     -0.87 
        #synlimits 
        #  5085.0  5320.0    0.02    1.00 
        #plotpars      1 
        #  5085.0  5320.0    0.00    1.00 
        #   0.0       0.0    0.00   1.0 
        #    g      1.0      0.0    0.00     0.0     0.0 
        #opacit        0 
         
        strwstart = str(string(wstart1,format='(F10.1)'),2) 
        strwend = str(string(wend1,format='(F10.1)'),2) 
        strdw = str(string(dw,format='(F10.3)'),2) 
        strfwhm = str(string(fwhm,format='(F10.2)'),2) 
         
        # Set up the parameter file 
        undefine,tparams 
        PUSH,tparams,"terminal       'xterm'" 
        PUSH,tparams,"model_in       '"+model+"'" 
        PUSH,tparams,"lines_in       '"+templist+"'" 
        PUSH,tparams,"atmosphere    1" 
        PUSH,tparams,"molecules     0"# NO MOLECULES FOR NOW 
        PUSH,tparams,"lines         1" 
        PUSH,tparams,"flux/int      0" 
        PUSH,tparams,"damping       0" 
        #PUSH,tparams,"abundances    1   1" 
        #PUSH,tparams,"        3     -0.87" 
        #***!!!!! NOT REQUIRED IF INPUT IN MODEL ATMOSPHERE!!!!! *** 
        # ADD ABUNDANCES OF ALPHA ELEMENTS 
        # For spectrum use these alpha elements: Mg, Si, S, Ar, Ca and Ti. 
        # Mg-12, Si-14, S-16, Ar-18, Ca-20, Ti-22 
        #if (float(alpha) ne 0.0) then begin 
        #  PUSH,tpaarms,'abundances    6   1' 
        #  PUSH,tparams,'       12     '+strtrim(string(logeps_mg,format='(F8.2)'),2) ; Mg 
        #  PUSH,tparams,'       14     '+strtrim(string(logeps_si,format='(F8.2)'),2) ; Si 
        #  PUSH,tparams,'       16     '+strtrim(string(logeps_s,format='(F8.2)'),2) ; S 
        #  PUSH,tparams,'       18     '+strtrim(string(logeps_ar,format='(F8.2)'),2) ; Ar 
        #  PUSH,tparams,'       20     '+strtrim(string(logeps_ca,format='(F8.2)'),2) ; Ca 
        #  PUSH,tparams,'       22     '+strtrim(string(logeps_ti,format='(F8.2)'),2) ; Ti 
        #end 
        PUSH,tparams,"synlimits" 
        PUSH,tparams,"  "+strwstart+"  "+strwend+"  "+strdw+"  1.00" 
        #PUSH,tparams,"  5085.0  5320.0    0.02    1.00" 
        PUSH,tparams,"plotpars      1" 
        PUSH,tparams,"  "+strwstart+"  "+strwend+"   0.0   1.00" 
        #PUSH,tparams,"  5085.0  5320.0    0.00    1.00" 
        PUSH,tparams,"   0.0       0.0    0.00   1.0" 
        PUSH,tparams,"    g     "+strfwhm+"      0.0    0.00     0.0     0.0   " 
        #PUSH,tparams,"    g      1.0      0.0    0.00     0.0     0.0   " 
        PUSH,tparams,"opacit        0" 
         
        # Make the temporary MOOGsynth input parameter file 
        tempfile = MKTEMP('moog') 
        WRITELINE,tempfile,tparams 
         
        # Run MOOGsynth 
        SPAWN,'./MOOGsynth '+tempfile,out,errout 
         
        # Load the spectrum 
        str = IMPORTASCII(tempfile+'.synout',fieldnames=['wave','flux'],fieldtypes=[4,4],skip=2,/noprint) 
        wave = str.wave 
        flux = str.flux 
        # Load the RAW spectrum 
        #LOADSUMOUT,tempfile+'.sumout',wave,flux,stp=stp 
         
        # Trim the edges 
        wave_orig = wave 
        flux_orig = flux 
        gg , = np.where(wave >= wstart and wave <= wend,ngg) 
        wave = wave[gg] 
        flux = flux[gg] 
         
        # Concatenate the spectra 
        if (i > 0): 
             
            # The first pixel will overlap, REMOVE it 
            REMOVE,0,wave,flux 
            # Concatenate with FINAL arrays 
            finalwave = [finalwave,wave] 
            finalflux = [finalflux,flux] 
             
        else: 
            finalwave = wave 
            finalflux = flux 
         
        plot,finalwave,finalflux,xr=[w0,w1],xs=1 
         
        # Gaussian convolution 
        #xx = scale_vector(findgen(1000),-100,100) 
        #gau = exp(-0.5*(xx^2.0)/(50.^2.0)^2.0) 
        #gau = gau/total(gau) 
        #flux2 = CONVOL(finaflux,gau,/center,/edge_truncate) 
         
        #stop 
         
        # Erase the temporary files 
        os.remove(tempfile,/allow,/quiet 
        os.remove(tempfile[0]+['.sumout','.stdout','.synout'],/allow,/quiet 
        os.remove(templist,/allow,/quiet 
         
        #stop 
         
 
     
    dt = time.time()-t0 
    print('Time = ',str(dt,2),' sec' 
     
    # Save the final spectrum 
    wave = finalwave 
    flux = finalflux 
    if not keyword_set(outfile) : 
        outfile = filename+'_spec.dat' 
    print('Saving final spectrum to ',outfile 
    save,wave,flux,file=outfile 
     
    import pdb; pdb.set_trace() 
     
 
