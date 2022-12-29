import os
import numpy as np
import shutil
import subprocess
import tempfile
import time
from . import utils, atomic, atmos


def synthesize(teff,logg,mh=0.0,am=0.0,cm=0.0,nm=0.0,vmicro=2.0,elems=None,
               wrange=[15000.0,17000.0],dw=0.1,atmod=None,atmos_type='kurucz',
               dospherical=True,linelist=None,solarisotopes=False,workdir=None,
               save=False,verbose=False):
    """
    Code to synthesize a spectrum with MOOG.
    
    Parameters
    ----------
    teff : float
       Effective temperature in K.
    logg : float
       Surface gravity.
    mh : float, optional
       Metallicity, [M/H].  Deftauls is 0.0 (solar).
    am : float, optional
       Alpha abundance, [alpha/M].  Default is 0.0 (solar).
    cm : float, optional
       Carbon abundance, [C/M].  Default is 0.0 (solar).
    nm : float, optional
       Nitrogen abundance, [N/M].  Default is 0.0 (solar).
    vmicro : float, optional
       Microturbulence in km/s.  Default is 2 km/s.
    solarisotopes : bool, optional
       Use solar isotope ratios, else "giant" isotope ratios ( default False ).
    elems : list, optional
       List of [element name, abundance] pairs.
    wrange : list, optional
       Two element wavelength range in A.  Default is [15000.0,17000.0].
    dw : float, optional
       Wavelength step.  Default is 0.1 A.
    atmod : str, optional
       Name of atmosphere model (default=None, model is determined from input parameters).
    atmos_type : str, optional
       Type of model atmosphere file.  Default is 'kurucz'.
    dospherical : bool, optional
       Perform spherically-symmetric calculations (otherwise plane-parallel).  Default is True.
    linelist : str
       Linelist filename.
    save : bool, optional
       Save temporary directory and files for synthesis.  Default=False.
    workdir : str, optional
       Directory to perform the work in.  By default a temporary directory is
         created and deleted after the work is done (unless save=True).
    verbose : bool, optional
       Verbose output to the screen.

    Returns
    -------
    flux : numpy array
       The fluxed synthetic spectrum.
    continuum : numpy array
       The continuum of the spectrum.
    wave : numpy array
       Wavelength array in A.

    Example
    -------

    flux,cont,wave = synthesize(5000.0,2.5,-1.0)

    """

    t0 = time.time()
    
    # Default abundances
    abundances = atomic.solar()
    abundances[2:] += mh
    abundances[6-1] += cm
    abundances[7-1] += nm
    for i in [8,10,12,14,16,18,20,22]:
        abundances[i-1] += am
    # Abundance overrides from els, given as [X/M]
    if elems is not None :
        for el in elems:
            atomic_num = atomic.periodic(el[0])
            abundances[atomic_num-1] = atomic.solar(el[0]) + mh + el[1]
    # Set negative abundances to zero
    for i in range(len(abundances)):
        if abundances[i]<0:
            abundances[i] = 0.0
    # leave off the last one, 99 for moog means to scale ALL abundances by this value
    abundances = np.delete(abundances,98)
    
    # Change to temporary directory
    if workdir is None:
        workdir = tempfile.mkdtemp(prefix='moog')
    cwd = os.getcwd()
    os.chdir(workdir)

    # Create the root name from the input parameters
    root = (atmos_type+'_t{:04d}g{:s}m{:s}a{:s}c{:s}n{:s}v{:s}').format(int(teff), atmos.cval(logg), 
                      atmos.cval(mh), atmos.cval(am), atmos.cval(cm), atmos.cval(nm),atmos.cval(vmicro))

    # Check that linelists and model atmosphere files exit
    #if type(linelists) is str:
    #    linelists = [linelists]
    #for l in linelists:
    if os.path.exists(linelist)==False:
        raise FileNotFoundError(l)
    if os.path.exists(atmod)==False:
        raise FileNotFoundError(atmod)

    if dospherical and ('marcs' in atmos_type) and logg <= 3.001:
        spherical= True
    else:
        spherical = False
    flux,cont,wave = do_moog(root,atmod,linelist,mh,am,abundances,wrange,dw,
                             save=save,solarisotopes=solarisotopes)

    os.chdir(cwd)
    if not save:
        shutil.rmtree(workdir)
    else:
        print('Saving temporary directory '+workdir)
        
    if verbose:
        print('dt = {:.3f}s'.format(time.time()-t0))
        
    return flux,cont,wave

    
def do_moog(root,atmod,linefile,mh,am,abundances,wrange,dw,save=False,
                 solarisotopes=False,babsma=None,atmos_type='marcs',
                 spherical=True,vmicro=2.0,tfactor=1.0,verbose=False):
    """
    Runs MOOG for specified input parameters.

    Parameters
    ----------
    root : str
       Root of filenames to use for this MOOG run.
    atmod : str, optional
       Name of atmosphere model (default=None, model is determined from input parameters).
    linefile : str
       Linelist filename.
    mh : float, optional
       Metallicity, [M/H].  Default is 0.0 (solar).
    am : float, optional
       Alpha abundance, [alpha/M].  Default is 0.0 (solar).
    abundances : list
       List of abundances.
    wrange : list, optional
       Two element wavelength range in A.  Default is [15000.0,17000.0].
    dw : float, optional
       Wavelength step.  Default is 0.1 A.
    solarisotopes : bool, optional
       Use solar isotope ratios, else "giant" isotope ratios.  Default is False.
    babsma : bool, optional
       The name of hte babsma output file of opacities.  Default is None to run babmsa.
    atmos_type : str, optional
       Model atmosphere type.  Default is 'marcs'.
    spherical : bool, optional
       Spherical atmosphere.  Default is True.
    vmicro : float, optional
       Microturbulent velocity in km/s.  Default is 2.0 km/s.
    save : bool, optional
       Save temporary directory and files for synthesis.  Default=False.
    verbose : bool, optional
       Verbose output to the screen.

    Returns
    -------
    flux : numpy array
       The fluxed synthetic spectrum.
    continuum : numpy array
       The continuum of the spectrum.
    wave : numpy array
       Wavelength array in A.

    Example
    -------

    flux,cont,wave = do_moog(root,atmod,linefile,-0.1,0.2,abund,wrange=[15000.0,17000.0],dw=0.1)

    """

    # MOOG setup
    #datadir = utils.datadir()
    #os.symlink(datadir,'./DATA')
    shutil.copy(atmod,'./'+os.path.basename(atmod))
    atmosfile = os.path.basename(atmod)
    alines = utils.readlines(atmosfile)
    
    # Prepare the end of the MOOG model atmosphere file 
    #NATOMS        3    0.10 
    #       6.0      7.70       8.0      8.30       7.0      7.50 
    #NMOL         21 
    #     606.0     106.0     607.0     608.0     107.0     108.0     112.0 
    #707.0 
    #     708.0     808.0      12.1   60808.0   10108.0     101.0   60606.0 
    #839.0 
    #     840.0     822.0      22.1      40.1      39.1 

    # MOOG internally stores the Asplund+2009 solar abundances
        
    # Create MOOG control file
    addon = []
    natoms = len(abundances)
    addon.append('NATOMS        {0:d}    {1:.3f}'.format(natoms,mh))
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
    #logeps_mg = alpha + metal + 7.58 
    #logeps_si = alpha + metal + 7.55 
    #logeps_s  = alpha + metal + 7.21 
    #logeps_ar = alpha + metal + 6.56 
    #logeps_ca = alpha + metal + 6.36 
    #logeps_ti = alpha + metal + 4.99
    #addon.append('      {0:4.1f}      {1:.3f}'.format(12.0,logeps_mg))  # Mg
    #addon.append('      {0:4.1f}      {1:.3f}'.format(14.0,logeps_si))  # Si 
    #addon.append('      {0:4.1f}      {1:.3f}'.format(16.0,logeps_s))   # S 
    #addon.append('      {0:4.1f}      {1:.3f}'.format(18.0,logeps_ar))  # Ar 
    #addon.append('      {0:4.1f}      {1:.3f}'.format(20.0,logeps_ca))  # Ca 
    #addon.append('      {0:4.1f}      {1:.3f}'.format(22.0,logeps_ti))  # Ti
    for iel,abun in enumerate(abundances):
        addon.append("      {0:4.1f}      {1:8.3f}".format(iel+1,abun))
    # NO MOLECULES FOR NOW!!!
    newalines = alines + addon
    utils.writelines(atmosfile,newalines,overwrite=True)
    
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
    w0 = wrange[0]
    w1 = wrange[1]
    fwhm = 0.01  # Gaussian broadening 
     
    # Read in the linelist
    linelist = utils.readlines(linefile,comment='#')
    nlinelist = len(linelist)
    lwave = np.array([float(l.split()[0]) for l in linelist]).astype(float)
    #READLINE,linefile,linelist,count=nlinelist,comment='#' 
    #dum = strsplitter(linelist,' ',/extract) 
    #lwave = float(reform(dum[0,:])) 
    wavemin = np.min(lwave) 
    wavemax = np.max(lwave) 
    
    # Make temporary linelist file 
    #templist = MKTEMP('line')
    tid,templist = tempfile.mkstemp(prefix='line')
    templist = os.path.basename(templist)
    gd, = np.where((lwave >= w0) & (lwave <= w1))
    tlinelist = np.char.array(linelist)[gd]
    utils.writelines(templist,tlinelist)
         
    # This is what the MOOG input file looks like 
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
    
    # Set up the parameter file 
    params = []
    params.append("synth")
    params.append("terminal       'xterm'")
    params.append("standard_out   out1")
    params.append("summary_out    out2")
    params.append("smoothed_out   out3")
    params.append("model_in       '{:s}'".format(atmosfile))
    params.append("lines_in       '{:s}'".format(templist))
    params.append("atmosphere    1")
    params.append("molecules     2")
    params.append("lines         1")
    params.append("flux/int      0")
    params.append("damping       0")
    #params.append("abundances    1   1" 
    #params.append("        3     -0.87" 
    #***!!!!! NOT REQUIRED IF INPUT IN MODEL ATMOSPHERE!!!!! *** 
    # ADD ABUNDANCES OF ALPHA ELEMENTS 
    # For spectrum use these alpha elements: Mg, Si, S, Ar, Ca and Ti. 
    # Mg-12, Si-14, S-16, Ar-18, Ca-20, Ti-22 
    #if (float(alpha) ne 0.0) then begin 
    #  PUSH,tpaarms,'abundances    6   1' 
    #  params.append('       12     '+strtrim(string(logeps_mg,format='(F8.2)'),2) ; Mg 
    #  params.append('       14     '+strtrim(string(logeps_si,format='(F8.2)'),2) ; Si 
    #  params.append('       16     '+strtrim(string(logeps_s,format='(F8.2)'),2) ; S 
    #  params.append('       18     '+strtrim(string(logeps_ar,format='(F8.2)'),2) ; Ar 
    #  params.append('       20     '+strtrim(string(logeps_ca,format='(F8.2)'),2) ; Ca 
    #  params.append('       22     '+strtrim(string(logeps_ti,format='(F8.2)'),2) ; Ti 
    #end 
    params.append("synlimits")
    params.append("  {0:10.3f}  {1:10.3f}  {2:10.3f}  {3:10.3f}".format(w0,w1,dw,fwhm))
    params.append("plotpars      1")
    params.append("  {0:10.3f}  {1:10.3f}   0.0   1.00".format(w0,w1))
    params.append("   0.0       0.0    0.00   1.0")
    params.append("    g     {:.3f}      0.0    0.00     0.0     0.0   ".format(fwhm))
    params.append("opacit        0")
    # MOOGSILENT is hardwired to use "batch.par"
    utils.writelines('batch.par',params)
    
    # Run MOOGSILENT
    ret = subprocess.check_output(['MOOGSILENT'],stderr=subprocess.STDOUT)    
    # Save the log file
    if type(ret) is bytes:
        ret = ret.decode()
    with open(root+'_MOOG.log','w') as f:
        f.write(ret)

    # Load the spectrum
    wave,flux = utils.read_synthfile('out2')
    cont = np.zeros(flux.shape,float)
    
    return flux,cont,wave
