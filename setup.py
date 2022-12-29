#from distutils.core import setup
#from distutils.sysconfig import get_python_lib

import os
import sys
import fileinput
import platform
import subprocess
from setuptools.command.install import install
from setuptools import setup, find_packages

from platform import system as current_platform
from shutil import copy, move, copytree, copyfile
from glob import glob

def get_bin_path():
    # Get environment bin/ directory
    bindir = None
    # check if the --user option was set
    userdir = None
    for u in install.user_options:
        if u[0]=='user':
            # "install in user site-package '/Users/nidever/.local/lib/python3.7/site-packages'"
            uline = u[2]
            userdir = uline[uline.find('site-package')+13:]
            userdir = userdir.replace("'","")
    if userdir is not None:
        # /Users/nidever/.local/bin
        bindir = os.path.dirname(os.path.dirname(os.path.dirname(userdir)))+'/bin'
        if os.path.exists(bindir)==False:
            bindir = False
    # Try virtual environment using sys.prefix
    if bindir is None:
        venv = get_virtualenv_path()
        if venv is not None:
            bindir = venv+'/bin'
            if os.path.exists(bindir)==False:
                bindir = None
    # Get bin/ directory from python executable
    if bindir is None:
        out = subprocess.run(['which','python'],shell=True)
        if type(out) is bytes:
            out = out.decode()
        bindir = os.path.dirname(out)
        if os.path.exists(bindir)==False:
            bindir = None

    if bindir is None:
        raise Exception('No bin/ directory found')

    return bindir


def compilemoog():
    """ Compile the MOOG Fortran code."""
    
    # Identify the platform
    platform = current_platform()

    # Check for platform first
    if platform not in ('Darwin', 'Linux'):
        sys.stderr.write("Platform '%s' not recognised!\n" % platform)
        sys.exit()

    # By default, we will use 32bit 
    is_64bits = False

    # Which system are we on?
    if platform == 'Darwin':
        run_make_files = ('Makefile.mac', 'Makefile.macsilent')
        #run_make_files = ('Makefile.macsynth',)
        machine = 'mac'
    elif platform == 'Linux':
        machine = 'pcl'
        is_64bits = sys.maxsize > 2**32
        if is_64bits:
            run_make_files = ('Makefile.rh64', 'Makefile.rh64silent')
            #run_make_files = ('Makefile.rh64',) 
        else:
            run_make_files = ('Makefile.rh', 'Makefile.rhsilent')
            #run_make_files = ('Makefile.rh',)

    # Check for gfortran or g77
    def system_call(command):
        """ Perform a system call with a subprocess pipe """
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        return process.communicate()[0]

    # Look for g77 and gfortran
    g77_exists = len(system_call("which g77")) > 0
    gfortran_exists = len(system_call("which gfortran")) > 0

    # If we have the choice, use gfortran
    if gfortran_exists:
        if is_64bits:
            fortran_vars = "FC = gfortran -m64\nFFLAGS = -Wall -O4 -ffixed-line-length-72 -ff2c"
        else:
            fortran_vars = "FC = gfortran\nFFLAGS = -Wall -O4 -ffixed-line-length-72 -ff2c"
    elif g77_exists:
        if platform == 'Linux':
            fortran_vars = 'FC = g77 -Wall'
        else:
            fortran_vars = 'FC = g77 -w'
    else:
        sys.stderr.write("Could not find g77 or gfortran on the system!\n")
        #sys.exit()
        return
        
    # Get our directories relative to the current path
    repository_dir = os.path.dirname(os.path.realpath(__file__))
    
    # We need a moog data directory
    data_dir = os.path.expanduser('~/.moog')
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

    # Copy files to data directory
    src_dir = os.path.join(repository_dir, 'src')
    data_files = glob('%s/*.dat' % src_dir)
    [copy(data_file, '%s/%s' % (data_dir, os.path.basename(data_file), )) for data_file in data_files]
    
    aqlib = "AQLIB = %s" % os.path.join(repository_dir, 'lib/aqlib')
    smlib = "SMLIB = %s" % os.path.join(repository_dir, 'lib/smlib')

    configuration = "\n".join([fortran_vars, aqlib, smlib])

    # Update the makefiles with the proper SMLIB and AQLIB
    run_make_files = [os.path.join(repository_dir, 'src', filename) for filename in run_make_files]
    #hardcoded_moog_files = [os.path.join(repository_dir, 'src', filename) for filename in ('Begin.f', 'Moog.f', 'Moogsilent.f')]
    hardcoded_moog_files = [os.path.join(repository_dir, 'src', filename) for filename in ('Begin.f', 'Moogsynth.f')]    
    
    # Setup: Move and create copies of the original
    for make_file in run_make_files:
        move(make_file, make_file + '.original')
        copy(make_file + '.original', make_file)

    for moog_file in hardcoded_moog_files:
        move(moog_file, moog_file + '.original')
        copy(moog_file + '.original', moog_file)

    # Update the run make files with the configuration
    for line in fileinput.input(run_make_files, inplace=True):
        line = line.replace('#$CONFIGURATION', configuration)

        sys.stdout.write(line)
        
    # Update the MOOG files
    for line in fileinput.input(hardcoded_moog_files, inplace=True):
        line = line.replace('$SRCDIR', src_dir)
        line = line.replace('$DATADIR', data_dir)
        line = line.replace('$MACHINE', machine)

        sys.stdout.write(line)

    # Run the appropriate make files
    print('Compiling the MOOG Fortran code')
    for make_file in run_make_files:
        os.system('cd src;make -f %s' % make_file)

    # Cleanup files: Replace with original files
    [move(moog_file + '.original', moog_file) for moog_file in hardcoded_moog_files if os.path.exists(moog_file + '.original')]
    [move(make_file + '.original', make_file) for make_file in run_make_files if os.path.exists(make_file + '.original')]

    # Copy to bin/ directory
    if os.path.exists('bin/')==False:
        os.mkdir('bin/')
    for f in ['MOOG','MOOGSILENT']:
    #for f in ['MOOGSYNTH']:        
        if os.path.exists('src/'+f):
            if os.path.exists('bin/'+f): os.remove('bin/'+f)
            copy('src/'+f,'bin/'+f)
        else:
            print('ERROR: '+f+' NOT found')
    
    # Get the path for the binaries
    bindir = get_bin_path()
        
    # Copy fortran binaries to bin/ directory
    for f in ['MOOG','MOOGSILENT']:
    #for f in ['MOOGSYNTH']:                
        if os.path.exists('bin/'+f):
            if os.path.exists(bindir+f): os.remove(bindir+f)
            print('Copying bin/'+f+' -> '+bindir+'/'+f)
            copyfile('bin/'+f,bindir+'/'+f)
            # Make executable
            os.chmod(bindir+'/'+f,0o755)
        else:
            print('bin/'+f+' NOT FOUND')

    
    ## Copy the AquaTerm framework
    #if not os.path.exists('/Library/Frameworks/AquaTerm.framework/'):
    #    try:
    #        system_call('cp -R %s /Library/Frameworks/AquaTerm.framework/' % os.path.join(repository_dir, 'lib/AquaTerm.framework/'))
    #    except:
    #        sys.stdout.write("AquaTerm framework could not be installed to /Library/Frameworks/AquaTerm.framework\n")
    #    else:
    #        sys.stdout.write("AquaTerm framework copied to /Library/Frameworks/AquaTerm.framework\n")

# We need to build MOOG and MOOGSILENT before they get moved to the scripts/
# directory so that they can be moved into the $PATH
if 'install' in sys.argv or 'develop' in sys.argv or 'bdist_wheel' in sys.argv:
    compilemoog()

setup(name='moogpy',
      version='1.0.8',
      description='MOOG spectral synthesis code and python wrapper',
      author='David Nidever',
      author_email='dnidever@montana.edu',
      url='https://github.com/dnidever/moogpy',
      #scripts=['bin/MOOG','bin/MOOGSILENT'],
      requires=['numpy','astropy(>=4.0)','scipy','matplotlib','dlnpyutils'],
      #cmdclass={'install': CustomInstall},
      #zip_safe = False,
      include_package_data=True,
      #packages=find_packages(where="python"),
      packages = ['moogpy'],
      #packages=find_namespace_packages(where="python"),
      package_dir={"": "python"}   
)
