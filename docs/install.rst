************
Installation
************


Important packages
==================
`moogpy` is a package to run the `MOOG <https://github.com/jsobeck/MOOG-SCAT_basic_git>`_ and `moog17scat <https://github.com/alexji/moog17scat>`_
spectral synthesis code by Chris Sneden and scattering improvements by Jennifer Sobeck and other updates by Alex Ji.
There is also a Python wrapper/driver based on some old IDL software I wrote and some code from Jon Holtzman in the
in the `APOGEE package <https://github.com/sdss/apogee>`_).

Installing MOOGPy
=================

The easiest way to install the code is with pip.  This will both compile and install the Fortran code as
well as the Python wrapper code.

.. code-block:: bash

    pip install moogpy

Fortran code
------------
    
The pip install will attempt to automatically compile the Fortran code and copy the binaries to your
Python scripts directory (which should be in your path).  If this fails for some reason, then you'll
need to compile it yourself.  You'll likely want to do a full git clone of the repository for this.
To compile the code you need either the ``g77` or the GNU Fortran compiler (``gfortran``).
The Fortran code lives in the `src/` directory.  All you should need to do is to cd into that
directory and type ``make -f Makefile.xxx`` where ``xxx`` is ``mac``, ``rh`` or ``rh64`` for linux/unix.
Copy the binaries ``MOOG`` and ``MOOGSILENT`` to a directory in your path (e.g., ~/bin/ or /usr/local/bin/).  

Dependencies
============

- numpy
- scipy
- astropy
- matplotlib
- `dlnpyutils <https://github.com/dnidever/dlnpyutils>`_
