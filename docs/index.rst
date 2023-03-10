.. turbospectrum documentation master file, created by
   sphinx-quickstart on Tue Feb 16 13:03:42 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

******
MOOGPy
******

Introduction
============
|moogpy| is a generical stellar spectral synthesis package that can be run from python.  It's essentially
a redistribution of the `MOOG <https://github.com/jsobeck/MOOG-SCAT_basic_git>`_ and `moog17scat <https://github.com/alexji/moog17scat>`_
Fortran spectral synthesis code by Chris Sneden with scattering updates from Jennifer Sobeck and other updates by Alex Ji.
In addition, there is Python wrapper/driver software where I reused some old IDL code I wrote and some code from Jon Holtzman
from the `APOGEE package <https://github.com/sdss/apogee>`_.  The setup.py file has also been modified to
automatically compile the Fortran code and copy them to the user's python scripts directory.  This reused code written by
Andy Casey in his `moog package <https://github.com/andycasey/moog>`_.

.. toctree::
   :maxdepth: 1

   install
   modules
	      

Description
===========
To run |moogpy| you need 1) a model atmosphere, 2) a linelist (or multiple), and 3) the set of stellar parameters
and elemental abundances that you want to run.

1) Model Atmospheres

   MOOG can read Kurucz/ATLAS or MARCS model atmospheres.  See pages 12-17 in the `MOOG manual <_static/WRITEnov2019.pdf>`_ for the format and examples.

2) Linelists

   MOOG requires a specific linelist format.  See pages 20-21 in the `MOOG manual <_static/WRITEnov2019.pdf>`_ for the format.
   
3) Stellar parameters and elemental abundances.

   The main stellar parameters are Teff, logg, [M/H], and [alpha/M].  These are the first four parameters in the
   main ``synthesis.synthesize()`` function.  The individual elements abundances are given in the ``elems`` parameters
   as a list of [element name, abundance] pairs, where abundance in the in [X/M] format relative to the overall metallicity.

More details on how the Fortran MOOG code works can be found in the `WRITEnov2019.pdf <_static/WRITEnov2019.pdf>`_ documentation.
   

Examples
========

.. toctree::
    :maxdepth: 1

    examples

*****
Index
*****

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
