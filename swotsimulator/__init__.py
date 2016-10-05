# #!/usr/bin/env python
# =======================================================================
#                        General Documentation

"""Utilities for SWOT Science Simulator for the ocean

   Some useful online help commands for the package:
   * help(swotsimulator):  Help for the package.  A list of all modules in
     this package is found in the "Package Contents" section of the
     help output.
   * help(swotsimulator.M):  Details of each module "M", where "M" is the
     module's name.

#-----------------------------------------------------------------------
#                       Additional Documentation
# Authors: Lucile Gaultier and Clement Ubelmann
#
# Modification History:
# - Jul 2014:  Original by Clement Ubelmann and Lucile Gaultier, JPL
# - Nov 2014: Beta version
# - Feb 2015: Version 1.0
# - Dec 2015: Version 2.0
#
# Notes:
# - Written for Python 2.3,  Python 3.4, tested with Python 2.7, Python 3.4
#
# Copyright (c)
# Copyright (c) 2002-2014, California Institute of Technology.
# All rights reserved. Based on Government Sponsored Research
# under contracts NAS7-1407 and/or NAS7-03001.
#
#-----------------------------------------------------------------------
"""
# -----------------------------------------------------------------------

# ---------------- Module General Import and Declarations ---------------

# - Set module version to package version:

import swotsimulator.package_version
__version__ = package_version.version
__author__  = package_version.author
__date__    = package_version.date
__email__   = package_version.email
__url__     = package_version.url
del package_version

# - If you're importing this module in testing mode, or you're running
#  pydoc on this module via the command line, import user-specific
#  settings to make sure any non-standard libraries are found:

import os
import sys
if (__name__ == "__main__") or \
   ("pydoc" in os.path.basename(sys.argv[0])):
    import user


# - Find python version number
__python_version__ = sys.version[:3]

# - Import numerical array formats
try:
    import numpy
except:
    print(''' Numpy is not available on this platform,
          ''')

# - Import scientific librairies
try:
    import scipy
except:
    print("""Scipy is not available on this platform,
          """)


# - Import netcdf reading librairies
try:
    import netCDF4
except:
    print(''' netCDF4 is not available on this machine,
          only netcdf3 can be read ''')
    # reading and writing netcdf functions in rw_data.py won't work'''
