# -----------------------#
# Files and directories
# -----------------------#
## -- Get the user home directory
from os.path import expanduser
import os
home = expanduser("~")
# ------ Directory that contains orbit file:
dir_setup = os.path.join(home, 'swotsimulator', 'data')
# ------ Directory that contains your own inputs:
indatadir = os.path.join(home, 'swotsimulator', 'example',
                         'input_fields')
# ------ Directory that contains your outputs:
outdatadir = os.path.join(home, 'swotsimulator', 'example',
                          'swot_output')
# ------ Orbit file:
# Order of columns (lon, lat, time) in the orbit file
# (default is (1, 2, 0) with order_orbit_col = None)
order_orbit_col = None
# Name of the orbit file
filesat = [dir_setup + os.sep + 'ephem_science_sept2015_ell.txt']
# , dir_setup+os.sep+'orbjason.txt', dir_setup+os.sep+'orbaltika.txt' ]
# ------ Name of the configuration (to build output files names) 
config="OREGON"

# -----------------------# 
# SWOT swath parameters 
# -----------------------# 
# ------ Satellite grid file root name:
# 	 (Final file name is root_name_[numberofpass].nc)
filesgrid = outdatadir + os.sep + config + '_grid'
# ------ Force the computation of the satellite grid:
makesgrid = True
# ------ Give a subdomain if only part of the model is needed:
#	 (modelbox=[lon_min, lon_max, lat_min, lat_max])
# 	 (If modelbox is None, the whole domain of the model is considered)
modelbox = None  # [230.144,234.598,42.27,47.8283] 
# ------ Along track resolution (in km):
delta_al = 6.
# ------ Shift longitude of the orbit file if no pass is in the domain 
#        (in degree): Default value is None (no shift)
shift_lon = None
# ------ Shift time of the satellite pass (in day):
#        Default value is None (no shift)
shift_time = None

# -----------------------#
# Model input parameters
# -----------------------#
# ------ List of model files:
#	 (The first file contains the grid and is not considered as model data)
#        To generate the noise alone, file_input=None and specify region 
#        in modelbox
file_input = indatadir + os.sep + 'list_of_file.txt'
# ------ Type of model data: 
#	 (Optional, default is NETCDF_MODEL and reads netcdf3 and netcdf4 files)
#	 (Other options are ROMS, NEMO and CLS to read Nemo, roms or CLS)
model = 'NETCDF_MODEL'
# ------ Type of grid: 
#        'regular' or 'irregular', if 'regular' only 1d coordinates 
#        are extracted from model       
grid = 'irregular'
# ------ Specify SSH variable:
var = 'H'
# ------ Specify factor to convert SSH values in m:
SSH_factor = 1.
# ------ Specify longitude variable:
lon = 'lon_rho'
# ------ Specify latitude variable:
lat = 'lat_rho'
# ------ Time step between two model outputs (in days):
timestep = 1.
# ------ Number of outputs to consider:
#        (timestep*nstep=total number of days)
nstep = 50.
# ------ Not a number value:
model_nan = 0.

# -----------------------# 
# NADIR output files  
# -----------------------# 
# ------ Output file root name:
#	 (Final file name is root_name_c[cycle].nc
file_output = outdatadir + os.sep + config
# ------ Interpolation of the SSH from the model (if grid is irregular and 
#         pyresample is not installed:
#        (either 'linear' or 'nearest', use 'nearest' for large region
#        as it is faster and use less memory.)
interpolation = 'nearest'
# -----------------------# 
# NADIR error parameters 
# -----------------------# 
# ------ File containing random coefficients to compute and save 
#	 random error coefficients so that runs are reproducible:
#        If file_coeff is specified and does not exist, file is created
#	 If you don't want runs to be reproducible, file_coeff is set to None
file_coeff = None  # outdatadir+os.sep+'Random_coeff.nc'
# ------ Number of random realisations for instrumental and geophysical error 
#        (recommended ncomp=2000), ncomp1d is used for 1D spectrum, and ncomp2d
#        is used for 2D spectrum (wet troposphere computation):
ncomp1d = 3000
ncomp2d = 2000

## -- Geophysical error
## ----------------------
# ------ Wet tropo error (True to compute it):
wet_tropo = True
# ------ Beam print size (in km):
#        Gaussian footprint of sigma km
sigma = 8.
# ------ Across track resolution of the beam for the correction of the 
#        wet tropo (in km):
delta_ac = 6.
