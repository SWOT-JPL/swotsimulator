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
working_directory = os.path.join(home, 'swotsimulator', 'example',
                           'swot_output')
# ------ Orbit file:
# Order of columns (lon, lat, time) in the orbit file
# (default is (1, 2, 0) with order_orbit_col = None)
ephemeris_cols = None
# Name of the orbit file
satname = "science"
ephemeris = os.path.join(dir_setup, 'ephem_science_sept2015_ell.txt')
# ------ Number of days in one cycle
cycle_duration = 20.86455
# ------ Satellite elevation
height = 891 * 10**3
# ------ Name of the configuration (to build output files names)
config = "OREGON"
#Number of processors to be used
proc_number = 1
# ------ Deactivate printing of progress bar to avoid huge log
progress_bar = True

# -----------------------#
# SWOT swath parameters
# -----------------------#
# ------ Satellite grid file root name:
# 	 (Final file name is root_name_[numberofpass].nc)
filesgrid = os.path.join(working_directory, '{}_{}_grid'.format(config,satname))
# ------ Force the computation of the satellite grid:
makesgrid = True
# ------ Give a subdomain if only part of the model is needed:
#	 (modelbox=[lon_min, lon_max, lat_min, lat_max])
# 	 (If modelbox is None, the whole domain of the model is considered)
modelbox = None  # [230.144,234.598,42.27,47.8283]
area = None
# ------ Distance between the nadir and the end of the swath (in km):
half_swath = 60.
# ------ Distance between the nadir and the beginning of the swath (in km):
half_gap = 10.
# ------ Along track resolution (in km):
delta_al = 2.
# ------ Across track resolution (in km):
delta_ac = 2.
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
#        To generate the noise alone, file_input = None
#        and specify region in modelbox
file_input = os.path.join(indatadir, 'list_of_file.txt')
# ------ Type of model data:
#	 (Optional, default is NETCDF_MODEL and reads netcdf3 and netcdf4 files)
#	 (Other options are ROMS, NEMO and CLS to read Nemo, roms or CLS)
model = 'NETCDF_MODEL'
# ------ First time of the model
first_time = '2011-11-15T00:00:00Z'
# ------ Grid file name
file_grid_model = os.path.join(indatadir, 'OREGON_grd.nc')
# ------ Type of grid:
#        'regular' or 'irregular', if 'regular' only 1d coordinates
#        are extracted from model
grid = 'irregular'
# ------ Specify list of variables, using the format: {key: [variable_name,
#        file_suffix], ...}, should contain at least the key 'ssh_true':
list_input_var = {'ssh_true': ['H', '']}
# ------ Specify factor to convert SSH values in m:
SSH_factor = 1.
# ------ Specify longitude variable:
lon = 'lon_rho'
# ------ Specify latitude variable:
lat = 'lat_rho'
# ------ Specify number of time in each file:
dim_time = 1
# ------ Time step between two model outputs (in days):
timestep = 1.
# ------ Number of outputs to consider:
#        (timestep*nstep=total number of days)
nstep = 25.
# ------ Not a number value:
model_nan = 0.

# -----------------------#
# SWOT output files
# -----------------------#
# ------ Output file root name:
#	 (Final file name is root_name_c[cycle]_p[pass].nc
file_output = os.path.join(working_directory, '{}_{}'.format(config, satname))
# ------ Interpolation of the SSH from the model (if grid is irregular and
#         pyresample is not installed:
#        (either 'linear' or 'nearest', use 'nearest' for large region
#        as it is faster and use less memory.)
interpolation = 'linear'
# ------ Save variables with all mockup variables ('all'), only swotsimulator
#        variables ('classic', default behaviour) or in expert mode ('expert')
product_type = 'all'

# -----------------------#
# SWOT error parameters
# -----------------------#
noise = ["altimeter", "karin", "roll_phase", "baseline_dilation", "timing",
         "wet_troposphere"]
# ------ KaRIN file containing spectrum for several SWH:
karin_noise = os.path.join(dir_setup, 'karin_noise.nc')
# "Seed for RandomState. Must be convertible to 32 bit "
nseed = 0
# ------ SWH for the region:
#        if swh greater than 7 m, swh is set to 7
swh = 2.0
# ------ Number of km of random coefficients for KaRIN noise (recommended nrandkarin=1000):
#nrandkarin = 1000

## -- Other instrument error (roll, phase, baseline dilation, timing)
## -----------------------------------------------------------------
# -- Compute nadir (True or False):
nadir = True
# ------ File containing spectrum of instrument error:
error_spectrum = os.path.join(dir_setup, "global_sim_instrument_error.nc")
# ------ Number of random realisations for instrumental and geophysical error
#        (recommended ncomp=2000), ncomp1d is used for 1D spectrum, and ncomp2d
#        is used for 2D spectrum (wet troposphere computation):
#ncomp1d = 4000
# ncomp2d = 2000
# ------ Cut off frequency:
#	 (Use lambda_cut=40000km for cross-calibration)
lambda_cut = 20000
lambda_max = lambda_cut
# ------ If savesignal is True, enter number of pseudo-period of superimposed
#        signals and repeat length
len_repeat = 40000 #*14*50.
# Roll-phase simulation of correction file
corrected_roll_phase_dataset =  os.path.join(dir_setup,
                                             'data_sim_slope_2cycles_v0.nc')


## -- Geophysical error
## ----------------------
# ------ Beam print size (in km):
#        Gaussian footprint of sigma km
sigma = 8.
# ------ Number of beam used to correct wet_tropo signal (1, 2 or 'both'):
nbeam = 2
# ------ Beam position if there are 2 beams (in km from nadir):
beam_position = [-35, 35]



central_pixel = False
