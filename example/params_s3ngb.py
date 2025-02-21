# -----------------------#
# Files and directories
# -----------------------#
## -- Get the user home directory
from os.path import expanduser
import os
_home = expanduser("~")
home = os.path.join(_home, "src")
home = '/mnt/data_12t'
# ------ Directory that contains orbit file:
dir_setup = os.path.join(home,  's3ng', 'data')
# ------ Directory that contains your own inputs:
indatadir = '/mnt/data_12t/llc2160_daily_latlon_SSH_notides'
# ------ Directory that contains your outputs:
working_directory = os.path.join(home, 's3ng')
# ------ Orbit file:
# Order of columns (lon, lat, time) in the orbit file
# (default is (1, 2, 0) with order_orbit_col = None)
ephemeris_cols = [0, 1, 2]
# Name of the orbit file
satname = "s3b"
ephemeris = os.path.join(dir_setup, f'orb_{satname}.txt')
# ------ Number of days in one cycle
cycle_duration = 27
# ------ Satellite elevation
height = 814.5 * 10**3
# ------ Name of the configuration (to build output files names)
config = "S3ng"
#Number of processors to be used
proc_number = 30
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
modelbox = [-180, 180, -90, 90] #None  # [230.144,234.598,42.27,47.8283]
area = [-180, 180, -90, 90]
# ------ Distance between the nadir and the end of the swath (in km):
half_swath = 80.5
# ------ Distance between the nadir and the beginning of the swath (in km):
half_gap = 5.5
# ------ Along track resolution (in km):
delta_al = 5.
# ------ Across track resolution (in km):
delta_ac = 5.
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
first_time = '2020-01-21T00:00:00Z'
# ------ Grid file name
file_grid_model = '/mnt/data/mitgcm/SSC/llc2160_2020-03-24T000000_SSU-SSV.nc'
# ------ Type of grid:
#        'regular' or 'irregular', if 'regular' only 1d coordinates
#        are extracted from model
grid = 'regular'
# ------ Specify list of variables, using the format: {key: [variable_name,
#        file_suffix], ...}, should contain at least the key 'ssh_true':
list_input_var = {'ssh_true': ['SSH_notides', 'SSH_notides']} # 'U': ['U', 'SSU-SSV'],
                  #'V': ['V', 'SSU-SSV']}
# ------ Specify factor to convert SSH values in m:
SSH_factor = 1.
# ------ Specify longitude variable:
lon = 'lon'
# ------ Specify latitude variable:
lat = 'lat'
# ------ Specify number of time in each file:
dim_time = 24
# ------ Time step between two model outputs (in days):
timestep = 1./24
# ------ Number of outputs to consider:
#        (timestep*nstep=total number of days)
nstep = 24*365
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
product_type = 'classic'

# -----------------------#
# SWOT error parameters
# -----------------------#
noise = ["altimeter", "karins3ng", "wet_troposphere", "systematic_errors3ng"]  # "roll_phase", "baseline_dilation", "timing",
         #"wet_troposphere"]
# ------ KaRIN file containing spectrum for several SWH:
#karin_noise = os.path.join(dir_setup, 'karin_noise.nc')
karin_noise = '/mnt/data_12t/swot/Input_syst_errors/SAOOH_random_error_optimB1_step2.nc'
file_systematic = '/mnt/data_12t/swot/Input_syst_errors/S3A_1year_baseline_and_phase_error.nc'

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
nbeam = 3
# ------ Beam position if there are 2 beams (in km from nadir):
beam_position = [-45, 0, 45]
central_pixel = False
