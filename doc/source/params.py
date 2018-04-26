# -----------------------#
# Files and directories
# -----------------------#

# ------ Directory that contains orbit file:
dir_setup = '[yourpath]/SWOT_simulator/data/'
# ------ Directory that contains your own inputs:
indatadir = '[yourpath_to_yourdata]/'
# ------ Directory that contains your outputs:
outdatadir = '[yourpath_to_outputs]/'
# ------ Orbit file:
# Order of columns (lon, lat, time) in the orbit file
# (default is (1, 2, 0) with order_orbit_col = None)
order_orbit_col = None
satname = [chosenorbit]
filesat = dir_setup+'/orbit292.txt'
# ------ Name of the configuration (to build output files names)
config = [yourconfigname]
# ------ Number of processors (should not be too large to avoid MemoryError)
proc_number = [number of proc]

# -----------------------#
# SWOT swath parameters
# -----------------------#

# ------ Satellite grid file root name:
#        (Final file name is root_name_p[numberofpass].nc)
filesgrid = os.path.join(outdatadir, '{}_{}_grid'.format(config,satname))
or filesgrid = outdatadir+'/'+'[your_grid_root_name]'
# ------ Force the computation of the satellite grid:
makesgrid = True or False
# ------ Give a subdomain if only part of the model is needed:
#        (modelbox=[lon_min, lon_max, lat_min, lat_max])
#        (If modelbox is None, the whole domain of the model is considered)
modelbox = None or [yourlon_min, yourlon_max, yourlat_min, yourlat_max]
# ------ Distance between the nadir and the end of the swath (in km):
halfswath = 60.
# ------ Distance between the nadir and the beginning of the swath (in km):
halfgap = 10.
# ------ Along track resolution (in km):
delta_al = 2.
# ------ Across track resolution (in km):
delta_ac = 2.
# ------ Shift longitude of the orbit file if no pass is in the domain (in degree):
#        Default value is None (no shift)
shift_lon = None
# ------ Shift time of the satellite pass (in day):
#        Default value is None (no shift)
shift_time = None

# -----------------------#
# Model input parameters
# -----------------------#
# ------ List of model files:
#        (The first file contains the grid and is not considered as model data)
#        To generate the noise alone, file_input = None and specify region in modelbox
file_input = os.path.join(indatadir, [your_list_of_file_name.txt]) or None
# ------ Type of model data:
#        (Optional, default is NETCDF_MODEL and reads netcdf3 and netcdf4 files)
model = 'NETCDF_MODEL'
# ------ Type of grid:
# 'regular' or 'irregular', if 'regular' only 1d coordinates are extracted from model
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
nstep = 25.
# ------ Not a number value:
model_nan = 0.

# -----------------------#
# SWOT output files
# -----------------------#
# ------ Output file root name:
#        (Final file name is root_name_c[cycle]_p[pass].nc
file_output = os.path.join(outdatadir, '{}_{}'.format(config, satname))
or file_output = os.path.join(outdatadir, [your_output_root_name])
# ------ Interpolation of the SSH from the model (if grid is irregular and 
#         pyresample is not installed:
#        (either 'linear' or 'nearest', use 'nearest' for large region
#        as it is faster and use less memory.)
interpolation = 'nearest' or 'linear'
# ------ Save variables with all mockup variables ('all'), only swotsimulator
#        variables ('classic', default behaviour) or in expert mode ('expert')
save_variables = 'classic' or 'all' or 'expert'


# -----------------------#
# SWOT error parameters
# -----------------------#
# ------ File containing random coefficients to compute and save
#        random error coefficients so that runs are reproducible:
#        If file_coeff is specified and does not exist, file is created
#        If you don't want runs to be reproducible, file_coeff is set to None
file_coeff = os.path.join(outdatadir,'Random_coeff.nc') or None
# ------ KaRIN noise (True to compute it):
karin = True or False
# ------ KaRIN file containing spectrum for several SWH:
karin_file = os.path.join(dir_setup, 'karin_noise.nc')
# ------ SWH for the region:
#        if swh greater than 7 m, swh is set to 7
swh = 2.0
# ------ Number of km of random coefficients for KaRIN noise 
#        (recommended nrandkarin=1000):
nrandkarin = 1000

## -- Other instrument error (roll, phase, baseline dilation, timing)
## -----------------------------------------------------------------
# -- Compute nadir (True or False):
nadir = True
# ------ File containing spectrum of instrument error:
file_inst_error = os.path.join(dir_setup, "global_sim_instrument_error.nc)
# ------ Number of random realisations for instrumental and geophysical error (recommended ncomp=2000), ncomp1d is used for 1D spectrum, and ncomp2d is used for 2D spectrum (wet troposphere computation):
ncomp1d = 2000
ncomp2d = 2000
# ------ Cut off frequency:
#	 (Use lambda_cut=40000km for cross-calibration)
lambda_cut = 20000
lambda_max = lambda_cut
# ------ Save entire rando signal instead of random coefficients. Enables a better randomness
savesignal = True or False
# ------ If savesignal is True, enter number of pseudo-period of superimposed
#        signals and repeat length
npseudoper = 30.
len_repeat = 40000*14*50.
# ------ Roll error (True to compute it):
roll = True or False
# ------ Phase error (True to compute it):
phase = True or False
# ------ Baseline dilation error (True to compute it):
baseline_dilation = True or False
# ------ Timing error (True to compute it):
timing = True or False

## -- Geophysical error
## ----------------------
# ------ Wet tropo error (True to compute it):
wet_tropo = True or False
# ------ Beam print size (in km):
#        Gaussian footprint of sigma km
sigma = 8.
# ------ Number of beam used to correct wet_tropo signal (1, 2 or 'both'):
nbeam = 1 or 2 or ‘both’
# ------ Beam position if there are 2 beams (in km from nadir):
beam_pos_l = -35.
beam_pos_r = 35.
# ------ Sea State Bias error (Not implemented yet):
ssb = False
