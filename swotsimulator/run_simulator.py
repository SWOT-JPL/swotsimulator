'''Main program:
Usage: run_simulator(file_param)  \n
If no param file is specified, the default one is exemple/params_exemple.txt \n
In the first part of the program, model coordinates are read and the
SWOT swath is computing accordingly. \n
The SWOT grid parameters are saved in netcdf files, if you don't want to
recompute them, set maksgrid (in params file) to False.\n

In the second part of the program, errors are computed on SWOT grid for
each pass, for each cycle. The error free SSH is the SSH interpolated
from the model at each timestep. Note that there is no temporal interpolation
between model files and thus if several files are used in the SSH
interpolation, some discontinuities may be visible. \n

OUTPUTS are netcdf files containing the requested errors, the error free
SSH and the SSH with errors. There is one file every pass and every cycle.

\n
\n
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
# All rights reserved. Based on Government Sponsored Research under
# contracts NAS7-1407 and/or NAS7-03001.
#
#-----------------------------------------------------------------------
'''
import os
import shutil
from scipy import interpolate
import numpy
import glob
import sys
import logging
# Define logger level for debug purposes
logger = logging.getLogger(__name__)
try:
    import params as p
except:
    if os.path.isfile('params.py'):
        logger.error("There is a wrong entry in your params file")
        import params
    else:
        logger.error("Error: No params.py module found")
        sys.exit(1)
import swotsimulator
import swotsimulator.build_swath as build_swath
import swotsimulator.rw_data as rw_data
import swotsimulator.build_error as build_error
import swotsimulator.mod_tools as mod_tools
import swotsimulator.const as const


# - Define global variables for progress bars
istep = 0
ntot = 1
ifile = 0


def run_simulator(file_param):

    # - Initialize some parameters values
    p.shift_lon = getattr(p, 'shift_lon', None)
    p.shift_time = getattr(p, 'p.shift_time', None)
    p.timeshift = getattr(p, 'p.timeshift', 0)
    if p.shift_time is None:
        p.timeshift = 0
    model = getattr(p, 'model', 'NETCDF_MODEL')
    p.model = model
    p.model_nan = getattr(p, 'model_nan', 0)
    p.SSH_factor = getattr(p, 'SSH_factor', 1.)
    p.nadir = getattr(p, 'nadir', True)
    p.grid = getattr(p, 'grid', 'regular')

    # - Progress bar variables are global
    global istep
    global ntot

    # - Read list of user model files """
    if p.file_input is not None:
        list_file = [line.strip() for line in open(p.file_input)]
    else:
        list_file = None
    # - Read model input coordinates '''
    # if no modelbox is specified (modelbox=None), the domain of the input
    # data is taken as a modelbox
    # coordinates from the region defined by modelbox are selected
    if p.file_input is not None:
        model_data_ctor = getattr(rw_data, model)
        model_data = model_data_ctor(p, file=os.path.join(p.indatadir,
                                     list_file[0]))
    if p.modelbox is not None:
        modelbox = numpy.array(p.modelbox, dtype='float')
        # Use convert to 360 data
        modelbox[0] = (modelbox[0]+360) % 360
        if modelbox[1] != 360:
            modelbox[1] = (modelbox[1]+360) % 360
    else:
        if p.file_input is not None:
            modelbox = model_data.calc_box()
        else:
            logger.error('modelbox should be provided if no model file is'\
                         'provided')
            sys.exit(1)
    if p.file_input is not None:
        model_data.read_coordinates()
        # Select model data in the region modelbox
        if p.grid == 'regular':
            if modelbox[0] < modelbox[1]:
                model_data.model_index_lon = numpy.where(((modelbox[0]-1) <= model_data.vlon)
                                                         & (model_data.vlon <= (modelbox[1]+1)))[0]
            else:
                model_data.model_index_lon = numpy.where(((modelbox[0]-1) <= model_data.vlon)
                                                         | (model_data.vlon <= (modelbox[1]+1)))[0]
            model_data.model_index_lat = numpy.where(((modelbox[2]-1) <= model_data.vlat)
                                                     & (model_data.vlat <= (modelbox[3]+1)))[0]
            model_data.vlon = model_data.vlon[model_data.model_index_lon]
            model_data.vlat = model_data.vlat[model_data.model_index_lat]

        else:
            if modelbox[0] < modelbox[1]:
                model_data.model_index = numpy.where(((modelbox[0]-1) <= model_data.vlon)
                                                     & (model_data.vlon <= (modelbox[1]+1))
                                                     & ((modelbox[2]-1) <= model_data.vlat)
                                                     & (model_data.vlat <= (modelbox[3]+1)))
            else:
                model_data.model_index = numpy.where(((modelbox[0]-1) <= model_data.vlon)
                                                     | (model_data.vlon <= (modelbox[1]+1))
                                                     & ((modelbox[2]-1) <= model_data.vlat)
                                                     & (model_data.vlat <= (modelbox[3]+1)))

            model_data.vlon = model_data.vlon[model_data.model_index]
            model_data.vlat = model_data.vlat[model_data.model_index]
        model_data.model = model
        model_data.vloncirc = numpy.rad2deg(numpy.unwrap(model_data.vlon))
    if modelbox[1] == 0:
        modelbox[1] = 359.99
    # - Make SWOT grid if necessary """
    if p.makesgrid==True:
        logger.info('\n Force creation of SWOT grid')
        orb = build_swath.makeorbit(modelbox, p, orbitfile=p.filesat)
        build_swath.orbit2swath(modelbox, p, orb)
        logger.info("\n SWOT Grids and nadir tracks have been written in "\
                    "{}".format(p.outdatadir))
        logger.info("-----------------------------------------------")

    # - Initialize random coefficients that are used to compute
    #   random errors following the specified spectrum
    err, errnad = load_error(p)

    # - Compute interpolated SSH and errors for each pass, at each
    #   cycle
    logger.info('Compute interpolated SSH and errors:')
    #   load all SWOT grid files (one for each pass)
    listsgridfile = sorted(glob.glob(p.filesgrid + '_p*.nc'))
    if not listsgridfile:
        logger.error('\n There is no SWOT grid file in {}, run simulator with'\
                    'option makesgrid set to true in your params'\
                    'file'.format(p.outdatadir))
        sys.exit(1)
    # Build model time steps from parameter file
    modeltime = numpy.arange(0, p.nstep*p.timestep, p.timestep)
    #   Remove the grid from the list of model files
    if p.file_input:
        list_file.remove(list_file[0])
        if len(modeltime) > len(list_file):
            logger.error('There is not enough model files in the list of files')
            sys.exit(1)
    #   Initialize progress bar variables
    istep = 0
    ntot = 1

    # - Loop on SWOT grid files
    for sgridfile in listsgridfile:
        #   Load SWOT grid files (Swath and nadir)
        sgrid = load_sgrid(sgridfile, p)
        sgrid.gridfile = sgridfile
        if p.nadir is True:
            ngrid = load_ngrid(sgridfile, p)
            ngrid.gridfile = sgridfile
        else:
            ngrid = None
    #   Select model data around the swath to reduce interpolation cost in
    #   griddata

    # - Generate SWOT like and nadir-like data:
    #   Compute number of cycles needed to cover all nstep model timesteps
        rcycle = (p.timestep * p.nstep)/float(sgrid.cycle)
        ncycle = int(rcycle)
    #   Loop on all cycles
        for cycle in range(0, ncycle+1):
            if ifile > (p.nstep/p.timestep + 1):
                break
            #   Create SWOT-like and Nadir-like data
            if not p.file_input:
                model_data = []

            SSH_true, SSH_true_nadir, vindice, vindice_nadir, time, progress = create_SWOTlikedata(
                    cycle, numpy.shape(listsgridfile)[0]*rcycle, list_file,
                    modelbox, sgrid, ngrid, model_data, modeltime, err, errnad,
                    p, progress_bar=True)
            #   Save outputs in a netcdf file
            if (~numpy.isnan(vindice)).any() or not p.file_input:
                save_SWOT(cycle, sgrid, err, p, time=time, vindice=vindice,
                          SSH_true=SSH_true)
                if p.nadir is True:
                    save_Nadir(cycle, ngrid, errnad, err, p, time=time,
                               vindice_nadir=vindice_nadir,
                               SSH_true_nadir=SSH_true_nadir)
            del time
            # if p.file_input: del index
        sgrid.lon = (sgrid.lon + 360) % 360
        if p.nadir is True:
            ngrid.lon = (ngrid.lon + 360) % 360
        if p.file_input is not None:
            model_data.vlon = (model_data.vlon + 360) % 360
        modelbox[0] = (modelbox[0] + 360) % 360
        modelbox[1] = (modelbox[1] + 360) % 360
        del sgrid
        if p.nadir is True:
            del ngrid
    if progress != 1:
        progress = mod_tools.update_progress(1,
                                             'All passes have been processed',
                                             '')
    # - Write Selected parameters in a txt file
    rw_data.write_params(p, os.path.join(p.outdatadir,
                         'swot_simulator.output'))
    logger.info("\n Simulated swot files have been written in {}".format(
                 p.outdatadir))
    logger.info("----------------------------------------------------------")


def run_nadir(file_param):
    if os.path.isfile(file_param):
        #    basedir=os.path.dirname(swotsimulator.__file__)
        shutil.copyfile(file_param, 'params.py')
    else:
        logger.error("Error: No such file: {}".format(file_param))
        sys.exit()
    try:
        import params as p
    except:
        if os.path.isfile('params.py'):
            logger.error("There is a wrong entry in your params file")
            import params
        else:
            logger.error("Error: No params.py module found")
            sys.exit()

    # - Initialize some parameters values
    p.shift_lon = getattr(p, 'shift_lon', None)
    p.shift_time = getattr(p, 'p.shift_time', None)
    model = getattr(p, 'model', 'NETCDF_MODEL')
    p.model = model
    p.model_nan = getattr(p, 'model_nan', 0)
    p.SSH_factor = getattr(p, 'SSH_factor', 1.)
    p.nadir = getattr(p, 'nadir', True)
    p.grid = getattr(p, 'grid', 'regular')
    p.karin = False
    p.phase = False
    p.roll = False
    p.baseline_dilation = False
    p.timing = False
    p.halfswath = 60.
    # - Progress bar variables are global
    global istep
    global ntot
    # - Read list of user model files """
    if p.file_input is not None:
        list_file = [line.strip() for line in open(p.file_input)]
        if len(modeltime) > len(list_file):
            logger.error('There is not enough model files in the list of files')
            sys.exit(1)
    else:
        list_file = None

    ##############################################
    # Select the spatial and temporal domains
    ##############################################

    # - Read model input coordinates '''
    # if no modelbox is specified (modelbox=None), the domain of the input data
    # is taken as a modelbox
    # coordinates from the region defined by modelbox are selected
    logger.debug('Read input')
    if p.file_input is not None:
        #model_data = eval('rw_data.' + model
        #                  + '(file=os.path.join(p.indatadir, list_file[0])'))
        model_data_ctor = getattr(rw_data, model)
        model_data = model_data_ctor(file=os.path.join(p.indatadir,
                                            list_file[0]))
    if p.modelbox is not None:
        modelbox = numpy.array(p.modelbox, dtype='float')
        # Use convert to 360 data
        modelbox[0] = (modelbox[0] + 360) % 360
        if modelbox[1] != 360:
            modelbox[1] = (modelbox[1] + 360) % 360
    else:
        if p.file_input is not None:
            modelbox = model_data.calc_box()
        else:
            logger.error('modelbox should be provided if no model file is'\
                         'provided')
            sys.exit()
    logger.debug(p.file_input)
    if p.file_input is not None:
        model_data.read_coordinates()
        # Select model data in the region modelbox
        if p.grid == 'regular':
            model_data.model_index_lon = numpy.where(((modelbox[0]-1) <= model_data.vlon)
                                                     & (model_data.vlon <= (modelbox[1]+1)))[0]
            model_data.model_index_lat = numpy.where(((modelbox[2]-1) <= model_data.vlat)
                                                     & (model_data.vlat <= (modelbox[3]+1)))[0]
            model_data.vlon = model_data.vlon[model_data.model_index_lon]
            model_data.vlat = model_data.vlat[model_data.model_index_lat]

        else:
            model_data.model_index = numpy.where(((modelbox[0]-1) <= model_data.vlon)
                                                 & (model_data.vlon <= (modelbox[1]+1))
                                                 & ((modelbox[2]-1) <= model_data.vlat)
                                                 & (model_data.vlat <= (modelbox[3]+1)))
        model_data.model = model
        model_data.vloncirc = numpy.rad2deg(numpy.unwrap(model_data.vlon))
    # Ugly trick when model box is [0 360] to avoid box being empty (360=0%360)
    if modelbox[1] == 0:
        modelbox[1] = 359.99

    # - Initialize random coefficients that are used to compute
    #   random errors following the specified spectrum
    err, errnad = load_error(p)

    # - Compute interpolated SSH and errors for each pass, at each
    #   cycle
    logger.info('Compute interpolated SSH and errors:')
    #   Model time step
    modeltime = numpy.arange(0, p.nstep*p.timestep, p.timestep)
    #   Remove the grid from the list of model files
    if p.file_input:
        list_file.remove(list_file[0])
    #   Initialize progress bar variables
    istep = 0
    ntot = 1
    #   Initialize list of satellites
    istring = len(p.dir_setup)
    if not isinstance(p.filesat, list):
        p.filesat = [p.filesat]
    for filesat in p.filesat:
        # Select satellite
        #ntmp, nfilesat = os.path.split(filesat[istring:-4])
        nfilesat = os.path.basename[os.path.split(filesat)[-1]]
        # Make satellite orbit grid
        if p.makesgrid is True:
            logger.warning('\n Force creation of satellite grid')
            ngrid = build_swath.makeorbit(modelbox, p, orbitfile=filesat)
            ngrid.file = '{}{}_grid.nc'.format((p.filesgrid).strip(),
                                               nfilesat.strip())
            ngrid.write_orb()
            ngrid.ipass = nfilesat
            ngrid.gridfile = '{}{}_grid.nc'.format((p.filesgrid).strip(),
                                                   nfilesat.strip())
        else:
            # To be replaced by load_ngrid
            ngrid.gridfile = '{}{}_grid.nc'.format((p.filesgrid).strip(),
                                                   nfilesat.strip())
            ngrid = rw_data.Sat_nadir(file=ngrid.gridfile)
            ngrid.ipass = nfilesat
            cycle = 0
            x_al = []
            al_cycle = 0
            timeshift = 0
            ngrid.load_orb(cycle=cycle, x_al=x_al, al_cycle=al_cycle,
                           timeshift=timeshift)
            ngrid.loncirc = numpy.rad2deg(numpy.unwrap(ngrid.lon))
            # ngrid=load_ngrid(sgridfile, p)
            # Select model data around the swath to reduce interpolation
            # cost in griddata
            # - Generate SWOT like and nadir-like data:
        #   Compute number of cycles needed to cover all nstep model timesteps
        if p.file_input is not None:
            model_index = numpy.where(((numpy.min(ngrid.lon)) <= model_data.vlon)
                             & (model_data.vlon <= (numpy.max(ngrid.lon)))
                             & ((numpy.min(ngrid.lat)) <= model_data.vlat)
                             & (model_data.vlat <= (numpy.max(ngrid.lat))))
        rcycle = (p.timestep * p.nstep)/float(ngrid.cycle)
        ncycle = int(rcycle)
        #   Loop on all cycles
        for cycle in range(0, ncycle+1):
            if ifile > (p.nstep/p.timestep + 1): break
            #   Create SWOT-like and Nadir-like data
            if p.file_input is None:
                model_data = []
            logger.debug('compute SSH nadir')
            SSH_true_nadir, vindice, time, progress = create_Nadirlikedata(
                           cycle, numpy.shape(p.filesat)[0]*rcycle, list_file,
                           modelbox, ngrid, model_data, modeltime, errnad, p,
                           progress_bar=True)
            # SSH_true_nadir, vindice_nadir=create_Nadirlikedata(cycle, sgrid,
            # ngrid, model_data, modeltime, err, errnad, p)
##   Save outputs in a netcdf file
            ngrid.gridfile = filesat
            if (~numpy.isnan(vindice)).any() or p.file_input is None:
                err=errnad
                err.wtnadir = numpy.zeros((1))
                err.wet_tropo2nadir = numpy.zeros((1))
                logger.debug('write file')
                save_Nadir(cycle, ngrid, errnad, err, p,time=time,
                           vindice_nadir=vindice,
                           SSH_true_nadir=SSH_true_nadir)
            del time
            # if p.file_input: del index
        ngrid.lon = (ngrid.lon + 360) % 360
        if p.file_input:
            model_data.vlon = (model_data.vlon + 360) % 360
        modelbox[0] = (modelbox[0] + 360) % 360
        modelbox[1] = (modelbox[1] + 360) % 360
        del ngrid
    if progress != 1:
        progress=mod_tools.update_progress(1, 'All passes have been processed',
                                           '')
    # - Write Selected parameters in a txt file
    rw_data.write_params(p, os.path.join(p.outdatadir,
                         'nadir_simulator.output'))
    logger.info("\nSimulated orbit files have been written in {}".format(
                p.outdatadir))
    logger.info("----------------------------------------------------------")


def load_error(p):
    '''Initialize random coefficients that are used to compute
    random errors following the specified spectrum. \n
    If a random coefficient file is specified, random coefficients
    are loaded from this file.
    '''
    err = build_error.error(p)
    if p.nadir is True:
        errnad = build_error.errornadir(p)
    else:
        errnad = None
    try:
        nhalfswath = int((p.halfswath-p.halfgap)/p.delta_ac) + 1
    except AttributeError:
        nhalfswath=60.
    if p.file_coeff:
        if os.path.isfile(p.file_coeff) and (not p.makesgrid):
            logger.warn('\n WARNING: Existing random coefficient file used')
            err.load_coeff(p)
        else:
            err.init_error(p, 2*nhalfswath)
            err.save_coeff(p, 2*nhalfswath)
        if p.nadir is True:
            if os.path.isfile(p.file_coeff[:-3] + '_nadir.nc') \
                    and (not p.makesgrid):
                logger.warn('WARNING: Existing random nadir coefficient file'\
                            'used')
                errnad.load_coeff(p)
            else:
                errnad.init_error(p)
                errnad.save_coeff(p)
    else:
        err.init_error(p, 2*nhalfswath)
        if p.nadir is True:
            errnad.init_error(p)
    return err, errnad


def load_sgrid(sgridfile, p):
    '''Load SWOT swath and Nadir data for file sgridfile '''

    # Load SWOT swath file
    sgrid = rw_data.Sat_SWOT(file=sgridfile)
    cycle = 0
    x_al = []
    x_ac = []
    al_cycle = 0
    timeshift = 0
    sgrid.load_swath(cycle=cycle ,x_al=x_al, x_ac=x_ac, al_cycle=al_cycle,
                     timeshift=timeshift)
    sgrid.loncirc = numpy.rad2deg(numpy.unwrap(sgrid.lon))
    # Extract the pass number from the file name
    ipass = int(sgridfile[-6: -3])
    sgrid.ipass=ipass
    return sgrid

def interpolate_regular_1D(lon_in, lat_in, var, lon_out, lat_out, Teval=None):
    ''' Interpolation of data when grid is regular and coordinate in 1D. '''
    # To correct for IDL issues
    if numpy.max(lon_in)>359 and numpy.min(lon_in)<1:
        lon_in[lon_in > 180] = lon_in[lon_in > 180] - 360
        #lon_in = np.mod(lon_in - (lref - 180), 360) + (lref - 180)
        lon_in = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(lon_in)))
    # Interpolate mask if it has not been done (Teval is None)
    if Teval is None:
        Teval = interpolate.RectBivariateSpline(lat_in, lon_in,
                                            numpy.isnan(var),
                                            kx=1, ky=1, s=0).ev(lat_out,
                                            lon_out)
    # Trick to avoid nan in interpolation
    var_mask = + var
    var_mask[numpy.isnan(var_mask)] = 0.
    # Interpolate variable
    var_out = interpolate.RectBivariateSpline(lat_in, lon_in, var_mask, kx=1,
                                              ky=1,
                                              s=0).ev(lat_out, lon_out)
    # Mask variable with Teval
    var_out[Teval > 0] = numpy.nan
    return var_out, Teval

def interpolate_irregular_pyresample(swath_in, var, grid_out, radius,
                                     interp_type='nearest'):
    ''' Interpolation of data when grid is irregular and pyresample is
    installed.'''
    import pyresample as pr
    if interp_type == 'nearest':
        var_out = pr.kd_tree.resample_nearest(swath_in, var, grid_out,
                                          radius_of_influence=radius
                                          *10**3, epsilon=100)
    else:
        var_out = pr.kd_tree.resample_gauss(swath_in, var, grid_out,
                                          radius_of_influence=3*radius *10**3,
                                          sigmas=radius*10**3, fill_value=None)
    return var_out

def load_ngrid(sgridfile, p):
    ipass = int(sgridfile[-6: -3])
    # Load Nadir track file
    ngrid = rw_data.Sat_nadir(file='{:s}nadir_p{:03d}.nc'.format(
                                   (p.filesgrid).strip(), ipass))
    cycle = 0
    x_al = []
    al_cycle = 0
    timeshift = 0
    ngrid.load_orb(cycle=cycle, x_al=x_al, al_cycle=al_cycle,
                   timeshift=timeshift)
    ngrid.loncirc = numpy.rad2deg(numpy.unwrap(ngrid.lon))
    ngrid.ipass = ipass
    return ngrid

def select_modelbox(sgrid, model_data):
    # mask=numpy.zeros((numpy.shape(model_data.vlon)))
    # nal=len(sgrid.x_al)
    # for kk in range(0,nal,10):
    # dlon1=model_data.vlon-sgrid.lon_nadir[kk]
    # dlon=numpy.minimum(numpy.mod(dlon1,360),numpy.mod(-dlon1,360))
    # ddist=(((dlon)*numpy.cos(sgrid.lat_nadir[kk]*numpy.pi/180.)*
    # const.deg2km)**2+((model_data.vlat-sgrid.lat_nadir[kk])*const.deg2km)**2 )
    #  mask[ddist<const.radius_interp**2]=1
    model_index = numpy.where(((numpy.min(sgrid.lon)) <= model_data.vlon)
                              & (model_data.vlon <= (numpy.max(sgrid.lon)))
                              & ((numpy.min(sgrid.lat)) <= model_data.vlat)
                              & (model_data.vlat <= (numpy.max(sgrid.lat))))
    model_data.model_index = model_index
    lonmodel1d = model_data.vlon[model_index].ravel()
    latmodel1d = model_data.vlat[model_index].ravel()

    if p.grid == 'regular':
        model_data.lon1d = lon1[numpy.where(((numpy.min(sgrid.lon)) <= lon1)
                                & (lon1 <= (numpy.max(sgrid.lon))))]
        model_data.lat1d = lat1[numpy.where(((numpy.min(sgrid.lat)) <= lat1)
                                & (lat1 <= (numpy.max(sgrid.lat))))]
    else:
        model_index = numpy.where(((numpy.min(sgrid.lon)) <= model_data.vlon)
                                  & (model_data.vlon <= (numpy.max(sgrid.lon)))
                                  & (numpy.min(sgrid.lat) <= model_data.vlat)
                                  & (model_data.vlat<=(numpy.max(sgrid.lat))))
        model_data.lon1d = model_data.vlon[model_index].ravel()
        model_data.lat1d = model_data.vlat[model_index].ravel()
        model_data.vlon = model_data.vlon[model_index].ravel()
        model_data.vlat = model_data.vlat[model_index].ravel()

    nx = len(lon1)
    ny = len(lat1)
    return model_index


def create_SWOTlikedata(cycle, ntotfile, list_file, modelbox, sgrid, ngrid,
                        model_data, modeltime, err, errnad, p,
                        progress_bar=True):
    '''Create SWOT and nadir errors err and errnad, interpolate model SSH model
    _data on swath and nadir track, compute SWOT-like and nadir-like data
    for cycle, SWOT grid sgrid and ngrid. '''
    # - Progress bar variables are global
    global istep
    global ntot
    #   Initialiaze errors and SSH
    progress = 0
    err.karin = numpy.zeros((numpy.shape(sgrid.lon)[0],
                            numpy.shape(sgrid.lon)[1]))
    err.roll = numpy.zeros((numpy.shape(sgrid.lon)[0],
                         numpy.shape(sgrid.lon)[1]))
    err.phase = numpy.zeros((numpy.shape(sgrid.lon)[0],
                          numpy.shape(sgrid.lon)[1]))
    err.baseline_dilation = numpy.zeros((numpy.shape(sgrid.lon)[0],
                                        numpy.shape(sgrid.lon)[1]))
    err.timing = numpy.zeros((numpy.shape(sgrid.lon)[0],
                             numpy.shape(sgrid.lon)[1]))
    err.wet_tropo1 = numpy.zeros((numpy.shape(sgrid.lon)[0],
                                 numpy.shape(sgrid.lon)[1]))
    err.wet_tropo2 = numpy.zeros((numpy.shape(sgrid.lon)[0],
                               numpy.shape(sgrid.lon)[1]))
    err.ssb = numpy.zeros((numpy.shape(sgrid.lon)[0],
                          numpy.shape(sgrid.lon)[1]))
    err.wt = numpy.zeros((numpy.shape(sgrid.lon)[0],
                         numpy.shape(sgrid.lon)[1]))
    SSH_true = numpy.zeros((numpy.shape(sgrid.lon)[0],
                           numpy.shape(sgrid.lon)[1]))  # , p.nstep))
    # Initialize nadir variable in case the simulator is run with nadir option
    # at False
    SSH_true_nadir = 0
    vindice_nadir = 0
    if p.nadir is True:
        err.wet_tropo1nadir = numpy.zeros((numpy.shape(ngrid.lon)[0]))
        err.wet_tropo2nadir = numpy.zeros((numpy.shape(ngrid.lon)[0]))
        err.wtnadir = numpy.zeros((numpy.shape(ngrid.lon)[0]))
        errnad.nadir = numpy.zeros((numpy.shape(ngrid.lon)[0]))
        vindice_nadir = numpy.zeros(numpy.shape(ngrid.lon)) * numpy.nan
        SSH_true_nadir = numpy.zeros((numpy.shape(ngrid.lon)[0]))
    vindice = numpy.zeros(numpy.shape(SSH_true)) * numpy.nan
    # Definition of the time in the model
    date1 = cycle * sgrid.cycle
    time = sgrid.time + date1
    # Look for satellite data that are beween step-p.timestep/2 and step+p.step/2
    if p.file_input:
        index_filemodel = numpy.where(((time[-1]-sgrid.timeshift) >= (modeltime-p.timestep/2.))
                                      & ((time[0]-sgrid.timeshift) < (modeltime+p.timestep/2.)) )
        # At each step, look for the corresponding time in the satellite data
        for ifile in index_filemodel[0]:
            progress = mod_tools.update_progress(float(istep)/float(ntot
                                      * ntotfile),
                                      'pass: {}'.format(sgrid.ipass),
                                      'model file: {}, cycle: {}'.format(
                                      list_file[ifile], cycle+1))
            # If there are satellite data, Get true SSH from model
            if numpy.shape(index_filemodel)[1] > 0:
                # number of file to be processed used in the progress bar
                ntot = ntot + numpy.shape(index_filemodel)[1]-1
                # if numpy.shape(index)[1]>1:
                # Select part of the track that corresponds to the time of the model (+-timestep/2)
                ind_time = numpy.where(((time-sgrid.timeshift) >= (modeltime[ifile]-p.timestep/2.))
                                       & ((time-sgrid.timeshift) < (modeltime[ifile]+p.timestep/2.)) )
                if p.nadir is True:
                    ind_nadir_time = numpy.where(((time-ngrid.timeshift) >= (modeltime[ifile]-p.timestep/2.))
                                                 & ((time-ngrid.timeshift) < (modeltime[ifile]+p.timestep/2.)) )
                # Load data from this model file
                model_step_ctor = getattr(rw_data, model_data.model)
                model_step = model_step_ctor(p, file=os.path.join(p.indatadir,
                                            list_file[ifile]), var=p.var)
                if p.grid == 'regular':
                    model_step.read_var()
                    SSH_model = model_step.vvar[model_data.model_index_lat, :]
                    SSH_model = SSH_model[:, model_data.model_index_lon]
                else:
                    model_step.read_var(index=model_data.model_index)
                    SSH_model = model_step.vvar
                # - Interpolate Model data on a SWOT grid and/or along the
                #   nadir track
                # if grid is regular, use interpolate.RectBivariateSpline to
                # interpolate
                if p.grid == 'regular' and \
                        len(numpy.shape(model_data.vlon)) == 1:
                    # ########################TODO
                    # To be moved to routine rw_data
                    indsorted = numpy.argsort(model_data.vlon)
                    model_data.vlon = model_data.vlon[indsorted]
                    SSH_model = SSH_model[:, indsorted]
                    # Flatten satellite grid and select part of the track
                    # corresponding to the model time
                    lonswot = sgrid.lon[ind_time[0], :].flatten()
                    latswot = sgrid.lat[ind_time[0], :].flatten()
                    SSH_true_ind_time, Teval = interpolate_regular_1D(
                                             model_data.vlon,
                                             model_data.vlat,
                                             SSH_model,
                                             sgrid.lon[ind_time[0], :].ravel(),
                                             sgrid.lat[ind_time[0], :].ravel())
                    nal, nac = numpy.shape(sgrid.lon[ind_time[0], :])
                    SSH_true[ind_time[0], :] = SSH_true_ind_time.reshape(nal,
                                                                         nac)
                    if p.nadir is True:
                        SSH_true_ind_time, Teval = interpolate_regular_1D(
                                          model_data.vlon,
                                          model_data.vlat,
                                          SSH_model,
                                          ngrid.lon[ind_nadir_time[0]].ravel(),
                                          ngrid.lat[ind_nadir_time[0]].ravel())
                        SSH_true_nadir[ind_nadir_time[0]] = SSH_true_ind_time
                else:
                    # Grid is irregular, interpolation can be done using
                    # pyresample module if it is installed or griddata
                    # function from scipy.
                    # Note that griddata is slower than pyresample functions.
                    try:
                        import pyresample as pr
                        model_data.vlon = pr.utils.wrap_longitudes(
                                                               model_data.vlon)
                        sgrid.lon = pr.utils.wrap_longitudes(sgrid.lon)
                        if len(numpy.shape(model_data.vlon)) <= 1:
                            model_data.vlon, model_data.vlat = numpy.meshgrid(
                                                               model_data.vlon,
                                                               model_data.vlat)
                        swath_def = pr.geometry.SwathDefinition(
                                    lons=model_data.vlon, lats=model_data.vlat)
                        grid_def = pr.geometry.SwathDefinition(lons=sgrid.lon,
                                                               lats=sgrid.lat)
                        SSH_true[ind_time[0], :] = interpolate_irregular_pyresample(
                                               swath_def, SSH_model, grid_def,
                                               max(p.delta_al, p.delta_ac),
                                               interp_type=p.interpolation)
                        if p.nadir is True:
                            ngrid.lon = pr.utils.wrap_longitudes(ngrid.lon)
                            ngrid_def = pr.geometry.SwathDefinition(lons=ngrid.lon, lats=ngrid.lat)
                            SSH_true_nadir[ind_time[0], :] = interpolate_irregular_pyresample(
                                               swath_def, SSH_model, ngrid_def,
                                               max(p.delta_al, p.delta_ac),
                                               interp_type=p.interpolation)
                    except:
                        SSH_true[ind_time[0], :] = interpolate.griddata(
                            (model_data.vlon.ravel(), model_data.vlat.ravel()),
                            SSH_model.ravel(), (sgrid.lon[ind_time[0], :],
                            sgrid.lat[ind_time[0], :]), method=p.interpolation)
                        if p.nadir is True:
                            SSH_true_nadir[ind_nadir_time[0]] = interpolate.griddata((model_data.vlon.ravel(), model_data.vlat.ravel()),
                                 SSH_model.ravel(),
                                 (ngrid.lon[ind_nadir_time[0]], ngrid.lat[ind_nadir_time[0]]),
                                 method=p.interpolation)
                        if p.interpolation == 'nearest':
                            if modelbox[0] > modelbox[1]:
                                SSH_true[numpy.where(((sgrid.lon < modelbox[0])
                                         & (sgrid.lon > modelbox[1]))
                                         | (sgrid.lat < modelbox[2])
                                         | (sgrid.lat > modelbox[3]))] = numpy.nan
                                if p.nadir is True:
                                    SSH_true_nadir[numpy.where((
                                                   (ngrid.lon < modelbox[0])
                                                   & (ngrid.lon > modelbox[1]))
                                                   | (ngrid.lat < modelbox[2])
                                                   | (ngrid.lat > modelbox[3])
                                                   )] = numpy.nan
                            else:
                                SSH_true[numpy.where((sgrid.lon < modelbox[0])
                                                    | (sgrid.lon > modelbox[1])
                                                    | (sgrid.lat < modelbox[2])
                                                    | (sgrid.lat > modelbox[3])
                                                    )] = numpy.nan
                                if p.nadir is True:
                                    SSH_true_nadir[numpy.where(
                                                    (ngrid.lon < modelbox[0])
                                                    | (ngrid.lon > modelbox[1])
                                                    | (ngrid.lat < modelbox[2])
                                                    | (ngrid.lat > modelbox[3])
                                                    )] = numpy.nan
                vindice[ind_time[0], :] = ifile
                if p.nadir is True:
                    vindice_nadir[ind_nadir_time[0]] = ifile
                    del ind_time, SSH_model, model_step, ind_nadir_time
                else:
                    del ind_time, SSH_model, model_step
            istep += 1
    else:
        istep += 1
        progress = mod_tools.update_progress(float(istep)/float(ntotfile*ntot),
                                             'pass: {}'.format(sgrid.ipass),
                                             'no model file provided, cycle: '\
                                             '{}'.format(cycle+1))
    err.make_error(sgrid, cycle, SSH_true, p)
    err.make_SSH_error(SSH_true, p)
    if p.nadir is True:
        errnad.make_error(ngrid, cycle, SSH_true_nadir, p)
        if p.nbeam == 1:
            errnad.SSH = SSH_true_nadir + errnad.nadir + err.wet_tropo1nadir
        else:
            errnad.SSH = SSH_true_nadir + errnad.nadir + err.wet_tropo2nadir
    # if p.file_input: del ind_time, SSH_model, model_step
    return SSH_true, SSH_true_nadir, vindice, vindice_nadir, time, progress


def create_Nadirlikedata(cycle, ntotfile, list_file, modelbox, ngrid,
                         model_data, modeltime,  errnad, p,
                         progress_bar=False):

    # - Progress bar variables are global
    global istep
    global ntot
    global ifile
    errnad.wet_tropo1nadir = numpy.zeros((numpy.shape(ngrid.lon)[0]))
    errnad.wt = numpy.zeros((numpy.shape(ngrid.lon)[0]))
    errnad.nadir = numpy.zeros((numpy.shape(ngrid.lon)[0]))
    SSH_true_nadir = numpy.zeros((numpy.shape(ngrid.lon)[0]))
    vindice = numpy.zeros(numpy.shape(ngrid.lon))*numpy.nan
    # Definition of the time in the model
    date1 = cycle * ngrid.cycle
    time = ngrid.time + date1
    # Look for satellite data that are beween step-p.timesetp/2 end setp+p.step/2
    if p.file_input is not None:
        index_filemodel = numpy.where(((time[-1]-ngrid.timeshift) >= (modeltime-p.timestep/2.))
                                      & ((time[0]-ngrid.timeshift) < (modeltime+p.timestep/2.)))  # [0]
        # At each step, look for the corresponding time in the satellite data
        for ifile in index_filemodel[0]:
            if progress_bar:
                progress = mod_tools.update_progress(float(istep)
                                          /float(ntotfile*ntot),
                                          'orbit: {}'.format(ngrid.ipass),
                                          'model file: {}, cycle: {}'.format(
                                          list_file[ifile], cycle+1))
            else:
                progress = None
            # If there are satellite data, Get true SSH from model
            if numpy.shape(index_filemodel)[1] > 0:
                ntot = ntot + numpy.shape(index_filemodel)[1] - 1
                ind_nadir_time = numpy.where(((time-ngrid.timeshift) >= (modeltime[ifile]-p.timestep/2.))
                                             & ((time-ngrid.timeshift) < (modeltime[ifile]+p.timestep/2.)))
                model_step_ctor = getattr(rw_data, model_data.model)
                model_step = model_data_ctor(file=os.path.join(p.indatadir,
                                            list_file[ifile]), var=p.var)
                if p.grid == 'regular':
                    model_step.read_var()
                    SSH_model = model_step.vvar[model_data.model_index_lat, :]
                    SSH_model = SSH_model[:, model_data.model_index_lon]
                else:
                    model_step.read_var(index=model_data.model_index)
                    SSH_model = model_step.vvar
            # - Interpolate Model data along the nadir
            # track
            # if grid is regular, use interpolate.RectBivariateSpline to
            # interpolate
            if p.grid == 'regular' and len(numpy.shape(model_data.vlon)) == 1:
                # ########################TODO
                # To be moved to routine rw_data
                indsorted = numpy.argsort(model_data.vlon)
                model_data.vlon = model_data.vlon[indsorted]
                SSH_model = SSH_model[:, indsorted]
                SSH_true_ind_time = interpolate_regular_1D(model_data.vlon,
                                          model_data.vlat, SSH_model,
                                          ngrid.lon[ind_nadir_time[0]].ravel(),
                                          ngrid.lat[ind_nadir_time[0]].ravel())
                SSH_true_nadir[ind_nadir_time[0]] = SSH_true_ind_time
            else:
                # Grid is irregular, interpolation can be done using pyresample
                # module if it is installed or griddata function from scipy.
                # Note that griddata is slower than pyresample functions.
                try:
                    import pyresample as pr
                    ngrid.lon = pr.utils.wrap_longitudes(ngrid.lon)
                    model_data.vlon = pr.utils.wrap_longitudes(model_data.vlon)
                    ngrid_def = pr.geometry.SwathDefinition(lons=ngrid.lon,
                                                            lats=ngrid.lat)
                    swath_def = pr.geometry.SwathDefinition(
                                                          lons=model_data.vlon,
                                                          lats=model_data.vlat)
                    SSH_true_nadir[ind_nadir_time[0]] = interpolate_irregular_pyresample(
                                                 swath_defu, u_model, grid_def,
                                                 max(p.delta_al, p.delta_ac),
                                                 interp_type=p.interpolation)
                except:
                    SSH_true_nadir[ind_nadir_time[0]] = interpolate.griddata(
                                                (model_data.vlon.ravel(),
                                                model_data.vlat.ravel()),
                                                SSH_model.ravel(),
                                                (ngrid.lon[ind_nadir_time[0]],
                                                ngrid.lat[ind_nadir_time[0]]),
                                                method=p.interpolation)
                    if p.interpolation == 'nearest':
                        if modelbox[0] > modelbox[1]:
                            SSH_true_nadir[numpy.where((
                                             (ngrid.lon < modelbox[0])
                                             & (ngrid.lon > modelbox[1]))
                                             | (ngrid.lat < modelbox[2])
                                             | (ngrid.lat > modelbox[3]))] = numpy.nan
                        else:
                            SSH_true_nadir[numpy.where(
                                             (ngrid.lon < modelbox[0])
                                             | (ngrid.lon > modelbox[1])
                                             | (ngrid.lat < modelbox[2])
                                             | (ngrid.lat > modelbox[3]))] = numpy.nan
            vindice[ind_nadir_time[0]] = ifile
            istep += 1
    else:
        if progress_bar:
            progress = mod_tools.update_progress(float(istep)/float(ntotfile
                                      * ntot), 'orbit: {}'.format(ngrid.ipass),
                                      'model file: {}, cycle: {}'.format(
                                      list_file[ifile], cycle+1))
        else:
            progress = None
        istep += 1
    errnad.make_error(ngrid, cycle, SSH_true_nadir, p)  # , ind[0])
    errnad.SSH = SSH_true_nadir + errnad.nadir + errnad.wet_tropo1nadir
    # del SSH_model, model_step, ind_nadir_time
    return SSH_true_nadir, vindice, time, progress


def save_SWOT(cycle, sgrid, err, p, time=[], vindice=[], SSH_true=[]):
    file_output = (p.file_output + '_c' + str(cycle+1).zfill(2) + '_p'
                   + str(sgrid.ipass).zfill(3) + '.nc')
    OutputSWOT = rw_data.Sat_SWOT(file=file_output, lon=(sgrid.lon+360) % 360,
                                  lat=sgrid.lat, time=time, x_ac=sgrid.x_ac,
                                  x_al=sgrid.x_al, cycle=sgrid.cycle,
                                  lon_nadir=(sgrid.lon_nadir+360) % 360,
                                  lat_nadir=sgrid.lat_nadir)
    OutputSWOT.gridfile = sgrid.gridfile
    OutputSWOT.write_data(SSH_model=SSH_true, index=vindice, roll_err=err.roll,
                          bd_err=err.baseline_dilation, phase_err=err.phase,
                          ssb_err=err.ssb, karin_err=err.karin,
                          pd_err_1b=err.wet_tropo1, pd_err_2b=err.wet_tropo2,
                          pd=err.wt, timing_err=err.timing, SSH_obs=err.SSH)
    return None


def save_Nadir(cycle, ngrid, errnad, err, p, time=[], vindice_nadir=[],
               SSH_true_nadir=[]):
    file_outputnadir = (p.file_output + 'nadir_c' + str(cycle+1).zfill(2)
                        + '_p' + str(ngrid.ipass).zfill(3) + '.nc')
    OutputNadir = rw_data.Sat_nadir(file=file_outputnadir,
                                    lon=(ngrid.lon+360) % 360,
                                    lat=ngrid.lat, time=time, x_al=ngrid.x_al,
                                    cycle=ngrid.cycle)
    OutputNadir.gridfile = ngrid.gridfile
    OutputNadir.write_data(SSH_model=SSH_true_nadir, index=vindice_nadir,
                           nadir_err=errnad.nadir, SSH_obs=errnad.SSH,
                           pd_err_1b=err.wet_tropo1nadir, pd=err.wtnadir,
                           pd_err_2b=err.wet_tropo2nadir)
    return None
