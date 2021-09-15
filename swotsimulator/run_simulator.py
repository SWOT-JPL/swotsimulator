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
# - Dec 2017: Version 3.0
# - Apr 2018: Version 3.1
#
# Notes:
# - Tested with Python 3.5, Python 3.6, Python 3.7
#
# Copyright (c)
# Copyright (c) 2002-2014, California Institute of Technology.
# All rights reserved. Based on Government Sponsored Research under
# contracts NAS7-1407 and/or NAS7-03001.
#
#-----------------------------------------------------------------------
'''
import os
import datetime
from scipy import interpolate
import numpy
import glob
import sys
import math
import time
import traceback
import swotsimulator.build_swath as build_swath
import swotsimulator.rw_data as rw_data
import swotsimulator.build_error as build_error
import swotsimulator.mod_tools as mod_tools
import swotsimulator.const as const
import swotsimulator.mod_run as mod
import swotsimulator.mod_parallel as parallel
import multiprocessing
import logging
# Define logger level for debug purposes
logger = logging.getLogger(__name__)
#logger = multiprocessing.log_to_stderr()
#logger.setLevel(logging.DEBUG)


# - Define global variables for progress bars
istep = 0
ntot = 1
ifile = 0


def run_simulator(p, die_on_error=False, nadir_alone=False):
    '''Main routine to run simulator, input is the imported parameter file,
    no outputs are returned but netcdf grid and data files are written as well
    as a skimulator.output file to store all used parameter.
    '''
    # - Initialize some parameters values
    timestart = datetime.datetime.now()
    mod_tools.initialize_parameters(p)
    mod_tools.check_path(p)

    # - Read list of user model files """
    model_data, list_file = mod.load_coordinate_model(p)
    ## - Read model input coordinates '''
    ## if no modelbox is specified (modelbox=None), the domain of the input
    ## data is taken as a modelbox
    ## coordinates from the region defined by modelbox are selected
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
            logger.error('modelbox should be provided if no model file is'
                         'provided')
            sys.exit(1)
    p.modelbox_calc = modelbox
    if p.file_input is not None:
        model_data.read_coordinates()
        # Select model data in the region modelbox
        model_data.len_coord = len(numpy.shape(model_data.vlon))
        if p.grid == 'regular' or model_data.len_coord == 1:
            if modelbox[0] < modelbox[1]:
                _i_lon = numpy.where(((modelbox[0]-1) <= model_data.vlon)
                                     & (model_data.vlon <= (modelbox[1]+1)))[0]
            else:
                _i_lon = numpy.where(((modelbox[0]-1) <= model_data.vlon)
                                     | (model_data.vlon <= (modelbox[1]+1)))[0]
            model_data.model_index_lon = _i_lon
            _i_lat = numpy.where(((modelbox[2]-1) <= model_data.vlat)
                                 & (model_data.vlat <= (modelbox[3]+1)))[0]
            model_data.model_index_lat = _i_lat
            model_data.vlon = model_data.vlon[model_data.model_index_lon]
            model_data.vlat = model_data.vlat[model_data.model_index_lat]

        else:
            if modelbox[0] < modelbox[1]:
                _i_box = numpy.where(((modelbox[0]-1) <= model_data.vlon)
                                     & (model_data.vlon <= (modelbox[1]+1))
                                     & ((modelbox[2]-1) <= model_data.vlat)
                                     & (model_data.vlat <= (modelbox[3]+1)))
            else:
                _i_box = numpy.where(((modelbox[0]-1) <= model_data.vlon)
                                     | (model_data.vlon <= (modelbox[1]+1))
                                     & ((modelbox[2]-1) <= model_data.vlat)
                                     & (model_data.vlat <= (modelbox[3]+1)))
            model_data.model_index = _i_box
            model_data.vlon = model_data.vlon[model_data.model_index]
            model_data.vlat = model_data.vlat[model_data.model_index]
        model_data.model = p.model
        model_data.vloncirc = numpy.rad2deg(numpy.unwrap(model_data.vlon))
    if modelbox[1] == 0:
        modelbox[1] = 359.99
    # - Make SWOT grid if necessary """
    if p.makesgrid is True:
        logger.info('\n Force creation of SWOT grid')
        # make nadir orbit
        orb = build_swath.makeorbit(modelbox, p, orbitfile=p.ephemeris)
        # build swath for this orbit
        if nadir_alone is True:
            build_swath.orbit2nadir(modelbox, p, orb, die_on_error)
            logger.info("\n Nadir tracks have been written in "
                        "{}".format(p.working_directory))
        else:
            build_swath.orbit2swath(modelbox, p, orb, die_on_error)
            logger.info("\n SWOT Grids and nadir tracks have been written in "
                        "{}".format(p.working_directory))
        logger.info("-----------------------------------------------")

    # - Initialize random coefficients that are used to compute
    #   random errors following the specified spectrum
    err, errnad = mod.load_error(p, nadir_alone=nadir_alone)

    # - Compute interpolated SSH and errors for each pass, at each
    #   cycle
    logger.info('Compute interpolated SSH and errors:')
    #   load all SWOT grid files (one for each pass)
    listsgridfile = sorted(glob.glob(p.filesgrid + '_p*.nc'))
    if not listsgridfile:
        logger.error('\n There is no SWOT grid file in {}, run simulator with'
                     ' option makesgrid set to true in your params'
                     ' file'.format(p.working_directory))
        sys.exit(1)
    # Build model time steps from parameter file
    modeltime = numpy.arange(0, p.nstep*p.timestep, p.timestep)
    #   Remove the grid from the list of model files
    if p.file_input and p.file_grid_model is None:
        logger.info("WARNING: the first file is not used to build data")
        list_file.remove(list_file[0])
        if len(modeltime) > len(list_file):
            logger.error('There is not enough model files in the list of'
                         'files')
            sys.exit(1)
    # - Loop on SWOT grid files to construct a list of jobs
    jobs = []
    p2 = mod_tools.todict(p)
    for sgridfile in listsgridfile:
        jobs.append([sgridfile, p2, listsgridfile, list_file, modelbox,
                     model_data, modeltime, err, errnad])
    ok = False
    # - Process list of jobs using multiprocessing
    try:
        ok = make_swot_data(p.proc_count, jobs, die_on_error, p.progress_bar)
    except parallel.DyingOnError:
        logger.error('An error occurred and all errors are fatal')
        sys.exit(1)
    # - Write Selected parameters in a txt file
    timestop = datetime.datetime.now()
    timestop = timestop.strftime('%Y%m%dT%H%M%SZ')
    timestart = timestart.strftime('%Y%m%dT%H%M%SZ')
    op_file = 'swot_simulator_{}_{}.output'.format(timestart, timestop)
    op_file = os.path.join(p.working_directory, op_file)
    rw_data.write_params(p, op_file)
    if ok is True:
        if p.progress_bar is True:
            __ = mod_tools.update_progress(1, 'All passes have been processed',
                                           '')
        else:
            __ = logger.info('All passes have been processed')
        logger.info("\n Simulated swot files have been written in {}".format(
                     p.working_directory))
        logger.info(''.join(['-'] * 61))
        #"----------------------------------------------------------")
        sys.exit(0)
    logger.error('\nERROR: At least one of the outputs was not saved.')
    sys.exit(1)


def exc_formatter(exc):
    """Format exception returned by sys.exc_info() as a string so that it can
    be serialized by pickle and stored in the JobsManager."""
    error_msg = traceback.format_exception(exc[0], exc[1], exc[2])
    return error_msg


def err_formatter(pid, grid, cycle, exc):
    """Transform errors stored by the JobsManager into readable messages."""
    msg = None
    if cycle < 0:
        msg = '/!\ Error occurred while processing grid {}'.format(grid)
    else:
        _msg = '/!\ Error occurred while processing cycle {} on grid {}'
        msg = _msg.format(cycle, grid)
    return msg


def run_nadir(p, die_on_error=False):

    # - Initialize some parameters values
    timestart = datetime.datetime.now()
    mod_tools.initialize_parameters(p)
    p.nadir = True
    p.karin = False
    p.phase = False
    p.roll = False
    p.baseline_dilation = False
    p.timing = False
    p.halfswath = 60.

    # - Progress bar variables are global
    global ntot

    # Build model time steps from parameter file
    modeltime = numpy.arange(0, p.nstep*p.timestep, p.timestep)
    # - Read list of user model files """
    if p.file_input is not None:
        list_file = [line.strip() for line in open(p.file_input)]
        if len(modeltime) > len(list_file):
            logger.error('There is not enough model files in the list of '
                         'files')
            sys.exit(1)
    else:
        list_file = None

    # ############################################
    # Select the spatial and temporal domains
    # ############################################

    # - Read model input coordinates '''
    model_data, list_file = mod.load_coordinate_model(p)
    # if no modelbox is specified (modelbox=None), the domain of the input data
    # is taken as a modelbox
    # coordinates from the region defined by modelbox are selected
    logger.debug('Read input')
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
            logger.error('modelbox should be provided if no model file is '
                         'provided')
            sys.exit()
    p.modelbox_calc = modelbox
    logger.debug(p.file_input)
    if p.file_input is not None:
        model_data.read_coordinates()
        # Select model data in the region modelbox
        model_data.len_coord = len(numpy.shape(model_data.vlon))
        if p.grid == 'regular' or model_data.len_coord == 1:
            if modelbox[0] < modelbox[1]:
                _i_lon = numpy.where(((modelbox[0] - 1) <= model_data.vlon)
                                     & (model_data.vlon <= (modelbox[1]+1)))[0]
            else:
                _i_lon = numpy.where(((modelbox[0]-1) <= model_data.vlon)
                                     | (model_data.vlon <= (modelbox[1]+1)))[0]
            model_data.model_index_lon = _i_lon
            _i_lat = numpy.where(((modelbox[2] - 1) <= model_data.vlat)
                                   & (model_data.vlat <= (modelbox[3]+1)))[0]

            model_data.model_index_lat = _i_lat
            model_data.vlon = model_data.vlon[model_data.model_index_lon]
            model_data.vlat = model_data.vlat[model_data.model_index_lat]

        else:
            if modelbox[0] < modelbox[1]:
                _i_box = numpy.where(((modelbox[0] - 1) <= model_data.vlon)
                                     & (model_data.vlon <= (modelbox[1]+1))
                                     & ((modelbox[2]-1) <= model_data.vlat)
                                     & (model_data.vlat <= (modelbox[3]+1)))
            else:
                _i_box = numpy.where(((modelbox[0]-1) <= model_data.vlon)
                                     | (model_data.vlon <= (modelbox[1]+1))
                                     & ((modelbox[2]-1) <= model_data.vlat)
                                     & (model_data.vlat <= (modelbox[3]+1)))

            model_data.model_index = _i_box
            model_data.vlon = model_data.vlon[model_data.model_index]
            model_data.vlat = model_data.vlat[model_data.model_index]
        model_data.model = p.model
        model_data.vloncirc = numpy.rad2deg(numpy.unwrap(model_data.vlon))
    # Ugly trick when model box is [0 360] to avoid box being empty (360=0%360)
    if modelbox[1] == 0:
        modelbox[1] = 359.99

    # - Initialize random coefficients that are used to compute
    #   random errors following the specified spectrum
    err, errnad = mod.load_error(p, nadir_alone=nadir_alone)

    # - Compute interpolated SSH and errors for each pass, at each
    #   cycle
    logger.info('Compute interpolated SSH and errors:')
    #   Remove the grid from the list of model files
    if p.file_input:
        list_file.remove(list_file[0])
        if len(modeltime) > len(list_file):
            logger.error('There is not enough model files in the list of'
                         ' files')
            sys.exit(1)
    #   Initialize progress bar variables
    ntot = 1

    #   Initialize list of satellites
    if not isinstance(p.ephemeris, list):
        p.ephemeris = [p.ephemeris]
    for filesat in p.ephemeris:
        # Select satellite
        # ntmp, nfilesat = os.path.split(filesat[istring:-4])
        nfilesat = os.path.basename(os.path.splitext(filesat)[0])
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
            gridfile = '{}{}_grid.nc'.format((p.filesgrid).strip(),
                                              nfilesat.strip())
            ngrid = rw_data.Sat_nadir(nfile=gridfile)
            ngrid.file = gridfile
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
        # if p.file_input is not None:
        #    _ind = numpy.where((numpy.min(ngrid.lon) <= model_data.vlon)
        #                       & (model_data.vlon <= numpy.max(ngrid.lon))
        #                       & (numpy.min(ngrid.lat) <= model_data.vlat)
        #                       & (model_data.vlat <= numpy.max(ngrid.lat)))
        #    model_index = _ind
        # - Generate and nadir-like data:
        #   Compute number of cycles needed to cover all nstep model timesteps
        rcycle = (p.timestep * p.nstep)/float(ngrid.cycle)
        ncycle = int(rcycle)
        #   Loop on all cycles
        for cycle in range(0, ncycle + 1):
            ### TODO move this line somwhere where we have ifile information
            #if ifile > (p.nstep/p.timestep + 1):
            #    break
            #   Create SWOT-like and Nadir-like data
            if p.file_input is None:
                model_data = []
            logger.debug('compute SSH nadir')
            nfile = numpy.shape(p.ephemeris)[0]*rcycle
            create = mod.create_Nadirlikedata(cycle, nfile, list_file,
                                              modelbox, ngrid, model_data,
                                              modeltime, errnad, p,
                                              progress_bar=True)
            SSH_true_nadir, vindice, time, progress = create
            # SSH_true_nadir, vindice_nadir=create_Nadirlikedata(cycle, sgrid,
            # ngrid, model_data, modeltime, err, errnad, p)
            #   Save outputs in a netcdf file
            ngrid.gridfile = filesat
            if (~numpy.isnan(vindice)).any() or p.file_input is None:
                err = errnad
                err.wtnadir = numpy.zeros((1))
                err.wet_tropo2nadir = numpy.zeros((1))
                logger.debug('write file')
                mod.save_Nadir(cycle, ngrid, errnad, err, p, time=time,
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
        str1 = 'All passes have been processed'
        progress = mod_tools.update_progress(1, str1, '')
    # - Write Selected parameters in a txt file
    timestop = datetime.datetime.now()
    timestop = timestop.strftime('%Y%m%dT%H%M%SZ')
    timestart = timestart.strftime('%Y%m%dT%H%M%SZ')
    op_file = 'nadir_simulator_{}_{}.output'.format(timestart, timestop)
    op_file = os.path.join(p.working_directory, op_file)
    rw_data.write_params(p, op_file)
    logger.info("\nSimulated orbit files have been written in {}".format(
                p.working_directory))
    logger.info("----------------------------------------------------------")


def make_swot_data(_proc_count, jobs, die_on_error, progress_bar):
    """ Compute SWOT-like data for all grids and all cycle, """
    # - Set up parallelisation parameters
    proc_count = min(len(jobs), _proc_count)

    status_updater = mod_tools.update_progress_multiproc
    jobs_manager = parallel.JobsManager(proc_count, status_updater,
                                        exc_formatter, err_formatter)

    ok = jobs_manager.submit_jobs(worker_method_swot, jobs, die_on_error,
                                  progress_bar)

    if not ok:
        # Display errors once the processing is done
        jobs_manager.show_errors()

    return ok


def worker_method_swot(*args, **kwargs):
    msg_queue, sgridfile, p2, listsgridfile = args[:4]

    list_file, modelbox, model_data, modeltime, err, errnad = args[4:]
    p = mod_tools.fromdict(p2)
    if err is None:
        nadir_alone = True
    else:
        nadir_alone = False
    compute_nadir = (('Altimeter' in p.noise) or (nadir_alone is True))
    #   Load SWOT grid files (Swath and nadir)
    if compute_nadir is True:
        ngrid = mod.load_ngrid(sgridfile, p, nadir_alone=nadir_alone)
        ngrid.gridfile = sgridfile
    else:
        ngrid = None
    if nadir_alone is False:
        sgrid = mod.load_sgrid(sgridfile, p)
        sgrid.gridfile = sgridfile
    else:
        sgrid = ngrid
    # Set Teval and nTeval to None to interpolate the mask once
    Teval = None
    nTeval = None
    #   Select model data around the swath to reduce interpolation cost in
    #   griddata

    # - Generate SWOT like and nadir-like data:
    #   Compute number of cycles needed to cover all nstep model timesteps
    rcycle = (p.timestep * p.nstep)/float(sgrid.cycle)
    ncycle = int(rcycle)
    #   Loop on all cycles
    for cycle in range(0, ncycle+1):
        # #TODO move this somwhere where we have ifile information
        #if ifile > (p.nstep/p.timestep + 1):
        #    break
        # Add a message to tell the main program that the cycle is being
        # processed
        msg_queue.put((os.getpid(), sgridfile, cycle + 1, None))

        #   Process_cycle: create SWOT-like and Nadir-like data
        if not p.file_input:
            model_data = []
        create = mod.create_SWOTlikedata(cycle, list_file, modelbox,
                                         sgrid, ngrid, model_data, modeltime,
                                         err, errnad, p, Teval=Teval,
                                         nTeval=nTeval)
        out_var, time, Teval, nTeval = create
        #   Save outputs in a netcdf file
        if nadir_alone is True:
            out_var['vindice'] = + out_var['vindice_nadir']
        if (~numpy.isnan(out_var['vindice'])).any() or not p.file_input:
            if nadir_alone is False:
                mod.save_SWOT(cycle, sgrid, err, p, out_var, time=time,
                              save_var=p.product_type)
            if compute_nadir is True:
                mod.save_Nadir(cycle, ngrid, errnad, err, p, out_var,
                               time=time)
    # Add a special message once a grid has been completely processed
    msg_queue.put((os.getpid(), sgridfile, None, None))
