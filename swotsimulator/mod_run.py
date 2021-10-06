''' Module to create SWOT data:

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
from scipy import interpolate
import numpy
import math
import glob
import sys
import time
import datetime
import logging
import swotsimulator.rw_data as rw_data
import swotsimulator.build_error as build_error
import swotsimulator.mod_tools as mod_tools
import swotsimulator.const as const
import multiprocessing
# Define logger level for debug purposes
logger = logging.getLogger(__name__)


def load_error(p, nadir_alone=False, seed=0):
    '''Initialize random coefficients that are used to compute
    random errors following the specified spectrum. \n
    If a random coefficient file is specified, random coefficients
    are ldir_alone=Falseoaded from this file.
    '''
    if ('Altimeter' in p.noise) or (nadir_alone is True):
        errnad = build_error.errornadir(p, nadir_alone=nadir_alone)
    else:
        errnad = None
    if nadir_alone is True:
        err = None
        errnad.init_error(p, seed=seed, nadir_alone=True)
    else:
        err = build_error.error(p)
        try:
            nhalfswath = int((p.half_swath - p.half_gap) / p.delta_ac) + 1
        except AttributeError:
            nhalfswath = 60.
        err.init_error(p, 2*nhalfswath, seed=seed)
        if 'Altimeter' in p.noise:
            errnad.init_error(p, seed=seed)
    return err, errnad


def load_sgrid(sgridfile, p):
    '''Load SWOT swath and Nadir data for file sgridfile '''

    # Load SWOT swath file
    sgrid = rw_data.Sat_SWOT(nfile=sgridfile)
    cycle = 0
    x_al = []
    x_ac = []
    al_cycle = 0
    timeshift = 0
    sgrid.load_swath(cycle=cycle, x_al=x_al, x_ac=x_ac, al_cycle=al_cycle,
                     timeshift=timeshift)
    sgrid.loncirc = numpy.rad2deg(numpy.unwrap(sgrid.lon))
    # Extract the pass number from the file name
    ipass = int(sgridfile[-6: -3])
    sgrid.ipass = ipass
    return sgrid


def load_coordinate_model(p):
    model = p.model
    if p.file_input is not None:
        list_file = [line.strip() for line in open(p.file_input)]
    else:
        list_file = None
    # - Read model input coordinates '''
    # If a list of model files are specified, read model file coordinates
    if p.file_input is not None:
        model_data_ctor = getattr(rw_data, model)
        if p.file_grid_model is not None:
            #_filename = list(p.file_grid_model)
            _filename = p.file_grid_model
        else:
            logger.info("WARNING: First file of list of files is used for"
                        "coordinates only")
            _filename = os.path.join(p.indatadir, list_file[0]) #.split(',')
            list_file.remove(list_file[0])
        filename = _filename
        model_data = model_data_ctor(p, nfile=filename, lon=p.lon, lat=p.lat,
                                     var=p.var)
    return model_data, list_file


def load_ngrid(sgridfile, p, nadir_alone=False):
    ipass = int(sgridfile[-6: -3])
    # Load Nadir track file
    if nadir_alone is True:
        nfile = '{:s}_p{:04d}.nc'.format((p.filesgrid).strip(), ipass)
    else:
        nfile = '{:s}nadir_p{:04d}.nc'.format((p.filesgrid).strip(), ipass)
    ngrid = rw_data.Sat_nadir(nfile=nfile)
    cycle = 0
    x_al = []
    al_cycle = 0
    timeshift = 0
    ngrid.load_orb(cycle=cycle, x_al=x_al, al_cycle=al_cycle,
                   timeshift=timeshift)
    ngrid.loncirc = numpy.rad2deg(numpy.unwrap(ngrid.lon))
    ngrid.ipass = ipass
    return ngrid


def interpolate_regular_1D(p, lon_in, lat_in, var, lon_out, lat_out,
                           Teval=None):
    ''' Interpolation of data when grid is regular and coordinate in 1D. '''
    #lon_in = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(lon_in)))
    interp = interpolate.RectBivariateSpline
    if Teval is None or p.ice_mask is True:
        _Teval = interp(lat_in, lon_in, numpy.isnan(var), kx=1, ky=1, s=0)
        Teval = _Teval.ev(lat_out, lon_out)
    # Trick to avoid nan in interpolation
    var_mask = + var
    var_mask._sharedmask=False
    var_mask[numpy.isnan(var_mask)] = 0.
    # Interpolate variable
    _var = interp(lat_in, lon_in, var_mask, kx=1, ky=1, s=0)
    var_out = _var.ev(lat_out, lon_out)
    # Mask variable with Teval
    var_out[Teval > 0] = numpy.nan
    return var_out, Teval


def interpolate_irregular_pyresample(swath_in, var, grid_out, radius,
                                     interp_type='nearest'):
    ''' Interpolation of data when grid is irregular and pyresample is
    installed.'''
    import pyresample as pr
    if interp_type == 'nearest':
        interp = pr.kd_tree.resample_nearest
        radius_n = radius * 10**3
        var_out = interp(swath_in, var, grid_out,
                         radius_of_influence=radius_n, epsilon=100)
    else:
        interp = pr.kd_tree.resample_gauss
        radius_g = radius * 3 * 10**3
        sigma_g = radius * 10**3
        var_out = interp(swath_in, var, grid_out, radius_of_influence=radius_g,
                         sigmas=sigma_g, fill_value=0)
    var_out[var_out == 0] = numpy.nan
    return var_out


def select_modelbox(sgrid, model_data, p):
    # mask=numpy.zeros((numpy.shape(model_data.vlon)))
    # nal=len(sgrid.x_al)
    # for kk in range(0,nal,10):
    # dlon1=model_data.vlon-sgrid.lon_nadir[kk]
    # dlon=numpy.minimum(numpy.mod(dlon1,360),numpy.mod(-dlon1,360))
    # ddist=(((dlon)*numpy.cos(sgrid.lat_nadir[kk]*numpy.pi/180.)*
    # const.deg2km)**2+((model_data.vlat-sgrid.lat_nadir[kk])
    # *const.deg2km)**2 )
    #  mask[ddist<const.radius_interp**2]=1
    model_data.len_coord = len(numpy.shape(model_data.vlon))
    if p.grid == 'regular' or model_data.len_coord == 1:
        _ind_lon = numpy.where((numpy.min(sgrid.lon) <= model_data.vlon)
                               & (model_data.vlon <= numpy.max(sgrid.lon)))
        model_data.lon1d = model_data.vlon[_ind_lon]
        _ind_lat = numpy.where((numpy.min(sgrid.lat) <= model_data.vlat)
                               & (model_data.vlat <= numpy.max(sgrid.lat)))
        model_data.lat1d = model_data.vlat[_ind_lat]
    else:
        model_index = numpy.where((numpy.min(sgrid.lon) <= model_data.vlon)
                                  & (model_data.vlon <= numpy.max(sgrid.lon))
                                  & (numpy.min(sgrid.lat) <= model_data.vlat)
                                  & (model_data.vlat <= numpy.max(sgrid.lat)))
        model_data.lon1d = model_data.vlon[model_index].ravel()
        model_data.lat1d = model_data.vlat[model_index].ravel()
        model_data.vlon = model_data.vlon[model_index]
        model_data.vlat = model_data.vlat[model_index]

    # nx = len(lon1)
    # ny = len(lat1)
    return None  # model_index


def create_SWOTlikedata(cycle, list_file, modelbox, sgrid, ngrid,
                        model_data, modeltime, err, errnad, p,
                        Teval=None, nTeval=None):
    '''Create SWOT and nadir errors err and errnad, interpolate model SSH model
    _data on swath and nadir track, compute SWOT-like and nadir-like data
    for cycle, SWOT grid sgrid and ngrid. '''
    #   Initialiaze errors and SSH
    if err is None:
        nadir_alone = True
    else:
        nadir_alone = False
    compute_nadir = (('Altimeter') in p.noise or (nadir_alone is True))
    if nadir_alone is False:
        shape_sgrid = (numpy.shape(sgrid.lon)[0], numpy.shape(sgrid.lon)[1])
        err.karin = numpy.zeros(shape_sgrid)
        err.roll = numpy.zeros(shape_sgrid)
        err.phase = numpy.zeros(shape_sgrid)
        err.corrected_roll_phase = numpy.zeros(shape_sgrid)
        err.baseline_dilation = numpy.zeros(shape_sgrid)
        err.timing = numpy.zeros(shape_sgrid)
        err.wet_tropo1 = numpy.zeros(shape_sgrid)
        err.wet_tropo2 = numpy.zeros(shape_sgrid)
        err.ssb = numpy.zeros(shape_sgrid)
        err.wt = numpy.zeros(shape_sgrid)
    if compute_nadir is True:
        shape_ngrid = (numpy.shape(ngrid.lon)[0])
    out_var = {}
    for key in model_data.input_var_list.keys():
        if nadir_alone is False:
            out_var[key] = numpy.zeros(shape_sgrid)
        nkey = '{}_nadir'.format(key)
        # Initialize nadir variable in case the simulator is run with nadir
        # option at False
        out_var[nkey] = 0
        if compute_nadir:
            out_var[nkey] = numpy.zeros(shape_ngrid)
    if nadir_alone is False:
        out_var['mask_land'] = numpy.zeros(shape_sgrid)
        out_var['vindice'] = numpy.full(shape_sgrid, numpy.nan)
        out_var['vindice_nadir'] = 0
    if compute_nadir:
        if nadir_alone is False:
            err.wet_tropo1nadir = numpy.zeros(shape_ngrid)
            err.wet_tropo2nadir = numpy.zeros(shape_ngrid)
            err.wtnadir = numpy.zeros(shape_ngrid)
        errnad.nadir = numpy.zeros(shape_ngrid)
        errnad.wet_tropo1nadir = numpy.zeros(shape_ngrid)
        out_var['vindice_nadir'] = numpy.full(shape_ngrid, numpy.nan)
    # Definition of the time in the model
    date1 = cycle * sgrid.cycle
    time = sgrid.time + date1
    # Look for satellite data that are beween step-p.timestep/2 and
    # step+p.step/2
    if p.file_input is not None:
        time_shift_end = time[-1] - sgrid.timeshift
        time_shift_start = time[0] - sgrid.timeshift
        model_tmin = modeltime - p.timestep/2.
        model_tmax = modeltime + p.timestep/2.
        index_filemodel = numpy.where((time_shift_end >= model_tmin)
                                      & (time_shift_start < model_tmax))
        # At each step, look for the corresponding time in the satellite data
        for ifile in index_filemodel[0]:
            # If there are satellite data, Get true SSH from model
            if numpy.shape(index_filemodel)[1] > 0:
                # if numpy.shape(index)[1]>1:
                # Select part of the track that corresponds to the time of the
                # model (+-timestep/2)
                model_tmin = modeltime[ifile] - p.timestep/2.
                model_tmax = modeltime[ifile] + p.timestep/2.
                if nadir_alone is False:
                    time_shift = time - sgrid.timeshift
                    ind_time = numpy.where((time_shift >= model_tmin)
                                           & (time_shift < model_tmax))
                if compute_nadir is True:
                    time_shift = time - ngrid.timeshift
                    ind_nadir_time = numpy.where((time_shift >= model_tmin)
                                                 & (time_shift < model_tmax))
                # Handle files with multiple time dimensions
                infile = int(ifile /p.dim_time)
                filetime = ifile - infile * p.dim_time

                # Load data from this model file
                model_step_ctor = getattr(rw_data, model_data.model)
                nfile = os.path.join(p.indatadir, list_file[infile])
                model_step = model_step_ctor(p, nfile=nfile,
                                             time=filetime)
                input_var = {}
                if p.grid == 'regular' or model_data.len_coord ==1:
                    model_step.read_var()
                    for key in model_step.input_var.keys():
                        _indlat = model_data.model_index_lat
                        _tmp = model_step.input_var[key][_indlat, :]
                        input_var[key] = +_tmp[:, model_data.model_index_lon]
                else:
                    model_step.read_var(index=model_data.model_index)
                    for key in model_step.input_var.keys():
                        input_var[key] = + model_step.input_var[key]
                # - Interpolate Model data on a SWOT grid and/or along the
                #   nadir track
                # Handle Greenwich line
                _greenwich = False
                if nadir_alone is False:
                    lon_grid = + sgrid.lon
                    if numpy.max(lon_grid) > 359:
                        _greenwich = True
                        lon_grid = numpy.mod(lon_grid + 180., 360.) - 180.
                if compute_nadir is True:
                    lon_ngrid = + ngrid.lon
                    if numpy.max(lon_ngrid) > 359:
                        _greenwich = True
                        lon_ngrid = numpy.mod(lon_ngrid + 180., 360.) - 180.
                lon_model = + model_data.vlon
                if _greenwich is True:
                    lon_model = numpy.mod(lon_model + 180., 360.) - 180.
                # if grid is regular, use interpolate.RectBivariateSpline to
                # interpolate

                if p.grid == 'regular' or model_data.len_coord == 1:
                    # ########################TODO
                    # To be moved to routine rw_data
                    indsorted = numpy.argsort(lon_model)
                    lon_model = lon_model[indsorted]
                    for key in input_var.keys():
                        input_var[key] = input_var[key][:, indsorted]
                    # Flatten satellite grid and select part of the track
                    # corresponding to the model time
                    interp = interpolate_regular_1D
                    if nadir_alone is False:
                        lonswot = lon_grid[ind_time[0], :].flatten()
                        latswot = sgrid.lat[ind_time[0], :].flatten()
                        for key in input_var.keys():
                            _ssh, Teval = interp(p, lon_model, model_data.vlat,
                                                 input_var[key], lonswot,
                                                 latswot, Teval=Teval)
                            nal, nac = numpy.shape(sgrid.lon[ind_time[0], :])
                            out_var[key][ind_time[0], :] = _ssh.reshape(nal, nac)
                            if key == 'ssh_true':
                                _mask_land = numpy.zeros((nal, nac))
                                Teval2d = Teval.reshape(nal, nac)
                                _mask_land[Teval2d > 0] = 3
                                out_var['mask_land'][ind_time[0], :] = _mask_land
                    if compute_nadir is True:
                        lonnadir = lon_ngrid[ind_nadir_time[0]].ravel()
                        latnadir = ngrid.lat[ind_nadir_time[0]].ravel()
                        for key in input_var.keys():
                            _ssh, nTeval = interp(p, lon_model,
                                                  model_data.vlat,
                                                  input_var[key], lonnadir,
                                                  latnadir, Teval=nTeval)
                            nkey = '{}_nadir'.format(key)
                            out_var[nkey][ind_nadir_time[0]] = _ssh
                else:
                    # Grid is irregular, interpolation can be done using
                    # pyresample module if it is installed or griddata
                    # function from scipy.
                    # Note that griddata is slower than pyresample functions.
                    try:
                        import pyresample as pr
                        wrap_lon = pr.utils.wrap_longitudes
                        geomdef = pr.geometry.SwathDefinition
                        interp = interpolate_irregular_pyresample
                        lon_model = wrap_lon(lon_model)
                        if model_data.len_coord <= 1:
                            logger.error('Model grid is irregular,'
                                         'coordinates should be in 2d')
                            sys.exit(1)
                        swath_def = geomdef(lons=lon_model,
                                            lats=model_data.vlat)
                        if nadir_alone is False:
                           lon_grid = wrap_lon(lon_grid)
                           grid_def = geomdef(lons=lon_grid[ind_time[0], :],
                                              lats=sgrid.lat[ind_time[0], :])
                           for key in input_var.keys():
                                _ssh = interp(swath_def, input_var[key], grid_def,
                                              max(p.delta_al, p.delta_ac),
                                              interp_type=p.interpolation)
                                out_var[key][ind_time[0], :] = _ssh
                                if key == 'ssh_true':
                                    out_var['mask_land'][numpy.isnan(_ssh)] = 3
                        if compute_nadir is True:
                            lon_ngrid = wrap_lon(lon_ngrid)
                            ngrid_def = geomdef(lons=lon_ngrid[ind_nadir_time[0]],
                                                lats=ngrid.lat[ind_nadir_time[0]])
                            for key in input_var.keys():
                                _ssh = interp(swath_def, input_var[key],
                                              ngrid_def,
                                              p.delta_al,
                                              interp_type=p.interpolation)
                                nkey = '{}_nadir'.format(key)
                                out_var[nkey][ind_nadir_time[0]] = _ssh
                    except ImportError:
                        interp = interpolate.griddata
                        model_ravel = (lon_model.ravel(),
                                       model_data.vlat.ravel())
                        if nadir_alone is False:
                            for key in input_var.keys():
                                _ssh = interp(model_ravel, input_var[key].ravel(),
                                              (lon_grid[ind_time[0], :],
                                              sgrid.lat[ind_time[0], :]),
                                              method=p.interpolation)
                                if key == 'ssh_true':
                                    out_var['mask_land'][numpy.isnan(_ssh)] = 3
                                out_var[key][ind_time[0], :] = _ssh
                        if compute_nadir is True:
                            for key in input_var.keys():
                                _ssh = interp(model_ravel,
                                              input_var[key].ravel(),
                                              (lon_ngrid[ind_nadir_time[0]],
                                              ngrid.lat[ind_nadir_time[0]]),
                                              method=p.interpolation)
                                nkey = '{}_nadir'.format(key)
                                out_var[nkey][ind_nadir_time[0]] = _ssh
                        if p.interpolation == 'nearest':
                            if modelbox[0] > modelbox[1]:
                                if nadir_alone is False:
                                    ind = numpy.where(((sgrid.lon < modelbox[0])
                                                      & (sgrid.lon > modelbox[1]))
                                                      | (sgrid.lat < modelbox[2])
                                                      | (sgrid.lat > modelbox[3]))
                                    for key in out_var.keys():
                                        if 'nadir' not in key:
                                            out_var[key][ind] = numpy.nan
                                if compute_nadir is True:
                                    lontmp = ngrid.lon
                                    lattmp = ngrid.lat
                                    ind = numpy.where(((lontmp < modelbox[0])
                                                      & (lontmp > modelbox[1]))
                                                      | (lattmp < modelbox[2])
                                                      | (lattmp > modelbox[3]))
                                    del lontmp, lattmp
                                    for key in out_var.keys():
                                        if 'nadir' in key:
                                            out_var[key][ind] = numpy.nan
                            else:
                                if nadir_alone is False:
                                    ind = numpy.where((sgrid.lon < modelbox[0])
                                                      | (sgrid.lon > modelbox[1])
                                                      | (sgrid.lat < modelbox[2])
                                                      | (sgrid.lat > modelbox[3]))
                                    for key in out_var.keys():
                                        if 'nadir' not in key:
                                            out_var[key][ind] = numpy.nan
                                if compute_nadir is True:
                                    lontmp = ngrid.lon
                                    lattmp = ngrid.lat
                                    ind = numpy.where((lontmp < modelbox[0])
                                                      | (lontmp > modelbox[1])
                                                      | (lattmp < modelbox[2])
                                                      | (lattmp > modelbox[3]))
                                    del lontmp, lattmp
                                    for key in out_var.keys():
                                        if 'nadir' in key:
                                            out_var[key][ind] = numpy.nan
                if nadir_alone is False:
                    out_var['vindice'][ind_time[0], :] = ifile
                if compute_nadir is True:
                    out_var['vindice_nadir'][ind_nadir_time[0]] = ifile
                    del input_var, model_step, ind_nadir_time
                else:
                    del ind_time, input_var, model_step
    if nadir_alone is False:
        err.make_error(sgrid, cycle, out_var['ssh_true'], p)
        if p.product_type != 'expert':
            err.reconstruct_2D(p, sgrid.x_ac)
            err.make_SSH_error(out_var['ssh_true'], p)
    if compute_nadir is True:
        errnad.make_error(ngrid, cycle, out_var['ssh_true_nadir'], p)
        if nadir_alone is False:
            if p.nbeam == 1:
                 errnad.SSH = (out_var['ssh_true_nadir'] + errnad.nadir
                              + err.wet_tropo1nadir)
            else:
                errnad.SSH = (out_var['ssh_true_nadir'] + errnad.nadir
                             + err.wet_tropo2nadir)
        else:
            errnad.SSH = (out_var['ssh_true_nadir'] + errnad.nadir
                          + errnad.wet_tropo1nadir)
    return out_var, time, Teval, nTeval


def make_flags(var, sgrid, modelbox, beam_pos):
    sgrid.lon = numpy.mod(sgrid.lon + 360, 360)
    if modelbox is not None:
        ind = numpy.where((sgrid.lon < modelbox[0]) | (sgrid.lon > modelbox[1])
                          | (sgrid.lat < modelbox[2])
                          | (sgrid.lat > modelbox[3]))
        var[ind] = 0
    flag_k = numpy.zeros(numpy.shape(var))
    flag_k[numpy.isnan(var)] = 3
    var_r = numpy.zeros((numpy.shape(var)[0], 2))
    var_r[:, 0] = var[:, beam_pos[0]]
    var_r[:, 1] = var[:, beam_pos[1]]
    flag_r = numpy.zeros(numpy.shape(var_r))
    flag_r[numpy.isnan(var_r)] = 1

    return flag_k

def save_SWOT(cycle, sgrid, err, p, out_var, time=[],
              save_var='all'):
    ofile = '{}_c{:03d}_p{:04d}.nc'.format(p.file_output, cycle + 1,
                                           sgrid.ipass)
    OutputSWOT = rw_data.Sat_SWOT(nfile=ofile, lon=(sgrid.lon+360) % 360,
                                  lat=sgrid.lat, time=time, x_ac=sgrid.x_ac,
                                  x_al=sgrid.x_al, cycle=sgrid.cycle,
                                  lon_nadir=(sgrid.lon_nadir+360) % 360,
                                  lat_nadir=sgrid.lat_nadir)
    OutputSWOT.gridfile = sgrid.gridfile
    OutputSWOT.first_time = p.first_time
    if p.product_type != 'expert':
        err.SSH2 = err.SSH
    if save_var == 'mockup':
        all_var = make_empty_vars(sgrid)
        beam_pos = (12, 47)
        tmp_var = + err.SSH
        flag = make_flags(tmp_var, sgrid, p.modelbox, beam_pos)
        all_var['karin_surf_type'] = out_var['mask_land']
        all_var['rad_surf_type'][:, 0] = out_var['mask_land'][:, beam_pos[0]]
        all_var['rad_surf_type'][:, 1] = out_var['mask_land'][:, beam_pos[1]]
        all_var['ssh_karin_swath'] = + err.SSH
        all_var['ssha_uncert'] = err.karin
    else:
        all_var = None
        err.SSH2 = None
        flag = None
    if save_var == 'expert':
        OutputSWOT.write_data(out_var,
                              roll_err_1d=err.roll1d, phase_err_1d=err.phase1d,
                              corrected_rollphase_err_1d = err.rollphase_est1d,
                              bd_err_1d=err.baseline_dilation1d,
                              ssb_err=err.ssb, karin_err=err.karin,
                              pd_err_1b=err.wet_tropo1,
                              pd_err_2b=err.wet_tropo2, pd=err.wt,
                              timing_err_1d=err.timing1d)
    elif save_var == 'mockup':
        OutputSWOT.write_data(empty_var=all_var)
    else:
        OutputSWOT.write_data(out_var, roll_err=err.roll,
                              bd_err=err.baseline_dilation,
                              corrected_rollphase_err=err.corrected_roll_phase,
                              phase_err=err.phase, ssb_err=err.ssb,
                              karin_err=err.karin, pd_err_1b=err.wet_tropo1,
                              pd_err_2b=err.wet_tropo2, pd=err.wt,
                              timing_err=err.timing, ssh_obs=err.SSH,
                              )
    return None


def save_Nadir(cycle, ngrid, errnad, err, p, out_var, time=[]):
    if type(ngrid.ipass) == str:
        ofile = '{}nadir_c{:03d}_{}.nc'.format(p.file_output, cycle + 1,
                                               ngrid.ipass)
    else:
        ofile = '{}nadir_c{:03d}_p{:04d}.nc'.format(p.file_output, cycle + 1,
                                                    ngrid.ipass)
    OutputNadir = rw_data.Sat_nadir(nfile=ofile,
                                    lon=(ngrid.lon+360) % 360,
                                    lat=ngrid.lat, time=time, x_al=ngrid.x_al,
                                    cycle=ngrid.cycle)
    OutputNadir.first_time = p.first_time
    OutputNadir.gridfile = ngrid.gridfile
    if err is None:
        OutputNadir.write_data(out_var,
                               nadir_err=errnad.nadir, SSH_obs=errnad.SSH)

    else:
        OutputNadir.write_data(out_var,
                               nadir_err=errnad.nadir, SSH_obs=errnad.SSH,
                               pd_err_1b=err.wet_tropo1nadir, pd=err.wtnadir,
                               pd_err_2b=err.wet_tropo2nadir)
    return None


def make_empty_vars(sgrid):
    nal, nac = numpy.shape(sgrid.lon)
    var = {}
    for key in const.list_var_mockup:
        var[key] = numpy.full((nal, nac), numpy.nan)
        if key == 'rad_surf_type':
            var[key] = numpy.full((nal, 2), numpy.nan)
    return var
