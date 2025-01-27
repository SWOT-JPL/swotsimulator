'''Modules to run the main program:
Contains the following functions:
load_error(param_file) : load or compute random coefficients to compute errors
\n
load_sgrid(sgridfile, p): load or compute SWOT grid
\n
load_ngrid(sgridfile, p): load or compute nadir grid
\n
select_modelbox(sgrid, model_data): Select region of interest with params file
or model data
\n
create_SWOTlikedata(cycle, ntotfile, list_file, sgrid, ngrid, model_data,
modeltime, err, errnad, p): interpolate data on SWOT swath
\n
CreateNadirlikedata(cycle, ntotfile, list_file, ngrid, model_data, modeltime,
errnad, p, progress_bar=False): Interpolate data on nadir track
\n
def save_SWOT(cycle, sgrid, err, p, time=[], vindice=[], SSH_true=[]):
Save SWOT-like data
\n
save_Nadir(cycle, ngrid, errnad, err , p, time=[], vindice_nadir=[],
SSH_true_nadir=[]): Save Nadir-like data
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
from scipy import interpolate
import numpy
import swotsimulator.rw_data as rw_data
istep = 0
ntot = 1
ifile = 0


def load_error(p):
    '''Initialize random coefficients that are used to compute
    random errors following the specified spectrum. \n
    If a random coefficient file is specified, random coefficients
    are loaded from this file.
    '''
    import swotsimulator.build_error as build_error
    err = build_error.error(p)
    if p.nadir:
        errnad = build_error.errornadir(p)
    nhalfswath = int((p.half_swath - p.half_gap)/p.delta_ac) + 1
    if p.file_coeff:
        if os.path.isfile(p.file_coeff) and (not p.makesgrid):
            print('\n WARNING: Existing random coefficient file used')
            err.load_coeff(p)
        else:
            err.init_error(p, 2*nhalfswath)
            err.save_coeff(p, 2*nhalfswath)
        if p.nadir:
            if os.path.isfile(p.file_coeff[:-3] + '_nadir.nc')  \
                    and (not p.makesgrid):
                print('WARNING: Existing random nadir coefficient file used')
                errnad.load_coeff(p)
            else:
                errnad.init_error(p)
                errnad.save_coeff(p)
    else:
        err.init_error(p, 2*nhalfswath)
        if p.nadir:
            errnad.init_error(p)
    return err, errnad


def load_sgrid(sgridfile, p):
    '''Load SWOT swath and Nadir data for file sgridfile '''
    import swotsimulator.rw_data as rw_data

    # Load SWOT swath file
    sgrid = rw_data.Sat_SWOT(file=sgridfile)
    cycle = 0
    x_al = []
    x_ac = []
    al_cycle = 0
    timeshift = 0
    sgrid.load_swath(cycle=cycle, x_al=x_al, x_ac=x_ac, al_cycle=al_cycle,
                     timeshift=timeshift)
    sgrid.loncirc = numpy.rad2deg(numpy.unwrap(sgrid.lon))
    # Extract the pass number from the file name
    ipass = int(sgridfile[-6:-3])
    sgrid.ipass = ipass
    return sgrid


def load_ngrid(sgridfile, p):
    ipass = int(sgridfile[-6:-3])
    # Load Nadir track file
    ngrid = rw_data.Sat_nadir(file=p.filesgrid + 'nadir_p'
                              + str(ipass).zfill(3) + '.nc')
    cycle = 0
    x_al = []
    al_cycle = 0
    timeshift = 0
    ngrid.load_orb(cycle=cycle, x_al=x_al, al_cycle=al_cycle,
                   timeshift=timeshift)
    ngrid.loncirc = numpy.rad2deg(numpy.unwrap(ngrid.lon))
    ngrid.ipass = ipass
    return ngrid


'''
def select_modelbox(sgrid, model_data, p):
    # Select data in modelbox from the model 
    model_index = numpy.where(((numpy.min(sgrid.lon)) <= model_data.vlon)
                              & (model_data.vlon <= (numpy.max(sgrid.lon)))
                              & ((numpy.min(sgrid.lat)) <= model_data.vlat)
                              & (model_data.vlat <= (numpy.max(sgrid.lat))))
    model_data.model_index = model_index
    # lonmodel1d = model_data.vlon[model_index].ravel()
    # latmodel1d = model_data.vlat[model_index].ravel()

    if p.grid == 'regular':
        maxlon = numpy.max(sgrid.lon)
        model_data.lon1d = lon1[numpy.where((numpy.min(sgrid.lon) <= lon1)
                                            & (lon1 <= maxlon))]
        model_data.lat1d = lat1[numpy.where((numpy.min(sgrid.lat) <= lat1)
                                            & (lat1 <= numpy.max(sgrid.lat)))]
    else:
        model_index = numpy.where((numpy.min(sgrid.lon) <= model_data.vlon)
                                  & (model_data.vlon <= numpy.max(sgrid.lon))
                                  & (numpy.min(sgrid.lat) <= model_data.vlat)
                                  & (model_data.vlat <= numpy.max(sgrid.lat)))
        model_data.lon1d = model_data.vlon[model_index].ravel()
        model_data.lat1d = model_data.vlat[model_index].ravel()
        model_data.vlon = model_data.vlon[model_index].ravel()
        model_data.vlat = model_data.vlat[model_index].ravel()

    # nx = len(lon1)
    # ny = len(lat1)
    return model_index
'''


def interpolate_regular_1D(lon_in, lat_in, var, lon_out, lat_out, Teval=None):
    ''' Interpolation of data when grid is regular and coordinate in 1D. '''
    # To correct for IDL issues
    ind_sort = numpy.argsort(lon_in)
    lon_in = lon_in[ind_sort]
    var = var[:, ind_sort]
    if numpy.max(lon_in) > 359 and numpy.min(lon_in) < 1:
        ind_in1 = numpy.where(lon_in <= 180)
        ind_in2 = numpy.where(lon_in > 180)
        # lon_in[lon_in > 180] = lon_in[lon_in > 180] - 360
        # lon_in = np.mod(lon_in - (lref - 180), 360) + (lref - 180)
        # lon_in = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(lon_in)))
        ind_out1 = numpy.where(lon_out <= 180)
        ind_out2 = numpy.where(lon_out > 180)
        # lon_out[lon_out > 180] = lon_out[lon_out > 180] - 360
        # lon_out = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(lon_out)))
        interp = interpolate.RectBivariateSpline
        mask_teval = (numpy.isnan(var) | numpy.ma.getmaskarray(var))
        if Teval is None:
            Teval = numpy.zeros(numpy.shape(lon_out))
            if ind_out1[0].any():
                _tmp = interp(lat_in, lon_in[ind_in1],
                              mask_teval[:, ind_in1[0]], kx=1, ky=1,
                              s=0)
                Teval[ind_out1] = _tmp.ev(lat_out[ind_out1], lon_out[ind_out1])
            if ind_out2[0].any():
                _tmp = interp(lat_in, lon_in[ind_in2],
                              mask_teval[:, ind_in2[0]], kx=1, ky=1,
                              s=0)
                Teval[ind_out2] = _tmp.ev(lat_out[ind_out2], lon_out[ind_out2])
        # Trick to avoid nan in interpolation
        var_mask = + var
        var_mask[numpy.isnan(var_mask) | numpy.ma.getmaskarray(var_mask)] = 0.
        # Interpolate variable
        var_out = numpy.full(numpy.shape(lon_out), numpy.nan)
        if ind_out1[0].any():
            _tmp = interp(lat_in, lon_in[ind_in1], var_mask[:, ind_in1[0]],
                          kx=1, ky=1, s=0)
            var_out[ind_out1] = _tmp.ev(lat_out[ind_out1], lon_out[ind_out1])
        if ind_out2[0].any():
            _tmp = interp(lat_in, lon_in[ind_in2], var_mask[:, ind_in2[0]],
                          kx=1, ky=1, s=0)
            var_out[ind_out2] = _tmp.ev(lat_out[ind_out2], lon_out[ind_out2])

    else:
        mask_teval = (numpy.isnan(var) | numpy.ma.getmaskarray(var))
        # Interpolate mask if it has not been done (Teval is None)
        interp = interpolate.RectBivariateSpline
        if Teval is None:
            _Teval = interp(lat_in, lon_in, mask_teval, kx=1, ky=1, s=0)
            Teval = _Teval.ev(lat_out, lon_out)
        # Trick to avoid nan in interpolation
        var_mask = + var
        var_mask[numpy.isnan(var_mask) | numpy.ma.getmaskarray(var_mask)] = 0.
        # Interpolate variable
        _var_out = interp(lat_in, lon_in, var_mask, kx=1, ky=1, s=0)
        var_out = _var_out.ev(lat_out, lon_out)
    # Mask variable with Teval
    var_out[Teval > 0] = numpy.nan
    #var_out[Teval > 0] = numpy.nan
    #var_out[abs(var_out) > 1000] = numpy.nan
    return var_out, Teval


def interpolate_irregular_pyresample(swath_in, var, grid_out, radius,
                                     interp_type='nearest'):
    ''' Interpolation of data when grid is irregular and pyresample is
    installed.'''
    import pyresample as pr
    interp = pr.kd_tree.resample_nearest
    if interp_type == 'nearest':
        interp = pr.kd_tree.resample_nearest
        radius_n = radius * 10**3
        var_out = interp(swath_in, var, grid_out, radius_of_influence=radius_n,
                         epsilon=100)
    else:
        interp = pr.kd_tree.resample_gauss
        radius_g = radius * 3 * 10**3
        sigma_g = radius * 10**3
        var_out = interp(swath_in, var, grid_out, radius_of_influence=radius_g,
                         sigmas=sigma_g, fill_value=None)
    return var_out


def create_SWOTlikedata(cycle, ntotfile, list_file, sgrid, ngrid, model_data,
                        modeltime, err, errnad, p):
    # import swotsimulator.rw_data as rw_data
    # import swotsimulator.build_error as build_error
    import swotsimulator.mod_tools as mod_tools
    # import swotsimulator.const as const
    '''Create SWOT and nadir errors err and errnad, interpolate model SSH
    model_data on swath and nadir track, compute SWOT-like and nadir-like data
    for cycle, SWOT grid sgrid and ngrid. '''
    # - Progress bar variables are global
    global istep
    global ntot
    #   Initialiaze errors and SSH
    progress = 0
    output_var = {}
    shape_i = numpy.shape(sgrid.lon)[0]
    for key in p.list_input_var.keys():
        output_var[key] = numpy.full(shape_i, numpy.nan)
    output_var['vindice'] = numpy.full(shape_i, numpy.nan)

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

    if p.nadir:
        output_var_nadir = {}
        for key in p.list_input_var.keys():
            output_var_nadir[key] = numpy.full((shape_i[0]), numpy.nan)
        output_var_nadir['vindice_nadir'] = numpy.full((shape_i), numpy.nan)
        err.wet_tropo1nadir = numpy.zeros((numpy.shape(ngrid.lon)[0]))
        err.wet_tropo2nadir = numpy.zeros((numpy.shape(ngrid.lon)[0]))
        err.wtnadir = numpy.zeros((numpy.shape(ngrid.lon)[0]))
        errnad.nadir = numpy.zeros((numpy.shape(ngrid.lon)[0]))

    date1 = cycle * sgrid.cycle

    # Definition of the time in the model
    time = sgrid.time  + date1
    # Look for satellite data that are beween step-p.timesetp/2 and
    # step+p.step/2
    if p.file_input is not None:
        lon2D = {}
        lat2D = {}
        # meshgrid in 2D for interpolation purposes
        for key in model_data.vlon.keys():
            if p.grid == 'irregular':
                lon2D[key], lat2D[key] = numpy.meshgrid(model_data.vlon[key],
                                                        model_data.vlat[key])
        index_filemodel = numpy.where(((time[-1] - sgrid.timeshift)
                                      >= (modeltime + p.timestep))
                                      & ((time[0] - sgrid.timeshift)
                                      <= (modeltime - p.timestep)))
        # local variable to find time record in file
        nfile = 0
        time_offset = 0
        # At each step, look for the corresponding time in the satellite data
        for ifile in index_filemodel[0]:
            fstep = float(istep) / float(ntot*ntotfile)
            spass = 'pass: {}'.format(sgrid.ipass)
            cy1 = cycle + 1
            sinfo = 'model file: {}, cycle: {}'.format(list_file[ifile], cy1)
            progress = mod_tools.update_progress(fstep, spass, sinfo)
            ifilenext = ifile + 1
            if len(modeltime) < ifilenext:
                ifilenext = ifile
            # If there are satellite data, Get true SSH from model
            if numpy.shape(index_filemodel)[1] > 0:
                # number of file to be processed used in the progress bar
                ntot = ntot + numpy.shape(index_filemodel)[1] - 1
                # if numpy.shape(index)[1]>1:
                # Select part of the track that corresponds to the time of the
                # model (+-timestep/2)
                ind_time = numpy.where(((time - sgrid.timeshift)
                                       >= modeltime[ifile])
                                       & ((time - sgrid.timeshift)
                                       < modeltime[ifile + 1]))
                #coefficient for temporal interpolation
                alpha = ((time[ind_time] - modeltime[ifile])
                         / (modeltime[ifile + 1] - modeltime[ifile]))
                if p.nadir:
                    ind_nadir_time = numpy.where(((time - ngrid.timeshift)
                                                 >= (modeltime[ifile] - p.timestep/2.))
                                                 & ((time - ngrid.timeshift)
                                                 < (modeltime[ifile] + p.timestep/2.)))
                    alpha_nadir = ((time[ind_nadir_time] - modeltime[ifile])
                                   / (modeltime[ifile + 1] - modeltime[ifile]))
                SSH_model_list = []
                for iifile in (ifile, ifilenext):
                    # pick file time and index
                    filetime = (iifile - time_offset)%p.dim_time
                    # read next file when the time dimension is reached
                    if filetime >= (time_offset + p.dim_time):
                        time_offset += p.dim_time
                        nfile += 1
                        filetime = (iifile - time_offset)%p.dim_time
                    nfile = int(iifile /p.dim_time)
                    filetime = iifile - nfile * p.dim_time

                    _tmpfilename = list_file[nfile]
                    filename = os.path.join(p.indatadir, _tmpfilename)

                    model_step_ctor = getattr(rw_data, model_data.model)
                    model_step = model_step_ctor(p, ifile=(filename, ),
						list_input_var=p.list_input_var,
						time=filetime)
                    input_var = {}

                    if p.grid == 'regular':
                        model_step.read_var(p, ind_lon=model_data.ind_lon)
                        for key in model_step.input_var.keys():
                            grid_key = model_step.numgrid[key]
                            _indlat = model_data.model_index_lat[grid_key]
                            _tmp = model_step.input_var[key][_indlat, :]
                            _indlon = model_data.model_index_lon[grid_key]
                            input_var[key] = +_tmp[:, _indlon]

                    else:
                        model_step.read_var(p, index=None)
                        for key in model_step.input_var.keys():
                            _ind = model_data.model_index[model_step.numgrid[key]]
                           input_var[key] = + model_step.input_var[key][_ind]

                    SSH_model_list.append(input_var)
                # - Interpolate Model data on a SWOT grid and/or along the
                # nadir track
                # if grid is regular, use interpolate.RectBivariateSpline to
                # interpolate
                if p.grid == 'regular' and len(numpy.shape(model_data.vlon)) == 1:
                    # Flatten satellite grid and select part of the track
                    # corresponding to the model time
                    lonswot = sgrid.lon[ind_time[0], :].flatten()
                    lonnadir = ngrid.lon[ind_nadir_time[0]].flatten()
                    latswot = sgrid.lat[ind_time[0], :].flatten()
                    latnadir = ngrid.lat[ind_nadir_time[0]].flatten()
                    interp = interpolate.RectBivariateSpline
                    nal, nac = numpy.shape(sgrid.lon[ind_time[0], :])
                    for key in model_step.input_var.keys():
                        mlon = model_data.vlon[model_step.numgrid[key]]
                        mlat = model_data.vlat[model_step.numgrid[key]]
                        swotlist = []
                        nadirlist = []
                        for input_var_i in SSH_model_list:
                            _tmp, Teval = interpolate_regular_1D(mlon, mlat,
                                                                 input_var_i[key],
                                                                 lonswot,
                                                                 latswot)
                            swotlist.append(_tmp.reshape(nal, nac))
                            if p.nadir:
                                _tmp, Teval = interpolate_regular_1D(mlon, mlat,
                                                                     input_var_i[key],
                                                                     lonnadir,
                                                                     latnadir)
                                nadirlist.append(_tmp)
                        _tmp2 = alpha * swotlist[0] + (1-alpha) * swotlist[1]
                        output_var[key][ind_time[0], :] = + _tmp2
                        if p.nadir:
                            _tmp2 = alpha_nadir * nadirlist[0] + (1-alpha_nadir) * nadirlist[1]
                            output_var_nadir[key][ind_time_nadir[0]] = + _tmp2

                else:
                    # Grid is irregular, interpolation can be done using
                    # pyresample module if it is installed or griddata function
                    # from scipy.
                    # Note that griddata is slower than pyresample functions.
                    try:
                        import pyresample as pr
                        lon_wrap = pr.utils.wrap_longitudes
                        geom = pr.geometry.SwathDefinition
                        sgrid.lon = lon_wrap(sgrid.lon)
                        sigma = max(p.delta_al, p.delta_ac)
                        grid_def = geom(lons=sgrid.lon, lats=sgrid.lat)
                        if p.nadir:
                            ngrid.lon = lon_wrap(ngrid.lon) 
                            ngrid_def = pr.geometry.SwathDefinition(lons=ngrid.lon, lats=ngrid.lat)
                        for key in model_step.input_var.keys():
                            grid_key = model_step.numgrid[key]
                            _ind =  model_data.model_index[grid_key]
                            _mlon = lon_wrap(model_data.vlon[grid_key])
                            _mlat = model_data.vlat[grid_key]
                            if len(_mlon[0]) <= 1:
                                _mlon = lon_wrap(lon2D[grid_key])
                                _mlat = lat2D[grid_key]
                            swath_def = geom(lons=_mlon[_ind], lats=_mlat[_ind])
                            swotlist = []
                            nadirlist = []
                            for input_var_i in SSH_model_list:
                                _tmp = interpolate_irregular_pyresample(
                                          swath_def, input_var_i[key],
                                          grid_def,
                                          sigma,
                                          interp_type=p.interpolation)
                                swotlist.append(_tmp)
                                if p.nadir:
                                    _tmp = interpolate_irregular_pyresample(
                                              swath__def, input_var_i[key],
                                              ngrid_def,
                                              sigma,
                                              interp_type=p.interpolation)
                                    nadirlist.append(_tmp)
                        _tmp2 = alpha * swotlist[0] + (1-alpha) * swotlist[1]
                        output_var[key][ind_time[0], :] = + _tmp2
                        if p.nadir:
                            _tmp2 = alpha_nadir * nadirlist[0] + (1-alpha_nadir) * nadirlist[1]
                            output_var_nadir[key][ind_time[0]] = + _tmp2
  
                    except ImportError:
                        for key in model_step.input_var.keys():
                            grid_key = model_step.numgrid[key]
                            _ind =  model_data.model_index[grid_key]
                            _mlon = model_data.vlon[grid_key]
                            _mlat = model_data.vlat[grid_key]
                            if len(_mlon) <= 1:
                                _mlon = lon2D[grid_key]
                                _mlat = lat2D[grid_key]
                            lonravel = + _mlon[_ind].ravel()
                            latravel = + _mlat[_ind].ravel()
                            swotlist = []
                            nadirlist = []
                            for input_var_i in SSH_model_list:
                                _tmp = + input_var_i[key].ravel()
                                interp = interpolate.griddata((lonravel, latravel),
                                                      _tmp, (sgrid.lon[ind_time[0], :],
                                                      sgrid.lat[ind_time[0], :]),
                                                      method=p.interpolation)
                                swotlist.append(interp)
                                if p.nadir:
                                    interp = interpolate.griddata((lonravel, latravel),
                                               _tmp,
                                               (ngrid.lon[ind_nadir_time[0]],
                                                ngrid.lat[ind_nadir_time[0]]),
                                               method=p.interpolation) 
                                    nadirlist.append(interp)                                   
                            _tmp2 = alpha * swotlist[0] + (1 - alpha) * swotlist[1]
                            output_var[key][ind_time[0], :] = + _tmp2
                            if p.nadir:
                                _tmp2 = alpha_nadir * nadirlist[0] + (1 - alpha_nadir) * nadirlist[1]
                                output_var_nadir[key][ind_time[0]] = + _tmp2
                    
                if p.interpolation == 'nearest':
                    if modelbox[0] > modelbox[1]:
                         for key in model_step.input_var.keys(): 
                             output_var[key][numpy.where(((sgrid.lon < modelbox[0])
                                                         & (sgrid.lon > modelbox[1]))
                                                         | (sgrid.lat < modelbox[2])
                                                         | (sgrid.lat > modelbox[3]))] = numpy.nan
                             if p.nadir:
                                 output_var_nadir[key][numpy.where(((ngrid.lon < modelbox[0])
                                                                    & (ngrid.lon > modelbox[1]))
                                                                    | (ngrid.lat < modelbox[2])
                                                                    | (ngrid.lat > modelbox[3]))] = numpy.nan
                     else:
                         for key in model_step.input_var.keys(): 
                             output_var[key][numpy.where((sgrid.lon < modelbox[0])
                                                         | (sgrid.lon > modelbox[1])
                                                         | (sgrid.lat < modelbox[2])
                                                         | (sgrid.lat > modelbox[3]))] = numpy.nan
                             if p.nadir:
                                  output_var_nadir[key][numpy.where((ngrid.lon < modelbox[0])
                                                                    | (ngrid.lon > modelbox[1])
                                                                    | (ngrid.lat < modelbox[2]) 
                                                                    | (ngrid.lat > modelbox[3]))] = numpy.nan
                output_var['vindice'][[ind_time[0], :] = ifile
                if p.nadir:
                    output_var_nadir['vindice'][ind_nadir_time[0]] = ifile
                del ind_time, model_step, ind_nadir_time
            istep += 1
    else:
        istep += 1
        progress = mod_tools.update_progress(float(istep)/float(ntotfile*ntot),
                                             'pass: ' + str(sgrid.ipass),
                                             'no model file provided' + ', cycle: ' + str(cycle+1))
    err.make_error(sgrid, cycle, SSH_true, p)
    err.make_SSH_error(SSH_true, p)
    if p.nadir:
        errnad.make_error(ngrid, cycle, SSH_true_nadir, p)
        if p.nbeam == 1:
            errnad.SSH = SSH_true_nadir + errnad.nadir + err.wet_tropo1nadir
        else:
            errnad.SSH = SSH_true_nadir + errnad.nadir + err.wet_tropo2nadir
    # if p.file_input: del ind_time, SSH_model, model_step
    return output_var, output_var_nadir, time, progress


def createNadirlikedata(cycle, ntotfile, list_file, ngrid, model_data,
                        modeltime,  errnad, p, progress_bar=False):
    import swotsimulator.rw_data as rw_data
    import swotsimulator.build_error as build_error
    import swotsimulator.mod_tools as mod_tools
    import swotsimulator.const as const

    # - Progress bar variables are global
    global istep
    global ntot
    global ifile
    shape_i = numpy.shape(ngrid.lon)[0]
    errnad.wet_tropo1nadir = numpy.zeros((shape_i))
    errnad.wt = numpy.zeros((shape_i))
    errnad.nadir = numpy.zeros((shape_i))
    output_var_nadir = {}
    for key in p.list_input_var.keys():
        output_var_nadir[key] = numpy.full((shape_i), numpy.nan)
    SSH_true = numpy.zeros((numpy.shape(ngrid.lon)[0]))
    output_var_nadir['vindice'] = numpy.zeros(numpy.shape(ngrid.lon))*numpy.nan
    # Definition of the time in the model
    time = ngrid.time + date1
    # Look for satellite data that are beween step-p.timesetp/2
    # and step+p.step/2
    if p.file_input:
        index_filemodel = numpy.where(((time[-1] - ngrid.timeshift) >= (modeltime-p.timestep))
                                        & ((time[0] - ngrid.timeshift) < (modeltime+p.timestep)))
        # At each step, look for the corresponding time in the satellite data
        for ifile in index_filemodel[0]:
            if progress_bar:
                progress = mod_tools.update_progress(float(istep)/float(ntotfile*ntot), 'pass: ' + str(ngrid.ipass), 'model file: ' + list_file[ifile] + ', cycle:' + str(cycle+1))
            # If there are satellite data, Get true SSH from model
            if numpy.shape(index_filemodel)[1] > 0:
                ntot = ntot + numpy.shape(index_filemodel)[1] - 1
                ind_nadir_time = numpy.where(((time - ngrid.timeshift)
                                                 >= (modeltime[ifile] - p.timestep/2.))
                                                 & ((time - ngrid.timeshift)
                                                 < (modeltime[ifile] + p.timestep/2.)))
                alpha_nadir = ((time[ind_nadir_time] - modeltime[ifile])
                                   / (modeltime[ifile + 1] - modeltime[ifile]))
                SSH_model_list = []
                for iifile in (ifile, ifilenext):
                    # pick file time and index
                    filetime = (iifile - time_offset)%p.dim_time
                    # read next file when the time dimension is reached
                    if filetime >= (time_offset + p.dim_time):
                        time_offset += p.dim_time
                        nfile += 1
                        filetime = (iifile - time_offset)%p.dim_time
                    nfile = int(iifile /p.dim_time)
                    filetime = iifile - nfile * p.dim_time

                    _tmpfilename = list_file[nfile]
                    filename = os.path.join(p.indatadir, _tmpfilename)

                    model_step_ctor = getattr(rw_data, model_data.model)
                    model_step = model_step_ctor(p, ifile=(filename, ),
						list_input_var=p.list_input_var,
						time=filetime)
                    input_var = {}

                    if p.grid == 'regular':
                        model_step.read_var(p, ind_lon=model_data.ind_lon)
                        for key in model_step.input_var.keys():
                            grid_key = model_step.numgrid[key]
                            _indlat = model_data.model_index_lat[grid_key]
                            _tmp = model_step.input_var[key][_indlat, :]
                            _indlon = model_data.model_index_lon[grid_key]
                            input_var[key] = +_tmp[:, _indlon]

                    else:
                        model_step.read_var(p, index=None)
                        for key in model_step.input_var.keys():
                            _ind = model_data.model_index[model_step.numgrid[key]]
                           input_var[key] = + model_step.input_var[key][_ind]

                    SSH_model_list.append(input_var)

            # - Interpolate Model data on a SWOT grid and/or along the nadir
            # track, if grid is regular, use interpolate.RectBivariateSpline to
            # interpolate
            if p.grid == 'regular' and len(numpy.shape(model_data.vlon[0])) == 1:
                # #######################TODO
                # To be moved to routine rw_data
                    # Flatten satellite grid and select part of the track
                    # corresponding to the model time
                    lonnadir = ngrid.lon[ind_nadir_time[0]].flatten()
                    latnadir = ngrid.lat[ind_nadir_time[0]].flatten()
                    interp = interpolate.RectBivariateSpline
                    nal, nac = numpy.shape(sgrid.lon[ind_time[0], :])
                    nadirlist = []
                    for key in model_step.input_var.keys():
                        mlon = model_data.vlon[model_step.numgrid[key]]
                        mlat = model_data.vlat[model_step.numgrid[key]]
                        nadirlist = []
                        for input_var_i in SSH_model_list::
                            _tmp, Teval = interpolate_regular_1D(mlon, mlat,
                                                                 input_var_i[key],
                                                                 lonnadir,
                                                                 latnadir)
                            nadirlist.append(_tmp)

                        _tmp2 = alpha_nadir * nadirlist[0] + (1-alpha_nadir) * nadirlist[1]
                        output_var_nadir[key][ind_time_nadir[0]] = + _tmp2

            else:
                # Grid is irregular, interpolation can be done using pyresample
                # module if it is installed or griddata function from scipy.
                # Note that griddata is slower than pyresample functions.
                try:
                    import pyresample as pr
                    lon_wrap = pr.utils.wrap_longitudes
                    geom = pr.geometry.SwathDefinition
                    sigma = max(p.delta_al, p.delta_ac)
                    ngrid.lon = lon_wrap(ngrid.lon) 
                    ngrid_def = pr.geometry.SwathDefinition(lons=ngrid.lon, lats=ngrid.lat)
                    for key in model_step.input_var.keys():
                        grid_key = model_step.numgrid[key]
                        _ind =  model_data.model_index[grid_key]
                        _mlon = lon_wrap(model_data.vlon[grid_key])
                        _mlat = model_data.vlat[grid_key]
                        if len(_mlon[0]) <= 1:
                            _mlon = lon_wrap(lon2D[grid_key])
                            _mlat = lat2D[grid_key]
                        swath_def = geom(lons=_mlon[_ind], lats=_mlat[_ind])
                        nadirlist = []
                        for input_var_i in SSH_model_list:
                            _tmp = interpolate_irregular_pyresample(
                                             swath__def, input_var_i[key],
                                              ngrid_def,
                                              sigma,
                                              interp_type=p.interpolation)
                            nadirlist.append(_tmp)
                        _tmp2 = alpha_nadir * nadirlist[0] + (1-alpha_nadir) * nadirlist[1]
                        output_var_nadir[key][ind_time[0]] = + _tmp2

                except ImportError:
                    for key in model_step.input_var.keys():
                        grid_key = model_step.numgrid[key]
                        _ind =  model_data.model_index[grid_key]
                        _mlon = model_data.vlon[grid_key]
                        _mlat = model_data.vlat[grid_key]
                        if len(_mlon) <= 1:
                            _mlon = lon2D[grid_key]
                            _mlat = lat2D[grid_key]
                        lonravel = + _mlon[_ind].ravel()
                        latravel = + _mlat[_ind].ravel()
                        nadirlist = []
                        for input_var_i in SSH_model_list:
                            interp = interpolate.griddata((lonravel, latravel),
                                               _tmp,
                                               (ngrid.lon[ind_nadir_time[0]],
                                                ngrid.lat[ind_nadir_time[0]]),
                                               method=p.interpolation) 
                            nadirlist.append(interp)                                   
                        _tmp2 = alpha_nadir * nadirlist[0] + (1 - alpha_nadir) * nadirlist[1]
                        output_var_nadir[key][ind_time[0]] = + _tmp2

                    if p.interpolation == 'nearest':
                        if modelbox[0] > modelbox[1]:
                            for key in model_step.input_var.keys():
                                output_var_nadir[key][numpy.where(((ngrid.lon < modelbox[0])
                                                              & (ngrid.lon > modelbox[1])) 
                                                              | (ngrid.lat < modelbox[2])
                                                              | (ngrid.lat > modelbox[3]))] = numpy.nan
                        else:
                            for key in model_step.input_var.keys():
                                output_var_nadir[key][numpy.where((ngrid.lon < modelbox[0]) 
                                                              | (ngrid.lon > modelbox[1])
                                                              | (ngrid.lat < modelbox[2])
                                                              | (ngrid.lat > modelbox[3]))] = numpy.nan
            vindice_nadir[ind_nadir_time[0]] = ifile
    errnad.make_error(ngrid, cycle, SSH_true_nadir, p)  # , ind[0])
    errnad.SSH = errnad.nadir+err.wet_tropo1nadir
    del SSH_model, model_step, ind_nadir_time
    return output_var_nadir, time, progress



def save_SWOT(cycle, sgrid, err, p, time=[], vindice=[], SSH_true=[]):
    file_output = (p.file_output + '_c' + str(cycle+1).zfill(2) + '_p'
                   + str(sgrid.ipass).zfill(3) + '.nc')
    OutputSWOT = rw_data.Sat_SWOT(file=file_output, lon=(sgrid.lon+360) % 360,
                                  lat=sgrid.lat, time=time, x_ac=sgrid.x_ac,
                                  x_al=sgrid.x_al, cycle=sgrid.cycle,
                                  lon_nadir=(sgrid.lon_nadir+360) % 360,
                                  lat_nadir=sgrid.lat_nadir)
    OutputSWOT.write_data(SSH_model=SSH_true, index=vindice, roll_err=err.roll,
                          bd_err=err.baseline_dilation, phase_err=err.phase,
                          ssb_err=err.ssb, karin_err=err.karin,
                          pd_err_1b=err.wet_tropo1, pd_err_2b=err.wet_tropo2,
                          pd=err.wt, timing_err=err.timing, SSH_obs=err.SSH)
    return None


def save_Nadir(cycle, ngrid, errnad, err , p, time=[], vindice_nadir=[],
               SSH_true_nadir=[]):
    file_outputnadir = (p.file_output + 'nadir_c' + str(cycle+1).zfill(2)
                        + '_p' + str(ngrid.ipass).zfill(3) + '.nc')
    OutputNadir = rw_data.Sat_nadir(file=file_outputnadir,
                                    lon=(ngrid.lon+360) % 360, lat=ngrid.lat,
                                    time=time, x_al=ngrid.x_al,
                                    cycle=ngrid.cycle)
    OutputNadir.write_data(SSH_model=SSH_true_nadir, index=vindice_nadir,
                           nadir_err=errnad.nadir, SSH_obs=errnad.SSH,
                           pd_err_1b=err.wet_tropo1nadir, pd=err.wtnadir,
                           pd_err_2b=err.wet_tropo2nadir)
    return None
