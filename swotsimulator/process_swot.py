import swotsimulator.rw_data as rw
import swotsimulator.build_error as build_error
import swotsimulator.settings as settings
import numpy
import scipy
import re
import os
import glob
import netCDF4
import swotsimulator.mod_run as mod
from typing import Optional
import pyresample as pr


class S3NGgrid():
    def __init__(self, mdata, out_var: Optional[dict] = None):
        if out_var is not None:
            self.lon = out_var['lon']
            self.lat = out_var['lat']
            self.lon_nadir = out_var['lon'][:, 0]
            self.lat_nadir = out_var['lat'][:, 0]
            self.lat = out_var['lat']
            self.x_al = out_var['xal']
            self.x_ac = out_var['xac']
            self.al_cycle = out_var['xal'][-1]
        self.gridfile = mdata.nfile
        self.cycle = mdata.cycle
        self.ipass = mdata.ipass


def interpolate_irregular_pyresample(swath_in, var, grid_out, radius,
                                     interp_type='nearest'):
    ''' Interpolation of data when grid is irregular and pyresample is
    installed.'''
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


def get_global_attr(nfile):
    pattern = "SWOT_L3_LR_SSH_Basic_(\d{3})_(\d{3})_(\d{8})T(\d{6})_(\d{8})T(\d{6})_v2.0.0.nc"
    #match = re.compile(pattern).search(nfile)
    match = re.findall(pattern, nfile)[0]
    cycle = int(match[0])
    ipass = int(match[1])
    return cycle, ipass


def lookforfile(mdata):
    inputdir = '/mnt/data_12t/swot_swh'
    pattern = 'SWOT_L2_LR_SSH_WindWave'
    swhpat = os.path.join(inputdir, f'cycle_{mdata.cycle:03d}',
                          f'{pattern}_{mdata.cycle:03d}_{mdata.ipass:03d}_*')
    listfile = glob.glob(swhpat)
    if len(listfile) == 0:
        nfile = None
    else:
        nfile = listfile[0]
    return nfile


def read_ssh_file(nfile):
    mdata = rw.SWOT_L3_REAL(p, nfile=nfile)
    mdata.read_var()
    mdata.read_coordinates()
    if (numpy.max(mdata.vlon) - numpy.min(mdata.vlon)) > 180:
        mdata.vlon = numpy.mod(mdata.vlon + 180, 360) - 180
    mdata.input_var['ssh_true'] = numpy.ma.masked_invalid(mdata.input_var['ssh_true'])
    mdata.x_ac = numpy.arange(-68, 69, 2) * 10**3
    #mdata.x_al = numpy.aranen(mdata.input_var['time'][:]), 1) * 2 * 10**3
    mdata.x_al = numpy.arange(0, len(mdata.input_var['time'][:]), 1) * 2 * 10**3
    if (numpy.max(mdata.vlon) - numpy.min(mdata.vlon)) > 180:
        mdata.lon = numpy.mod(mdata.lon + 180, 360) - 180
    mdata.cycle, mdata.ipass = get_global_attr(nfile)
    mdata.nfile = nfile
    return mdata


def read_swh_file(nfile: str, mdata, var: Optional[str] = 'swh_karin'):
    swhdata = S3NGgrid(mdata)
    fid = netCDF4.Dataset(nfile, 'r')
    swhdata.vlon = fid['longitude'][:]
    if (numpy.max(swhdata.vlon) - numpy.min(swhdata.vlon)) > 180:
        swhdata.vlon = numpy.mod(swhdata.vlon + 180, 360) - 180
    swhdata.vlat = fid['latitude'][:]
    swhdata.input_var = {'swh': fid[var][:]}
    swhmodel = fid['swh_model'][:]
    swhdata.input_var['time'] = fid['time'][:]
    fid.close()
    swhdata.input_var['swh'][abs(swhdata.input_var['swh'].data) > 20] = numpy.nan
    swhdata.input_var['swh'] = numpy.ma.masked_invalid(swhdata.input_var['swh'])
    _mask = swhdata.input_var['swh'].mask
    swhdata.input_var['swh'][_mask] = swhmodel[_mask]
    swhdata.input_var['swh'][abs(swhdata.input_var['swh'].data) > 60] = numpy.nan
    swhdata.input_var['swh'] = swhdata.input_var['swh'].data
    # numpy.ma.masked_invalid(swhdata.input_var['swh'])
    swhdata.x_ac = numpy.arange(-68, 69, 2) * 10**3
    swhdata.x_al = numpy.arange(0, len(swhdata.input_var['time'][:]),
                                1) * 2 * 10**3
    return swhdata


def interpolate(mdata, var: str, len_xal: float):
    s3ng_xac = numpy.arange(-70 * 10**3, 71 * 10**3, 5 * 10**3)
    s3ng_xal = numpy.arange(0, len_xal, 5 * 10**3)

    funclon = scipy.interpolate.RectBivariateSpline(mdata.x_al, mdata.x_ac,
                                                    mdata.vlon)
    funclat = scipy.interpolate.RectBivariateSpline(mdata.x_al, mdata.x_ac,
                                                    mdata.vlat)

    lon_final = funclon(s3ng_xal, s3ng_xac)
    lat_final = funclat(s3ng_xal, s3ng_xac)
    time_final = []
    if 'time' in mdata.input_var.keys():
        functim = scipy.interpolate.interp1d(mdata.x_al,
                                             mdata.input_var['time'].data)
        time_final = functim(s3ng_xal)
    wrap_lon = pr.utils.wrap_longitudes
    geomdef = pr.geometry.SwathDefinition
    interp = interpolate_irregular_pyresample
    lon_model = wrap_lon(mdata.vlon)
    sigm = 1.
    swath_def = geomdef(lons=lon_model, lats=mdata.vlat)
    lon_grid = wrap_lon(lon_final)
    grid_def = geomdef(lons=lon_grid[:, :],
                       lats=lat_final[:, :])
    _ssh = interp(swath_def, mdata.input_var[var], grid_def,
                  max(p.delta_al, p.delta_ac) * sigm,
                  interp_type='gaussian')
    return {'time': time_final,
            'xac': s3ng_xac,
            'xal': s3ng_xal,
            'lon': lon_final,
            'lat': lat_final,
            var: _ssh}


def load_error(p, nadir_alone=False, seed=0):
    '''Initialize random coefficients that are used to compute
    random errors following the specified spectrum. \n
    If a random coefficient file is specified, random coefficients
    are ldir_alone=Falseoaded from this file.
    '''
    err = build_error.error(p)
    err.init_error2(p)  # , 2*nhalfswath, seed=seed)
    return err


def create_error(err, out_var, sgrid, p, time):
    if 'swh' in out_var.keys():
        swh = out_var['swh']
    else:
        swh = None
    #for key in (p.noise):
    #    print('makeerr', key)
    #import pdb ; pdb.set_trace()
    err.make_error(sgrid, sgrid.cycle, out_var['ssh_true'], p, time, swh=swh)
    if p.product_type != 'expert':
        # TODO err.reconstruct_2D(p, sgrid.x_ac)
        err.make_SSH_error(out_var['ssh_true'], p)
    return err


def create_SWOTlikedata(nfile, p, err):
    mdata = read_ssh_file(nfile)
    swhfile = lookforfile(mdata)
    print(swhfile)
    if swhfile is not None:
        swhdata = read_swh_file(swhfile, mdata)
    else:
        print('missing SWH file')
        return None
    out_var = interpolate(mdata, 'ssh_true', mdata.x_al[-1])
    if swhfile is not None:
        swh_var = interpolate(swhdata, 'swh', mdata.x_al[-1])
    out_var['swh'] = swh_var['swh']
    sgrid = S3NGgrid(mdata, out_var=out_var)
    sgrid.x_ac = sgrid.x_ac / 1000
    sgrid.x_al = sgrid.x_al / 1000
    out_var['time'] = (out_var['time'] - (numpy.datetime64(p.first_time)
                            - numpy.datetime64("2000-01-01 00:00:00")).astype("float")) / 86400
    #out_var['time'] = out_var['time'] / 86400.

    err = create_error(err, out_var, sgrid, p, out_var['time'])
    #err.make_SSH_error(out_var['ssh_true'], p))
    print(sgrid.cycle, sgrid.ipass)
    mod.save_SWOT(sgrid.cycle - 1, sgrid, err, p, out_var, time=out_var['time'],
                  save_var=p.product_type)


class par():
    def __init__(self):
        self.var = 'ssha_unfiltered' # 'ssha_noiseless'
        self.list_input_var = {'ssh_true': [self.var, '']} #, 'ssha': ['ssha', '']}
        self.lon = 'longitude'
        self.lat = 'latitude'
        self.model_nan = -2147483647.0
        self.delta_al = 5
        self.delta_ac = 5
        self.seed = 1
        self.file_systematic = '/mnt/data_12t/swot/Input_syst_errors/S3A_1year_baseline_and_phase_error.nc'
        self.karin_noise = '/mnt/data_12t/swot/Input_syst_errors/SAOOH_random_error_optimB1_step2.nc'
        self.product_type = 'Classic'
        self.file_output = 'SWOT_S3NGT'
        # par.noise = ["karins3ng", "wet_troposphere", "systematic_errors3ng"]
        #self.noise = ["Karins3ng", "systematic_errors3ng", "SystematicErrors3ng", "WetTroposphere", "karins3ng", "wet_troposphere", ]
        self.noise = ["SystematicErrors3ng", "WetTroposphere", "Karins3ng",]
        #self.noise = ["WetTroposphere", "Karins3ng",]
        self.file_input = 'S3NG'
        self.first_time = '2023-07-26T00:00:00Z'
        self.nbeam = 3
        self.sigma = 10
        self.beam_position = [-45, 0, 45]
        self.len_repeat = 40000
        self.lambda_cut = 20000
        self.lambda_max = self.lambda_cut
        self.rng = settings.Seed(self.seed)
        self.height = 891 * 10**3
        self.add_systematic_error = False



if '__main__' == __name__:
    p = par()
    print('load_error')
    err = load_error(p)
    cycle = 8
    #for ipass in range(1, 29):
    for cycle in range(cycle, cycle + 1):
    #for cycle in range(1, 19):
        path = os.path.join('/mnt/data_12t/swot/', f'cycle_{cycle:03d}')
        for ipass in range(244, 585):
            print(ipass)
            patt = os.path.join(path,
                                f'SWOT_L3_LR_SSH_Basic_{cycle:03d}_{ipass:03d}*v2.0.0.nc')
            if len(glob.glob(patt)) == 0:
                print(f'File {patt} not found')
                continue
            nfile = glob.glob(patt)[0]
            print(nfile)
            # nfile = '/mnt/data_12t/swot/cycle_013/SWOT_L3_LR_SSH_Basic_013_001_20240327T143444_20240327T152611_v1.0.2.nc'
            create_SWOTlikedata(nfile, p, err)
