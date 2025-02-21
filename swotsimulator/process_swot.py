import swotsimulator.rw_data as rw
import swotsimulator.build_error as build_error
import numpy
import scipy
import pyresample as pr

class S3NGgrid(mdata):
    def __init__():
        

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

    
def read_ssh_file(nfile):
    mdata = rw.SWOT_L3_REAL(p, nfile=nfile)
    mdata.read_var()
    mdata.read_coordinates()
    mdata.input_var['ssh_true'] = numpy.ma.masked_invalid(mdata.input_var['ssh_true'])
    mdata.x_ac = numpy.arange(-68, 69, 2)*10**3
    mdata.x_al = numpy.arange(0, len(mdata.input_var['time'][:]), 1) *2 *10**3
    mdata.x_al = numpy.arange(0, len(mdata.input_var['time'][:]), 1) *2 *10**3
    if (numpy.max(mdata.vlon) - numpy.min(mdata.vlon)) > 180:
        mdata.lon = numpy.mod(mdata.lon + 180, 360) - 180
    return mdata

def interpolate(mdata):
    swh = 2 * numpy.ones(mdata.input_var['ssh_true'].shape)
    s3ng_xac = numpy.arange(-70*10**3, 71*10**3, 5*10**3)
    s3ng_xal = numpy.linspace(0, mdata.x_al[-1], 5*10**3)

    funclon = scipy.interpolate.RectBivariateSpline(mdata.x_al, mdata.x_ac, 
                                             (numpy.mod(mdata.vlon +180, 360) -180))
    funclat = scipy.interpolate.RectBivariateSpline(mdata.x_al, mdata.x_ac, 
                                             mdata.vlat)
    functim = scipy.interpolate.interp1d(mdata.x_al, mdata.input_var['time'].data)
    
    lon_final = funclon(s3ng_xal, s3ng_xac)
    lat_final = funclat(s3ng_xal, s3ng_xac)
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
    _ssh = interp(swath_def, mdata.input_var['ssh_true'], grid_def,
                  max(p.delta_al, p.delta_ac)*sigm,
                  interp_type='gaussian')
    sgrid = 
    return {'time': time_final,
            'xac': s3ng_xac,
            'xal': s3ng_xal,
            'lon': lon_final,
            'lat': lat_final,
            'ssh_true': _ssh}


def load_error(p, nadir_alone=False, seed=0):
    '''Initialize random coefficients that are used to compute
    random errors following the specified spectrum. \n
    If a random coefficient file is specified, random coefficients
    are ldir_alone=Falseoaded from this file.
    '''
    err = build_error.error(p)
    err.init_error2(p) #, 2*nhalfswath, seed=seed)
    return err


def create_error(out_var, sgrid, p, time):
    if 'swh' in out_var.keys():
        swh = out_var['swh']
    else:
        swh = None
    for key in (p.noise):
        print('makeerr', key)
    err.make_error(sgrid, cycle, out_var['ssh_true'], p, time, swh=swh)
    print(err.karin)
    if p.product_type != 'expert':
        # TODO err.reconstruct_2D(p, sgrid.x_ac)
        err.make_SSH_error(out_var['ssh_true'], p)


def create SWOTlikedata(nfile, p):
    mdata = read_ssh_file(nfile)
    out_var = interpolate(mdata)
    sgrid = S3NGgrid()
