''' Module to read and write data
Contains functions to read variables and coordinates from a netcdf files. \n
Contains model classes: \n
                       -ROMS \n
                       -NEMO \n
                       -NETCDF_MODEL \n
Contains satellite class: Sat_SWOT \n
Contains file instrumentation class: file_instr \n
'''
version = '2.21'
from netCDF4 import Dataset
import numpy
import sys, os
import time as ti
import datetime
import logging
import swotsimulator.const as const
logger = logging.getLogger(__name__)


def read_params(params_file):
    """ Read parameters from parameters file and store it in p.\n
    This program is not used in swot_simulator."""
    import imp
    with open(params_file, 'r') as f:
        p = imp.load_source('p', '', f)
    return p


def write_params(params, pfile):
    """ Write parameters that have been selected to run swot_simulator. """
    with open(pfile, 'w') as f:
        for key in dir(params):
            if not key[0:2] == '__':
                f.write('{} = {}\n'.format(key, params.__dict__[key]))


def read_coordinates(ifile,  nlon, nlat, twoD=True):
    ''' General routine to read coordinates in a netcdf file. \n
    Inputs are file name, longitude name, latitude name. \n
    Outputs are longitudes and latitudes (2D arrays).'''

    # - Open Netcdf file
    try:
        fid = Dataset(ifile, 'r')
    except IOError:
        logger.error('There was an error opening the file {}'.format(ifile))
        sys.exit(1)

    # - Read 1d or 2d coordinates
    try:
        vartmp = fid.variables[nlat]
    except:
        logger.error('Coordinates {} not found in file {}'.format(nlat, ifile))
        sys.exit(1)
    try:
        vartmp = fid.variables[nlon]
    except:
        logger.error('Coordinates {} not found in file {}'.format(nlon, ifile))
        sys.exit(1)
    if len(vartmp.shape) == 1:
        lon_tmp = numpy.array(fid.variables[nlon][:]).squeeze()
        lat_tmp = numpy.array(fid.variables[nlat][:]).squeeze()
        if twoD:
            lon, lat = numpy.meshgrid(lon_tmp, lat_tmp)
        else:
            lon = lon_tmp
            lat = lat_tmp
    elif len(vartmp.shape) == 2:
        lon = numpy.array(fid.variables[nlon][:, :]).squeeze()
        lat = numpy.array(fid.variables[nlat][:, :]).squeeze()
        if not twoD:
            lon = lon[0, :]
            lat = lat[:, 0]
    else:
        logger.warn('unknown dimension for lon and lat')
    fid.close()
    return lon, lat


def read_var(ifile, var, index=None, time=0, depth=0, model_nan=None):
    ''' General routine to read variables in a netcdf file. \n
    Inputs are file name, variable name, index=index to read part
    of the variables, time=time to read a specific time, depth=depth to read a
    specific depth, model_nan=nan value '''

    # - Open Netcdf file
    try:
        fid = Dataset(ifile, 'r')
    except IOError:
        logger.error('There was an error opening the file {}'.format(ifile))
        sys.exit(1)

    # - Check dimension of variable
    try:
        vartmp = fid.variables[var]
    except:
        logger.error('Variable {} not found in file {}'.format(var, ifile))
        sys.exit(1)
    # - Read variable
    if index is None:
        if len(vartmp.shape) == 1:
            T = numpy.array(fid.variables[var][:]).squeeze()
        elif len(vartmp.shape) == 2:
            T = numpy.array(fid.variables[var][:, :]).squeeze()
        elif len(vartmp.shape) == 3:
            if time is None:
                T = numpy.array(fid.variables[var][:, :, :]).squeeze()
            else:
                T = numpy.array(fid.variables[var][time, :, :]).squeeze()
        elif len(vartmp.shape) == 4:
            if time is None:
                if depth is None:
                    T = numpy.array(fid.variables[var][:, :, :, :]).squeeze()
                else:
                    T = numpy.array(fid.variables[var][:, depth, :, :]).squeeze()
            elif depth is None:
                T = numpy.array(fid.variables[var][time, :, :, :]).squeeze()
            else:
                T = numpy.array(fid.variables[var][time, depth, :, :]).squeeze()
        else:
            logger.error('wrong dimension in variables {}'.format(var))
            sys.exit(1)
    else:
        if len(vartmp.shape) == 1:
            Ttmp = numpy.array(fid.variables[var][:]).squeeze()
            T = Ttmp[index]
        elif len(vartmp.shape) == 2:
            Ttmp = numpy.array(fid.variables[var][:, :]).squeeze()
            T = Ttmp[index]
        elif len(vartmp.shape) == 3:
            if time is None:
                U = numpy.array(fid.variables[var][:, :, :]).squeeze()
                T = U[:, index]
            else:
                U = numpy.array(fid.variables[var][time, :, :]).squeeze()
                T = U[index]
        elif len(vartmp.shape) == 4:
            if time is None:
                if depth is None:
                    U = numpy.array(fid.variables[var][:, :, :, :]).squeeze()
                    T = U[:, :, index]
                else:
                    U = numpy.array(fid.variables[var][:, depth, :, :]).squeeze()
                    T = U[:, index]
            elif depth is None:
                U = numpy.array(fid.variables[var][time, :, :, :]).squeeze()
                T = U[:, index]
            else:
                U = numpy.array(fid.variables[var][time, depth, :, :]).squeeze()
                T = U[index]
        else:
            logger.error('wrong dimension in variables {}'.format(var))

    fid.close()
    # - Mask value that are NaN
    if model_nan is not None:
        T[numpy.where(T == model_nan)] = numpy.nan
    T = numpy.ma.masked_invalid(T)
    T[T.mask] = numpy.nan

    return T


def write_var(fid, dim, key, value):
    ''' Write variable key with value value in netcdf file fid'''
    dim_tim, dim_ac, dim_rad = dim
    # Get attribute from dictionary
    vdic = getattr(const, key)
    scale = vdic['scale']
    fill_value = vdic['fill_value']
    # Write variables in 1d
    if len(value.shape) == 1:
        var = fid.createVariable(vdic['varname'], vdic['type'],
                                 (dim_tim,), fill_value=fill_value)
        value = value / scale
        if 'f' not in vdic['type']:
            value = numpy.rint(value)
        value[numpy.isnan(value)] = fill_value
        var[:] = value
    # Write variables in 2d
    if len(value.shape) == 2:
        if value.shape[1] == 2:
            var = fid.createVariable(vdic['varname'], vdic['type'],
                                     (dim_tim, dim_rad),
                                     fill_value=fill_value)
        else:
            var = fid.createVariable(vdic['varname'], vdic['type'],
                                     (dim_tim, dim_ac),
                                     fill_value=fill_value)

        value = value / scale
        if 'f' not in vdic['type']:
            value = numpy.rint(value)
        value[numpy.isnan(value)] = fill_value
        var[:, :] = value
    # Set attributes to variable
    vmin = vdic['min_value']
    if vmin is not None:
        var.valid_min = vmin
    vmax = vdic['max_value']
    if vmax is not None:
        var.valid_max = vmax
    var.units = vdic['unit']
    var.long_name = vdic['longname']
    if scale != 1.:
        var.scale_factor = scale
    if 'comment' in vdic.keys():
        var.comment = vdic['comment']
    if 'standard_name' in vdic.keys():
        var.standard_name = vdic['standard_name']
    if 'calendar' in vdic.keys():
        var.calendar = vdic['calendar']


class Sat_SWOT():
    ''' Sat_SWOT class: to read and write data that has been
    created by SWOT simulator '''
    def __init__(self, nfile=None, lon=None, lat=None, lon_nadir=None,
                 lat_nadir=None, time=None, cycle=None, al_cycle=None,
                 x_al=None, x_ac=None, timeshift=None):
        self.file = nfile
        self.lon = lon
        self.lat = lat
        self.lon_nadir = lon_nadir
        self.lat_nadir = lat_nadir
        self.time = time
        self.cycle = cycle
        self.x_al = x_al
        self.x_ac = x_ac
        self.al_cycle = al_cycle
        self.timeshift = timeshift

    def write_swath(self, **kwargs):
        '''Write swath location in Satellite grid file sgridfile.\n
        Dimensions are  time (i.e. along track), x_ac (across
        track) and cycle (1). \n
        Variables are longitude, latitude, number of days in a cycle,
        distance crossed in a cycle, time, along track and across track
        distances are stored.'''
        # - Open Netcdf file in write mode
        fid = Dataset(self.file, 'w')
        # - Create Global attribute
        fid.title = 'SWOT swath grid simulated by SWOT simulator'
        fid.keywords = 'check keywords'  # Check keywords
        fid.Conventions = "CF-1.6"
        fid.summary = 'SWOT grid data produced'
        fid.description = "SWOT fixed grid"
        fid.Metadata_Conventions = "Unidata Dataset Discovery v1.0"
        fid.history = 'Grid File created by swotsimulator version ' + version
        fid.processing_level = 'L2'
        fid.standard_name_vocabulary = "CF-1.6"
        fid.creator_name = "Lucile Gaultier and Clement Ubelmann"
        fid.creator_email = "lucile.gaultier@gmail.com"
        fid.publisher_url = "github/SWOTSimulator/"
        fid.time_coverage_start = self.time[0]
        # p.date0+"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.time_coverage_end = self.time[-1]
        # p.date0 +"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.geospatial_lat_min = "{:.2f}".format(numpy.min(self.lat))
        fid.geospatial_lat_max = "{:.2f}".format(numpy.max(self.lat))
        fid.geospatial_lat_units = "degrees_north"
        fid.geospatial_lon_max = "{:.2f}".format(numpy.max(self.lon))
        fid.geospatial_lon_min = "{:.2f}".format(numpy.min(self.lon))
        fid.geospatial_lon_units = "degrees_east"
        fid.project = "SWOT"
        fid.date_created = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.date_modified = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.keywords_vocabulary = "NASA"
        fid.references = 'Gaultier, L., C. Ubelmann, and L.-L. Fu, 2016: The '\
                         'Challenge of Using Future SWOT Data for Oceanic '\
                         'Field Reconstruction. J. Atmos. Oceanic Technol.,'\
                         '33, 119-126, doi:10.1175/jtech-d-15-0160.1.'\
                         'http://dx.doi.org/10.1175/JTECH-D-15-0160.1.'
        fid.cycle = "{0:d}".format(int(self.al_cycle))
        # - Create dimensions
        fid.createDimension('time', numpy.shape(self.lon)[0])
        fid.createDimension('x_ac', numpy.shape(self.lon)[1])
        fid.createDimension('cycle', 1)

        # - Create and write Variables
        dim_tim = 'time'
        dim_1 = 'cycle'
        dim_ac = 'x_ac'
        vtime = fid.createVariable('time', 'f8', (dim_tim,))
        vtime_sec = fid.createVariable('time_sec', 'f8', (dim_tim,))
        vlon_nadir = fid.createVariable('lon_nadir', 'f4', (dim_tim,))
        vlat_nadir = fid.createVariable('lat_nadir', 'f4', (dim_tim,))
        vlon = fid.createVariable('lon', 'f4', (dim_tim, dim_ac))
        vlat = fid.createVariable('lat', 'f4', (dim_tim, dim_ac))
        vcycle = fid.createVariable('cycle', 'f4', (dim_1,))
        valcycle = fid.createVariable('al_cycle', 'f4', (dim_1,))
        vtimeshift = fid.createVariable('timeshift', 'f4', (dim_1,))
        vx_al = fid.createVariable('x_al', 'f4', (dim_tim,))
        vx_ac = fid.createVariable('x_ac', 'f4', (dim_ac,))
        vtime[:] = self.time
        vtime.axis = "T"
        vtime.units = "days since the beginning of the sampling"
        vtime.long_name = "Time"
        vtime.standard_name = "time"
        vtime.calendar = "gregorian"
        vtime_sec[:] = self.time * 86400
        vtime.axis = "T"
        vtime.units = "seconds since the beginning of the sampling"
        vtime.long_name = "Time in seconds"
        vtime.standard_name = "time"
        vtime.calendar = "gregorian"
        vlon[:, :] = self.lon
        vlon.axis = "X"
        vlon.long_name = "Longitude"
        vlon.standard_name = "longitude"
        vlon.units = "degrees_east"
        vlat[:, :] = self.lat
        vlat.axis = "Y"
        vlat.long_name = "Latitude"
        vlat.standard_name = "latitude"
        vlat.units = "degrees_north"
        vlon_nadir[:] = self.lon_nadir
        vlon_nadir.axis = "X"
        vlon_nadir.long_name = "Longitude of satellite nadir point"
        vlon_nadir.standard_name = "longitude"
        vlon_nadir.units = "degrees_east"
        vlat_nadir[:] = self.lat_nadir
        vlat_nadir.axis = "Y"
        vlat_nadir.long_name = "Latitude of satellite nadir point"
        vlat_nadir.standard_name = "latitude"
        vlat_nadir.units = "degrees_north"
        vcycle[:] = self.cycle
        vcycle.units = "days during a cycle"
        vcycle.long_name = "Cycle"
        valcycle[:] = self.al_cycle
        valcycle.units = "km"
        valcycle.long_name = " Distance travelled during the pass"
        vtimeshift[:] = self.timeshift
        vtimeshift.units = "day"
        vtimeshift.long_name = "Shift time to match model time"
        vx_al[:] = self.x_al
        vx_al.units = "km"
        vx_al.long_name = "Along track distance from the beginning of the pass"
        vx_ac[:] = self.x_ac
        vx_ac.units = "km"
        vx_ac.long_name = "Across track distance from nadir"
        fid.close()
        return None

    def write_data(self, **kwargs):
        '''Write SWOT data in output file file_output
        Dimensions are x_al (along track distance), x_ac (across
        track distance). \n
        Variables are longitude, latitude, index (file number),
        error-free SSH (SSH interpolated from the model), selected
        errors (karin, wet_tropo, roll, baseline_dilation, phase,
        timing) and SSH with errors. \n
        '''
        # - Open netcdf file in write mode
        fid = Dataset(self.file, 'w') #, format='NC_NETCDF4')
        fid.description = "Ouptut from SWOT simulator"
        if hasattr(self, 'gridfile'):
            fid.corresponding_grid = self.gridfile
        if hasattr(self, 'fid.corresponding_grid'):
            fid.corresponding_grid = self.corresponding_grid

        fid.title = 'SWOT-like data simulated by SWOT simulator'
        fid.keywords = 'SWOT, altimetry, SSH, satellite, remote sensing'
        fid.Conventions = "CF-1.6"
        fid.summary = 'SWOT grid data produced'
        fid.description = "SWOT fixed grid"
        fid.Metadata_Conventions = "Unidata Dataset Discovery v1.0"
        fid.history = 'Grid File created by swotsimulator version ' + version
        fid.processing_level = 'L2'
        fid.standard_name_vocabulary = "CF-1.6"
        fid.creator_name = "Lucile Gaultier and Clement Ubelmann"
        fid.creator_email = "lucile.gaultier@gmail.com"
        fid.publisher_url = "github/SWOTSimulator/"
        fid.time_coverage_start = self.time[0]
        # p.date0+"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.time_coverage_end = self.time[-1]
        # p.date0 +"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.geospatial_lat_min = "{:.2f}".format(numpy.min(self.lat))
        fid.geospatial_lat_max = "{:.2f}".format(numpy.max(self.lat))
        fid.geospatial_lat_units = "degrees_north"
        fid.geospatial_lon_max = "{:.2f}".format(numpy.max(self.lon))
        fid.geospatial_lon_min = "{:.2f}".format(numpy.min(self.lon))
        fid.geospatial_lon_units = "degrees_east"
        fid.project = "SWOT"
        fid.date_created = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.date_modified = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.keywords_vocabulary = "NASA"
        fid.references = 'Gaultier, L., C. Ubelmann, and L.-L. Fu, 2016: The '\
                         'Challenge of Using Future SWOT Data for Oceanic '\
                         'Field Reconstruction. J. Atmos. Oceanic Technol., '\
                         '33, 119-126, doi:10.1175/jtech-d-15-0160.1. '\
                         'http://dx.doi.org/10.1175/JTECH-D-15-0160.1.'
        # fid.cycle = "{0:d}".format(int(self.al_cycle))
        # - Create dimensions
        dim_tim = 'time'
        dim_ac = 'nC'
        dim_rad = 'n2'
        fid.createDimension(dim_tim, numpy.shape(self.lon)[0])
        fid.createDimension(dim_ac, numpy.shape(self.lon)[1])
        fid.createDimension(dim_rad, 2)

        # - Create and write variables
        vlon_nadir = fid.createVariable('lon_nadir', 'i4', (dim_tim,))
        vlat_nadir = fid.createVariable('lat_nadir', 'i4', (dim_tim,))
        vlon = fid.createVariable('lon', 'i4', (dim_tim, dim_ac),
                                  fill_value=2147483647)
        vlat = fid.createVariable('lat', 'i4', (dim_tim, dim_ac),
                                  fill_value=2147483647)
        vx_al = fid.createVariable('x_al', 'f4', (dim_tim,))
        vx_ac = fid.createVariable('x_ac', 'f4', (dim_ac,))
        scale = 1e-6
        vlon[:, :] = numpy.rint(self.lon / scale)
        vlon.units = "deg"
        vlon.long_name = "Longitude"
        vlon.scale_factor = scale
        vlon.valid_min = 0
        vlon.valid_max = 359999999
        vlat[:, :] = numpy.rint(self.lat / scale)
        vlat.units = "deg"
        vlat.long_name = "Latitude"
        vlat.scale_factor = scale
        vlat.valid_min = -80000000
        vlat.valid_max = 80000000
        vlon_nadir[:] = numpy.rint(self.lon_nadir / scale)
        vlon_nadir.units = "deg"
        vlon_nadir.scale_factor = scale
        vlon_nadir.valid_min = 0
        vlon_nadir.valid_max = 359999999
        vlat_nadir[:] = numpy.rint(self.lat_nadir / scale)
        vlat_nadir.units = "deg"
        vlat_nadir.scale_factor = scale
        vlat_nadir.valid_min = -80000000
        vlat_nadir.valid_max = 80000000
        vx_al[:] = self.x_al
        vx_al.units = "km"
        vx_al.long_name = "Along track distance from the beginning of the pass"
        vx_ac[:] = self.x_ac
        vx_ac.units = "km"
        vx_ac.long_name = "Across track distance from nadir"
        dim = [dim_tim, dim_ac, dim_rad]
        if 'empty_var' in kwargs:
            dformat = '%Y-%m-%d %H:%M:%S'
            start_date = datetime.datetime.strptime(self.start_date, dformat)
            first_date = datetime.datetime(2000, 1, 1)
            if start_date < first_date:
                logger.info('start_date has been replaced by 1 January 2000'
                            ' as it is anterior to this date')
            time_day = []
            time_sec = []
            for itime in self.time:
                td = datetime.timedelta(itime) + start_date - first_date
                time_day.append(td.days)
                time_sec.append(td.seconds + td.microseconds * 10**(-6))
            write_var(fid, dim, 'time_sec', numpy.asarray(time_sec))
            write_var(fid, dim, 'time_day', numpy.asarray(time_day))
        else:
            vtime = fid.createVariable('time', 'u8', (dim_tim,))
            vtime[:] = numpy.rint(self.time * 86400 / scale)
            vtime.units = "s"
            vtime.scale_factor = scale
            vtime.valid_min = 0
            vtime.long_name = "Time from beginning of simulation (in s)"


        for key, value in kwargs.items():
            if key == 'empty_var'and value is not None:
                for key2, value2 in value.items():
                    write_var(fid, dim, key2, value2)
            elif key != 'empty_var' and value.any():
                write_var(fid, dim, key, value)
        fid.close()
        return None

    def load_swath(self, **kwargs):
        '''Load swath variables stored in Satellite grid file sgridfile. \n
        (longitude, latitude, number of days in a cycle, crossed distance
        during a cycle, time, along track and across track position).'''

        # - Open Netcdf file
        try:
            fid = Dataset(self.file, 'r')
        except IOError:
            logger.error('There was an error opening the file '
                         '{}'.format(self.file))
            sys.exit(1)
        stime = []
        stime_sec = []
        slon = []
        slat = []
        slon_nadir = []
        slat_nadir = []
        listvar = {'time': stime, 'lon': slon, 'lat': slat,
                   'time_sec': stime_sec,
                   'lon_nadir': slon_nadir, 'lat_nadir': slat_nadir}

        # - Read variables in listvar and return them
        for stringvar in listvar:
            if stringvar not in fid.variables.keys():
                continue
            var = fid.variables[stringvar]
            if len(var.shape) == 1:
                _tmpvar = fid.variables[stringvar][:]
                listvar[stringvar] = numpy.array(_tmpvar).squeeze()
            elif len(var.shape) == 2:
                _tmpvar = fid.variables[stringvar][:, :]
                listvar[stringvar] = numpy.array(_tmpvar).squeeze()
            else:
                logger.error('Wrong number of dimension for variable'
                             ' {}, should be less than 2'.format(var))
                sys.exit(1)
            setattr(self, stringvar, listvar[stringvar])
        # - Read variables in arguments
        for key, value in kwargs.items():
            if key not in fid.variables.keys():
                continue
            var = fid.variables[key]
            if len(var.shape) == 1:
                value = numpy.array(fid.variables[key][:]).squeeze()
            elif len(var.shape) == 2:
                value = numpy.array(fid.variables[key][:, :]).squeeze()
            elif len(var.shape) == 3:
                value = numpy.array(fid.variables[key][:, :, :]).squeeze()
            else:
                logger.error('Wrong number of dimension for variable'
                             ' {}, should be less than 3'.format(var))
                sys.exit(1)
            # value[value == var.fill_value] = numpy.nan
            setattr(self, key, value)
        if hasattr(fid, 'corresponding_grid'):
            self.corresponding_grid = fid.corresponding_grid
        fid.close()
        return None

    def write_expert_data(self, **kwargs):
        '''Write SWOT data in output file file_output
        Dimensions are x_al (along track distance), x_ac (across
        track distance). \n
        Variables are longitude, latitude, index (file number),
        error-free SSH (SSH interpolated from the model), selected
        errors (karin, wet_tropo, roll, baseline_dilation, phase,
        timing) and SSH with errors. \n
        '''
        # - Open netcdf file in write mode
        fid = Dataset(self.file, 'w') #, format='NC_NETCDF4')
        fid.description = "Ouptut from SWOT simulator"
        if hasattr(self, 'gridfile'):
            fid.corresponding_grid = self.gridfile
        if hasattr(self, 'fid.corresponding_grid'):
            fid.corresponding_grid = self.corresponding_grid

        fid.title = 'SWOT-like data simulated by SWOT simulator'
        fid.keywords = 'SWOT, altimetry, SSH, satellite, remote sensing'
        fid.Conventions = "CF-1.6"
        fid.summary = 'SWOT grid data produced'
        fid.description = "SWOT fixed grid"
        fid.Metadata_Conventions = "Unidata Dataset Discovery v1.0"
        fid.history = 'Grid File created by swotsimulator version ' + version
        fid.processing_level = 'L2'
        fid.standard_name_vocabulary = "CF-1.6"
        fid.creator_name = "Lucile Gaultier and Clement Ubelmann"
        fid.creator_email = "lucile.gaultier@gmail.com"
        fid.publisher_url = "github/SWOTSimulator/"
        fid.time_coverage_start = self.time[0]
        # p.date0+"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.time_coverage_end = self.time[-1]
        # p.date0 +"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.geospatial_lat_min = "{:.2f}".format(numpy.min(self.lat))
        fid.geospatial_lat_max = "{:.2f}".format(numpy.max(self.lat))
        fid.geospatial_lat_units = "degrees_north"
        fid.geospatial_lon_max = "{:.2f}".format(numpy.max(self.lon))
        fid.geospatial_lon_min = "{:.2f}".format(numpy.min(self.lon))
        fid.geospatial_lon_units = "degrees_east"
        fid.project = "SWOT"
        fid.date_created = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.date_modified = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.keywords_vocabulary = "NASA"
        fid.references = 'Gaultier, L., C. Ubelmann, and L.-L. Fu, 2016: The '\
                         'Challenge of Using Future SWOT Data for Oceanic '\
                         'Field Reconstruction. J. Atmos. Oceanic Technol., '\
                         '33, 119-126, doi:10.1175/jtech-d-15-0160.1. '\
                         'http://dx.doi.org/10.1175/JTECH-D-15-0160.1.'
        # fid.cycle = "{0:d}".format(int(self.al_cycle))
        # - Create dimensions
        dim_tim = 'time'
        dim_ac = 'nC'
        dim_rad = 'n2'

        fid.createDimension(dim_tim, numpy.shape(self.lon)[0])
        fid.createDimension(dim_ac, numpy.shape(self.lon)[1])
        fid.createDimension(dim_rad, 2)

        # - Create and write variables
        vtime = fid.createVariable('time_sec', 'u8', (dim_tim,))
        vlon_nadir = fid.createVariable('lon_nadir', 'i4', (dim_tim,))
        vlat_nadir = fid.createVariable('lat_nadir', 'i4', (dim_tim,))
        vlon = fid.createVariable('lon', 'i4', (dim_tim, dim_ac),
                                  fill_value=2147483647)
        vlat = fid.createVariable('lat', 'i4', (dim_tim, dim_ac),
                                  fill_value=2147483647)
        vx_al = fid.createVariable('x_al', 'f4', (dim_tim,))
        vx_ac = fid.createVariable('x_ac', 'f4', (dim_ac,))
        scale = 1e-6
        vtime[:] = numpy.rint(self.time * 86400 / scale)
        vtime.units = "s"
        vtime.scale_factor = scale
        vtime.valid_min = 0
        vtime.long_name = "Time from beginning of simulation (in s)"
        vlon[:, :] = numpy.rint(self.lon / scale)
        vlon.units = "deg"
        vlon.long_name = "Longitude"
        vlon.scale_factor = scale
        vlon.valid_min = 0
        vlon.valid_max = 359999999
        vlat[:, :] = numpy.rint(self.lat / scale)
        vlat.units = "deg"
        vlat.long_name = "Latitude"
        vlat.scale_factor = scale
        vlat.valid_min = -80000000
        vlat.valid_max = 80000000
        vlon_nadir[:] = numpy.rint(self.lon_nadir / scale)
        vlon_nadir.units = "deg"
        vlon_nadir.scale_factor = scale
        vlon_nadir.valid_min = 0
        vlon_nadir.valid_max = 359999999
        vlat_nadir[:] = numpy.rint(self.lat_nadir / scale)
        vlat_nadir.units = "deg"
        vlat_nadir.scale_factor = scale
        vlat_nadir.valid_min = -80000000
        vlat_nadir.valid_max = 80000000
        vx_al[:] = self.x_al
        vx_al.units = "km"
        vx_al.long_name = "Along track distance from the beginning of the pass"
        vx_ac[:] = self.x_ac
        vx_ac.units = "km"
        vx_ac.long_name = "Across track distance from nadir"
        dim = [dim_tim, dim_ac, dim_rad]
        for key, value in kwargs.items():
            if key == 'empty_var'and value is not None:
                for key2, value2 in value.items():
                    write_var(fid, dim, key2, value2)
            elif key != 'empty_var' and value.any():
                write_var(fid, dim, key, value)
        fid.close()
        return None


class Sat_nadir():
    def __init__(self, nfile=None, lon=None, lat=None, time=None, cycle=None,
                 al_cycle=None, x_al=None, timeshift=None):
        self.file = nfile
        self.lon = lon
        self.lat = lat
        self.time = time
        self.cycle = cycle
        self.x_al = x_al
        self.al_cycle = al_cycle
        self.timeshift = timeshift

    def load_orb(self, **kwargs):
        '''Load swath variables stored in Satellite grid file sgridfile. \n
        (longitude, latitude, number of days in a cycle, crossed distance
        during a cycle, time, along track and across track position).'''

        # - Open Netcdf file
        try:
            fid = Dataset(self.file, 'r')
        except IOError:
            logger.error('There was an error opening the file'
                         '{}'.format(self.file))
            sys.exit(1)
        stime = []
        slon = []
        slat = []
        listvar = {'time': stime, 'lon': slon, 'lat': slat}

        # - Read variables in listvar and return them
        for stringvar in listvar:
            var = fid.variables[stringvar]
            if len(var.shape) == 1:
                _tmpvar = fid.variables[stringvar][:]
                listvar[stringvar] = numpy.array(_tmpvar).squeeze()
            elif len(var.shape) == 2:
                _tmpvar = fid.variables[stringvar][:, :]
                listvar[stringvar] = numpy.array(_tmpvar).squeeze()
            else:
                logger.error('Wrong number of dimension for variable'
                             ' {}, should be less than 2'.format(var))
                sys.exit(1)
            setattr(self, stringvar, listvar[stringvar])
        # - Read variables in arguments
        for key, value in kwargs.items():
            var = fid.variables[key]
            if len(var.shape) == 1:
                value = numpy.array(fid.variables[key][:]).squeeze()
            elif len(var.shape) == 2:
                value = numpy.array(fid.variables[key][:, :]).squeeze()
            elif len(var.shape) == 3:
                value = numpy.array(fid.variables[key][:, :, :]).squeeze()
            else:
                logger.error('Wrong number of dimension for variable'
                             ' {}, should be less than 3'.format(var))
                sys.exit(1)
            setattr(self, key, value)
        fid.close()
        return None

    def write_data(self, **kwargs):
        '''Write SWOT data in output file file_output
        Dimensions are x_al (along track distance), x_ac (across
        track distance). \n
        Variables are longitude, latitude, index (file number),
        error-free SSH (SSH interpolated from the model), selected
        errors (karin, wet_tropo, roll, baseline_dilation, phase,
        timing) and SSH with errors. \n
        '''
        # - Open netcdf file in write mode
        fid = Dataset(self.file, 'w')
        fid.description = "Orbit computed by SWOT simulator"
        if hasattr(self, 'gridfile'):
            fid.corresponding_grid = self.gridfile
        if hasattr(self, 'corresponding_grid'):
            fid.corresponding_grid = self.corresponding_grid

        fid.title = 'Altimeter like data simulated by SWOT simulator'
        fid.keywords = 'check keywords'  # Check keywords
        fid.Conventions = "CF-1.6"
        fid.summary = 'SWOT grid data produced'
        fid.description = "SWOT fixed grid"
        fid.Metadata_Conventions = "Unidata Dataset Discovery v1.0"
        fid.history = 'Grid File created by swotsimulator version ' + version
        fid.processing_level = 'L2'
        fid.standard_name_vocabulary = "CF-1.6"
        fid.creator_name = "Lucile Gaultier and Clement Ubelmann"
        fid.creator_email = "lucile.gaultier@gmail.com"
        fid.publisher_url = "github/SWOTSimulator/"
        fid.time_coverage_start = self.time[0]
        # p.date0+"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.time_coverage_end = self.time[-1]
        # p.date0 +"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.geospatial_lat_min = "{:.2f}".format(numpy.min(self.lat))
        fid.geospatial_lat_max = "{:.2f}".format(numpy.max(self.lat))
        fid.geospatial_lat_units = "degrees_north"
        fid.geospatial_lon_max = "{:.2f}".format(numpy.max(self.lon))
        fid.geospatial_lon_min = "{:.2f}".format(numpy.min(self.lon))
        fid.geospatial_lon_units = "degrees_east"
        fid.project = "SWOT"
        fid.date_created = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.date_modified = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.keywords_vocabulary = "NASA"
        fid.references = 'Gaultier, L., C. Ubelmann, and L.-L. Fu, 2016: The '\
                         'Challenge of Using Future SWOT Data for Oceanic ' \
                         'Field Reconstruction. J. Atmos. Oceanic Technol., '\
                         '33, 119-126, doi:10.1175/jtech-d-15-0160.1. '\
                         'http://dx.doi.org/10.1175/JTECH-D-15-0160.1.'
        # fid.cycle = "{0:d}".format(int(self.al_cycle))
        # - Create dimensions
        dim_tim = 'time'
        dim_1 = 'cycle'
        fid.createDimension(dim_tim, numpy.shape(self.lon)[0])
        fid.createDimension(dim_1, 1)

        # - Create and write variables
        vtime = fid.createVariable('time', 'u8', (dim_tim,))
        vlon = fid.createVariable('lon', 'i4', (dim_tim,))
        vlat = fid.createVariable('lat', 'i4', (dim_tim,))
        vxal = fid.createVariable('x_al', 'f4', (dim_tim,))
        vcycle = fid.createVariable('ncycle', 'f4', (dim_1,))
        scale = 1e-6
        vtime[:] = numpy.rint(self.time * 86400 / scale)
        vtime.units = "s"
        vtime.long_name = "Time from beginning of simulation (in s)"
        vtime.scale_factor = scale
        vtime.valid_min = 0
        vlon[:] = numpy.rint(self.lon / scale)
        vlon.units = "deg"
        vlon.long_name = "Longitude"
        vlon.scale_factor = scale
        vlon.valid_min = 0
        vlon.valid_max = 359999999
        vlat[:] = numpy.rint(self.lat / scale)
        vlat.units = "deg"
        vlat.long_name = "Latitude"
        vlat.scale_factor = scale
        vlat.valid_min = -80000000
        vlat.valid_max = 80000000
        vxal[:] = self.x_al
        vxal.units = "km"
        vxal.long_name = "Along track distance"
        vcycle[0] = self.cycle
        vcycle.units = "days"
        vcycle.long_name = "Number of days in a cycle"
        dim = [dim_tim, 0, 0]
        for key, value in kwargs.items():
            if value.any():
                write_var(fid, dim, key, value)
        fid.close()
        return None

    def write_orb(self, **kwargs):
        '''Write swath location in Satellite grid file sgridfile.\n
        Dimensions are  x_al (along track distance), x_ac (across
        track distance) and cycle (1). \n
        Variables are longitude, latitude, number of days in a cycle,
        distance crossed in a cycle, time, along track and across track
        distances are stored.'''
        # - Open Netcdf file in write mode
        fid = Dataset(self.file, 'w')
        fid.title = 'Satellite orbit grid computed by SWOT simulator'
        fid.keywords = 'check keywords'  # Check keywords
        fid.Conventions = "CF-1.6"
        fid.summary = 'SWOT grid data produced'
        fid.description = "SWOT fixed grid"
        fid.Metadata_Conventions = "Unidata Dataset Discovery v1.0"
        hist = 'Grid File created by swotsimulator version {}'.format(version)
        fid.history = hist
        fid.processing_level = 'L2'
        fid.standard_name_vocabulary = "CF-1.6"
        fid.creator_name = "Lucile Gaultier and Clement Ubelmann"
        fid.creator_email = "lucile.gaultier@gmail.com"
        fid.publisher_url = "github/SWOTSimulator/"
        fid.time_coverage_start = self.time[0]
        # p.date0+"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.time_coverage_end = self.time[-1]
        # p.date0 +"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.geospatial_lat_min = "{:.2f}".format(numpy.min(self.lat))
        fid.geospatial_lat_max = "{:.2f}".format(numpy.max(self.lat))
        fid.geospatial_lat_units = "degrees_north"
        fid.geospatial_lon_max = "{0:.2f}".format(numpy.max(self.lon))
        fid.geospatial_lon_min = "{0:.2f}".format(numpy.min(self.lon))
        fid.geospatial_lon_units = "degrees_east"
        fid.project = "SWOT"
        fid.date_created = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.date_modified = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.keywords_vocabulary = "NASA"
        fid.references = 'Gaultier, L., C. Ubelmann, and L.-L. Fu, 2016: The '\
                         'Challenge of Using Future SWOT Data for Oceanic '\
                         'Field Reconstruction. J. Atmos. Oceanic Technol., '\
                         '33, 119-126, doi:10.1175/jtech-d-15-0160.1. '\
                         'http://dx.doi.org/10.1175/JTECH-D-15-0160.1.'
        fid.cycle = "{0:d}".format(int(self.al_cycle))
        # - Create dimensions
        # if (not os.path.isfile(self.file)):
        fid.createDimension('time', numpy.shape(self.lon)[0])
        fid.createDimension('cycle', 1)

        # - Create and write Variables
        vtime = fid.createVariable('time', 'f8', ('time',))
        vlon = fid.createVariable('lon', 'f4', ('time',))
        vlat = fid.createVariable('lat', 'f4', ('time',))
        vcycle = fid.createVariable('cycle', 'f4', ('cycle',))
        valcycle = fid.createVariable('al_cycle', 'f4', ('cycle',))
        vtimeshift = fid.createVariable('timeshift', 'f4', ('cycle',))
        vx_al = fid.createVariable('x_al', 'f4', ('time',))
        vtime[:] = self.time
        vtime.units = "days"
        vlon[:] = self.lon
        vlon.units = "deg"
        vlat[:] = self.lat
        vlat.units = "deg"
        vcycle[:] = self.cycle
        vcycle.units = "days"
        vcycle.long_name = "Number of days during a cycle"
        valcycle[:] = self.al_cycle
        valcycle.units = "km"
        valcycle.long_name = " Distance travelled during the pass"
        vtimeshift[:] = self.timeshift
        vtimeshift.units = "km"
        vtimeshift.long_name = "Shift time to match model time"
        vx_al[:] = self.x_al
        vx_al.units = "km"
        vx_al.long_name = "Along track distance from the beginning of the pass"
        fid.close()
        return None


class file_instr():
    '''Class to open file containing instrumentation errors. \n
    USAGE: file_instr(file=file name) \n
    Mandatory argument is the file name. \n'''
    def __init__(self, file=None,):
        self.file = file

    def read_var(self, **kwargs):
        '''Read variables in instrumentation file. \n
        Possible arguments are the variable names in
        the instrument file.'''
        try:
            fid = Dataset(self.file, 'r')
        except IOError:
            logger.error('There was an error opening the file'
                         '{}'.format(self.file))
            sys.exit(1)
        for key, value in kwargs.items():
            var = fid.variables[key]
            if len(var.shape) == 1:
                value = numpy.array(fid.variables[key][:]).squeeze()
            else:
                logger.error('Wrong dimension in instrumentation file')
                sys.exit(1)
            setattr(self, key, value)
        fid.close()
        return None


class file_karin():
    '''Class to open file containing SWH and corresponding karin noise. \n
    USAGE: file_karin(file=file name) \n
    Mandatory argument is the file name. \n'''
    def __init__(self,
                 file=None,):
        self.file = file

    def read_karin(self, swh):
        '''Read and interpolate karin noise. \n
        Possible arguments is SWH. \n
        '''
        try:
            fid = Dataset(self.file, 'r')
        except IOError:
            logger.error('There was an error opening the file '
                         '{}'.format(self.file))
            sys.exit(1)
        var = fid.variables['height_sdt']
        if len(var.shape) == 2:
            hsdt = numpy.array(fid.variables['height_sdt'][:, :]).squeeze()
        else:
            logger.error('Wrong dimension in height_sdt in Karin noise file')
            sys.exit(1)
        var = fid.variables['cross_track']
        if len(var.shape) == 1:
            self.x_ac = numpy.array(fid.variables['cross_track'][:]).squeeze()
        else:
            logger.error('Wrong dimension in x_ac variable in Karin noise'
                         'file')
            sys.exit(1)
        var = fid.variables['SWH']
        if len(var.shape) == 1:
            swh_tmp = numpy.array(fid.variables['SWH'][:]).squeeze()
        else:
            logger.info('Wrong dimension in SWH variable in Karin noise file')
            sys.exit(1)
        i = numpy.argmin(abs(swh_tmp - swh))
        if swh_tmp[i] > swh:
            i += 1
        if numpy.max(swh_tmp) <= swh:
            self.hsdt = hsdt[-1, :]
            logger.warn('WARNING: swh={} is greater than the maximum value in'
                        ' {}, therefore swh is set to the file maximum value:'
                        ' swh='.format(swh, self.file, numpy.max(swh_tmp)))
        else:
            rswh = swh - swh_tmp[i]
            self.hsdt = hsdt[i, :] * (1 - rswh) + rswh * hsdt[i + 1, :]
        fid.close()
        return None


class NEMO():
    '''Class to read NEMO data \n
    USAGE is NEMO(nfile=name of file ,var= variable name,
    lon=longitude name, lat=latitude name, depth= depth name,
    time=time name).\n
    Argument file is mandatory, other arguments have default
    values var='sossheig', lon='nav_lon', lat='nav_lat', depth='depth',
    time='time. \n'''
    def __init__(self, p, nfile=None, var='sossheig', lon='nav_lon',
                 lat='nav_lat', time='time', depth='depth'):
        self.nvar = var
        self.nlon = lon
        self.nlat = lat
        self.ntime = time
        self.nfile = nfile
        self.ndepth = depth
        self.model_nan = getattr(p, 'model_nan', 0)
        p.model_nan = self.model_nan
        self.SSH_factor = getattr(p, 'SSH_factor', 1)
        p.SSH_factor = self.SSH_factor
        self.grid = p.grid

    def read_var(self, index=None):
        '''Read variables from NEMO file \n
        Argument is index=index to load part of the variable.'''
        self.vvar = read_var(self.nfile, self.nvar, index=index, time=0,
                             depth=0, model_nan=self.model_nan)
        self.vvar = self.vvar * self.SSH_factor
        return None

    def read_coordinates(self, index=None):
        '''Read coordinates from NEMO file \n
        Argument is index=index to load part of the variable.'''
        if self.grid == 'regular':
            lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat,
                                        twoD=False)
        else:
            lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat)

        self.vlat = lat
        self.vlon = (lon + 360) % 360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from NEMO file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if (numpy.min(self.vlon) < 1.) and (numpy.max(self.vlon) > 359.):
            self.vlon[numpy.where(self.vlon > 180.)] -= 360
            lon1 = (numpy.min(self.vlon) + 360) % 360
            lon2 = (numpy.max(self.vlon) + 360) % 360
            if lon1 == lon2:
                lon1 = 0
                lon2 = 360
        else:
            lon1 = numpy.min(self.vlon)
            lon2 = numpy.max(self.vlon)
        return [lon1, lon2, numpy.min(self.vlat), numpy.max(self.vlat)]


class ROMS():
    '''Class to read ROMS data \n
    USAGE is ROMS(nfile=name of file ,var= variable name,
    lon=longitude name, lat=latitude name, depth= depth name,
    time=time name).\n
    Argument file is mandatory, other arguments have default
    values var='rho', lon='x_rho', lat='y_rho', depth='depth',
    time='time. \n
    Variable units is specified in params file and default value
    is True (coordinates in degree). \n
    If units is False (coordinates in km), specify left
    low corner of the domain (lon0, lat0) in params file.'''
    def __init__(self, p, nfile=None, var='rho', depth='depth', time='time',
                 lon='x_rho', lat='y_rho'):
        self.nvar = var
        self.nlon = lon
        self.nlat = lat
        self.ntime = time
        self.nfile = nfile
        self.ndepth = depth
        self.model_nan = getattr(p, 'model_nan', 0)
        p.model_nan = self.model
        self.SSH_factor = getattr(p, 'SSH_factor', 1)
        p.SSH_factor = self.SSH_factor
        self.grid = p.grid

    def read_var(self, index=None):
        '''Read variables from ROMS file\n
        Argument is index=index to load part of the variable.'''
        self.vvar = read_var(self.nfile, self.nvar, index=index, time=0,
                             depth=0, model_nan=self.model_nan)
        self.vvar = self.vvar * self.SSH_factor
        return None

    def read_coordinates(self, index=None):
        '''Read coordinates from ROMS file \n
        Argument is index=index to load part of the variable.'''
        if self.grid == 'regular':
            lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat,
                                        twoD=False)
        else:
            lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat)
        self.vlat = lat
        self.vlon = (lon + 360) % 360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from ROMS file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if (numpy.min(self.vlon) < 1.) and (numpy.max(self.vlon) > 359.):
            self.vlon[numpy.where(self.vlon > 180.)] -= 360
            lon1 = (numpy.min(self.vlon) + 360) % 360
            lon2 = (numpy.max(self.vlon) + 360) % 360
        else:
            lon1 = numpy.min(self.vlon)
            lon2 = numpy.max(self.vlon)
        return [lon1, lon2, numpy.min(self.vlat), numpy.max(self.vlat)]


class NETCDF_MODEL():
    '''Class to read any netcdf data.\n
    USAGE is NETCDF_MODEL(nfile=name of file ,var= variable name,
    lon=variable longitude, lat=variable latitude, units=).\n
    Argument file is mandatory, arguments var, lon, lat
    are specified in params file. \n
    '''
    def __init__(self, p, nfile=None, var=None, lon=None, lat=None, depth=0,
                 time=0):
        if var is None:
            self.nvar = p.var
        else:
            self.nvar = var
        if lon is None:
            self.nlon = p.lon
        else:
            self.nlon = lon
        if lat is None:
            self.nlat = p.lat
        else:
            self.nlat = lat
        self.nfile = nfile
        self.depth = depth
        self.time = time
        self.model_nan = getattr(p, 'model_nan', 0)
        p.model_nan = self.model_nan
        self.SSH_factor = getattr(p, 'SSH_factor', 1.)
        p.SSH_factor = self.SSH_factor
        self.grid = p.grid

    def read_var(self, index=None):
        '''Read variables from netcdf file \n
        Argument is index=index to load part of the variable.'''
        self.vvar = read_var(self.nfile, self.nvar, index=index,
                             time=self.time, depth=self.depth,
                             model_nan=self.model_nan)
        self.vvar = self.vvar * self.SSH_factor
        # self.vvar[numpy.where(numpy.isnan(self.vvar))]=0
        return None

    def read_coordinates(self, index=None):
        '''Read coordinates from netcdf file \n
        Argument is index=index to load part of the variable.'''
        if self.grid == 'regular':
            lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat,
                                        twoD=False)
        else:
            lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat)
        self.vlat = lat
        self.vlon = (lon + 360) % 360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from netcdf file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if (numpy.min(self.vlon) < 1.) and (numpy.max(self.vlon) > 359.):
            self.vlon[numpy.where(self.vlon > 180.)] -= 360
            lon1 = (numpy.min(self.vlon) + 360) % 360
            lon2 = (numpy.max(self.vlon) + 360) % 360
        else:
            lon1 = numpy.min(self.vlon)
            lon2 = numpy.max(self.vlon)
        return [lon1, lon2, numpy.min(self.vlat), numpy.max(self.vlat)]


class CLS_MODEL():
    '''Class to read CLS data model type.\n
    USAGE is NETCDF_MODEL(nfile=name of file ,var= variable name,
    lon=variable longitude, lat=variable latitude, units=).\n
    Argument file is mandatory, arguments var, lon, lat
    are specified in params file. \n
    '''
    def __init__(self, p, nfile=None, var=None, lon=None, lat=None, depth=0,
                 time=0):
        if var is None:
            self.nvar = p.var
        if lon is None:
            self.nlon = p.lon
        if lat is None:
            self.nlat = p.lat
        self.nfile = nfile
        self.depth = depth
        self.time = time
        self.p = p
        self.model_nan = getattr(p, 'model_nan', 0.)
        p.model_nan = self.model_nan
        self.SSH_factor = getattr(p, 'SSH_factor', 1)
        p.SSH_factor = self.SSH_factor
        self.grid = p.grid

    def read_var(self, index=None):
        '''Read variables from netcdf file \n
        Argument is index=index to load part of the variable.'''
        self.vvar = numpy.transpose(read_var(self.nfile, self.nvar,
                                    index=index, time=self.time,
                                    depth=self.depth,
                                    model_nan=self.model_nan))
        self.vvar = self.vvar * self.SSH_factor
        # self.vvar[numpy.where(numpy.isnan(self.vvar))]=0
        return None

    def read_coordinates(self, index=None):
        '''Read coordinates from netcdf file \n
        Argument is index=index to load part of the variable.'''
        if self.grid == 'regular':
            lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat,
                                        twoD=False)
        else:
            lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat)
        self.vlat = lat
        self.vlon = (lon + 360) % 360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from netcdf file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if (numpy.min(self.vlon) < 1.) and (numpy.max(self.vlon) > 359.):
            self.vlon[numpy.where(self.vlon > 180.)] -= 360
            lon1 = (numpy.min(self.vlon) + 360) % 360
            lon2 = (numpy.max(self.vlon) + 360) % 360
        else:
            lon1 = numpy.min(self.vlon)
            lon2 = numpy.max(self.vlon)
        return [lon1, lon2, numpy.min(self.vlat), numpy.max(self.vlat)]


class MITgcm():
    '''Class to read MITgcm binary data.\n
    USAGE is MITgcm(nfile=name of file ,var= variable name,
    lon=variable longitude, lat=variable latitude, units=).\n
    Argument file is mandatory, arguments var, lon, lat
    are specified in params file. \n
    '''
    def __init__(self, p, nfile=None, var='Eta', lon='XC', lat='YC', depth=0,
                 time=0):
        self.nvar = var
        self.nlon = lon
        self.nlat = lat
        self.nfile = nfile
        self.depth = depth
        self.time = time
        self.model_nan = p.model_nan = 0
        self.SSH_factor = getattr(p, 'SSH_factor', 1)
        p.SSH_factor = self.SSH_factor
        self.indir = p.indatadir
        self.shape = (p.Nx, p.Ny)

    def read_var(self, index=None):
        '''Read variables from binary file \n
        Argument is index=index to load part of the variable.'''
        self.vvar = (numpy.fromfile(self.nfile, '>f4', shape=self.shape,
                                    mode='r')[:] * self.SSH_factor)
        return None

    def read_coordinates(self, index=None):
        '''Read MITgcm output coordinates saved in XC.data and YC.data \n
        Argument is index=index to load part of the variable.'''
        lon = numpy.memmap(os.path.join(self.indir, 'XC.data'), '>f4',
                           shape=self.shape, mode='r')
        lat = numpy.memmap(os.path.join(self.indir, 'YC.data'), '>f4',
                           shape=self.shape, mode='r')
        self.vlat = lat[:, :]
        self.vlon = (lon[:, :] + 360) % 360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from netcdf file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if (numpy.min(self.vlon) < 1.) and (numpy.max(self.vlon) > 359.):
            self.vlon[numpy.where(self.vlon > 180.)] -= 360
            lon1 = (numpy.min(self.vlon) + 360) % 360
            lon2 = (numpy.max(self.vlon) + 360) % 360
        else:
            lon1 = numpy.min(self.vlon)
            lon2 = numpy.max(self.vlon)
        return [lon1, lon2, numpy.min(self.vlat), numpy.max(self.vlat)]
