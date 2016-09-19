'''
Module to read and write data
Contains functions to read variables and coordinates from a netcdf files. \n
Contains model classes: \n
                       -ROMS \n 
                       -NEMO \n 
                       -NETCDF_MODEL \n
Contains satellite class: Sat_SWOT \n
Contains file instrumentation class: file_instr \n
'''
netcdf4=True
try:
    from netCDF4 import Dataset
except ImportError:
    print('WARNING: package netCDF4 missing, scipy.io.netcdf is used instead of netCDF4')
    from scipy.io.netcdf import netcdf_file as Dataset
    netcdf4=False
#from scipy.io.netcdf import netcdf_file as Dataset
import swotsimulator.const as const
import numpy
import sys, os
import time as ti
try: import params as p
except: 
    print('params.py not found')
    sys.exit()
#from scipy.io import netcdf as nc
#from scipy.io.netcdf import netcdf as nc

def read_params(params_file):
    """ Read parameters from parameters file and store it in p.\n
    This program is not used in swot_simulator."""
    import imp
    f=open(params_file)
    #global p
    p = imp.load_source('p', '',f)
    f.close()
    return p

def write_params(params, pfile):
    """ Write parameters that have been selected to run swot_simulator. """
    f=open(pfile, 'w')
    for key in dir(params):
        if not key[0:2]=='__':
            f.write(key+ ' = '+str(params.__dict__[key])+'\n')
    f.close()

def read_coordinates(file,  nlon, nlat, twoD=True ):
    ''' General routine to read coordinates in a netcdf file. \n
    Inputs are file name, longitude name, latitude name. \n
    Outputs are longitudes and latitudes (2D arrays).'''

## - Open Netcdf file
    try : 
        fid = Dataset(file, 'r')
    except IOError: 
        print('There was an error opening the file '+file)
        sys.exit()
    #fid = Dataset(file, 'r') 

## - Read 1d or 2d coordinates
    try:
        vartmp = fid.variables[nlat]
    except: 
        sys.exit('Coordinates '+ nlat+ ' not found in file ' +file)

    try:
        vartmp = fid.variables[nlon]
    except: 
        sys.exit('Coordinates '+ nlon+ ' not found in file ' +file)
        #sys.exit()
    if  len(vartmp.shape) == 1:
        lon_tmp = numpy.array(fid.variables[nlon][:]).squeeze()
        lat_tmp = numpy.array(fid.variables[nlat][:]).squeeze()
        if twoD:
          lon, lat = numpy.meshgrid(lon_tmp, lat_tmp)
        else:
            lon=lon_tmp ; lat=lat_tmp
    elif len(vartmp.shape) == 2:
        lon = numpy.array(fid.variables[nlon][:,:]).squeeze()
        lat = numpy.array(fid.variables[nlat][:,:]).squeeze()
        if not twoD:
            lon=lon[0,:]
            lat=lat[:,0]
    else:
        print('unknown dimension for lon and lat')
    fid.close()
    return lon, lat




def read_var(file, var, index=None, time=0, depth=0, model_nan=None):
    ''' General routine to read variables in a netcdf file. \n
    Inputs are file name, variable name, index=index to read part 
    of the variables, time=time to read a specific time, depth=depth to read a 
    specific depth, model_nan=nan value '''

## - Open Netcdf file
    try : 
    #fid = nc.netcdf_file(file, 'r')
        fid = Dataset(file, 'r')
    except IOError: 
        print('There was an error opening the file '+file)
        sys.exit()
    #fid = Dataset(file, 'r') 

## - Check dimension of variable
    try : vartmp = fid.variables[var]
    except: 
        sys.exit('Variable '+ var + ' not found in file ' +file)
## - Read variable
    if index is None:
        if len(vartmp.shape) == 1:
            T = numpy.array(fid.variables[var][:]).squeeze()
        elif len(vartmp.shape) == 2 :
            T = numpy.array(fid.variables[var][:,:]).squeeze()
        elif len(vartmp.shape) == 3 :
            if time is None:
                T = numpy.array(fid.variables[var][:,:,:]).squeeze()
            else:
                T = numpy.array(fid.variables[var][time,:,:]).squeeze()
        elif len(vartmp.shape) == 4 :
            if time is None:
                if depth is None:
                    T = numpy.array(fid.variables[var][:,:,:,:]).squeeze()
                else:
                    T = numpy.array(fid.variables[var][:,depth,:,:]).squeeze()
            elif depth is None:
                T = numpy.array(fid.variables[var][time,:,:,:]).squeeze()
            else:
                T = numpy.array(fid.variables[var][time,depth,:,:]).squeeze()
    else :
            if len(vartmp.shape) == 1:
                Ttmp = numpy.array(fid.variables[var][:]).squeeze()
                T = Ttmp[index] #numpy.array(fid.variables[var][index]).squeeze()
            elif len(vartmp.shape) == 2 :
                Ttmp = numpy.array(fid.variables[var][:,:]).squeeze()
                T=Ttmp[index]
            elif len(vartmp.shape) == 3 :
                if time is None:
                    U = numpy.array(fid.variables[var][:,:,:]).squeeze()
                    T=U[:,index]
                else:
                    U = numpy.array(fid.variables[var][time,:,:]).squeeze()
                    T=U[index]
            elif len(vartmp.shape) == 4 :
                if time is None:
                    if depth is None:
                        U = numpy.array(fid.variables[var][:,:,:,:]).squeeze()
                        T=U[:,:,index]
                    else:
                        U = numpy.array(fid.variables[var][:,depth,:,:]).squeeze()
                        T=U[:,index]
                elif depth is None:
                    U = numpy.array(fid.variables[var][time,:,:,:]).squeeze()
                    T=U[:,index]
                else:
                    U = numpy.array(fid.variables[var][time,depth,:,:]).squeeze()
                    T=U[index]

            #print 'Unsupported number of dimensions'
            #sys.exit()
    fid.close()
## - Mask value that are NaN
    if not model_nan is None:
        T[numpy.where(T==model_nan)]=numpy.nan
    return T

class Sat_SWOT():
    ''' Sat_SWOT class: to read and write data that has been 
    created by SWOT simulator '''
    def __init__(self, 
                file=None, 
                lon=None,
                lat=None,
                lon_nadir=None,
                lat_nadir=None,
                time=None,
                cycle=None,
                al_cycle=None,
                x_al=None,
                x_ac=None, 
                timeshift=None):
        self.file=file
        self.lon=lon
        self.lat=lat
        self.lon_nadir=lon_nadir
        self.lat_nadir=lat_nadir
        self.time=time
        self.cycle=cycle
        self.x_al=x_al
        self.x_ac=x_ac
        self.al_cycle=al_cycle
        self.timeshift=timeshift

    def write_swath(self, **kwargs): 
        '''Write swath location in Satellite grid file sgridfile.\n
        Dimensions are  time (i.e. along track), x_ac (across 
        track) and cycle (1). \n
        Variables are longitude, latitude, number of days in a cycle, 
        distance crossed in a cycle, time, along track and across track 
        distances are stored.'''
## - Open Netcdf file in write mode
        if netcdf4:
          fid = Dataset(self.file, 'w', format='NETCDF4_CLASSIC')
        else:
          fid = Dataset(self.file, 'w')
        ## - Create Global attribute
        fid.title = 'SWOT grid'
        fid.keywords = 'check keywords'  # Check keywords
        fid.Conventions = "CF-1.6"
        fid.summary = 'SWOT grid data produced'
        fid.description = "SWOT fixed grid"
        fid.Metadata_Conventions = "Unidata Dataset Discovery v1.0"
        fid.history = 'Grid File created by swotsimulator version '  # Add version
        fid.processing_level = 'L2'
        fid.standard_name_vocabulary = "CF-1.6"
        fid.creator_name = "Lucile Gaultier and Clement Ubelmann"
        fid.creator_email = "lucile.gaultier@gmail.com"
        fid.publisher_url = "github/SWOTSimulator/"
        fid.time_coverage_start = "YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.time_coverage_end = "YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.geospatial_lat_min = numpy.min(lat).format(%2f)
        fid.geospatial_lat_max = numpy.max(lat)
        fid.geospatial_lat_units = "degrees_north" ;
        fid.geospatial_lon_max = numpy.max(lon) ;
        fid.geospatial_lon_min = numpy.min(lon) ;
        fid.geospatial_lon_units = "degrees_east";
        fid.project = "SWOT"
        fid.date_created = "today"  # to be modified
        fid.date_modified = "??"  # tbm
        fid.keywords_vocabulary = "NASA"
        fid.references = "Gaultier, L., C. Ubelmann, and L.-L. Fu, 2016: The Challenge of Using Future SWOT Data for Oceanic Field Reconstruction. J. Atmos. Oceanic Technol., 33, 119–126, doi:10.1175/jtech-d-15-0160.1. http://dx.doi.org/10.1175/JTECH-D-15-0160.1."
        fid.cycle = str(int(self.al_cycle))
        ## - Create dimensions
        #if (not os.path.isfile(self.file)):
        fid.createDimension('time', numpy.shape(self.lon)[0])
        #fid.createDimension('time_nadir', numpy.shape(self.lon)[0])
        fid.createDimension('x_ac', numpy.shape(self.lon)[1])
        #fid.createDimension('cycle', 1)

## - Create and write Variables
        vtime = fid.createVariable('time', 'f', ('time',))
        vlon_nadir = fid.createVariable('lon_nadir', 'f4', ('time',))
        vlat_nadir = fid.createVariable('lat_nadir', 'f4', ('time',))
        #vlon_nadir = fid.createVariable('lon_nadir', 'f4', ('time_nadir',))
        #vlat_nadir = fid.createVariable('lat_nadir', 'f4', ('time_nadir',))
        vlon = fid.createVariable('lon', 'f4', ('time','x_ac'))
        vlat = fid.createVariable('lat', 'f4', ('time','x_ac'))
        #vcycle = fid.createVariable('cycle', 'f4', ('cycle',))
        valcycle = fid.createVariable('al_cycle', 'f4', ('cycle',))
        vtimeshift = fid.createVariable('timeshift', 'f4', ('cycle',))
        vx_al = fid.createVariable('x_al', 'f4', ('time',))
        vx_ac = fid.createVariable('x_ac', 'f4', ('x_ac',))
        vtime[:]=self.time
        vtime.axis="T"
        vtime.units="days since the beginning of the sampling"
        vtime.long_name="Time"
        vtime.standard_name="time"
        vtime.calendar="gregorian"
        vlon[:,:]=self.lon
        vlon.axis="X"
        vlon.long_name="Longitude"
        vlon.standard_name="longitude"
        vlon.units="degrees_east"
        vlat[:,:]=self.lat
        vlat.axis="Y"
        vlat.long_name="Latitude"
        vlat.standard_name="latitude"
        vlat.units="degrees_north"
        vlon_nadir[:]=self.lon_nadir
        vlon.axis="X"
        vlon.long_name="Longitude"
        vlon.standard_name="longitude"
        vlon.units="degrees_east"
        vlat_nadir[:]=self.lat_nadir
        vlat.axis="Y"
        vlat.long_name="Latitude"
        vlat.standard_name="latitude"
        vlat.units="degrees_north"
        vcycle[:]=self.cycle
        vcycle.units="days during a cycle"
        vcycle.long_name="Cycle"
        valcycle[:]=self.al_cycle
        valcycle.units="km"
        valcycle.long_name=" Distance travelled during the pass"
        vtimeshift[:]=self.timeshift
        vtimeshift.units="km"
        vtimeshift.long_name="Shift time to match model time"
        vx_al[:]=self.x_al
        vx_al.units="km"
        vx_al.long_name="Along track distance from the beginning of the pass"
        vx_ac[:]=self.x_ac
        vx_ac.units="km"
        vx_ac.long_name="Across track distance from nadir"
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
## - Open netcdf file in write mode
        if netcdf4: 
          fid = Dataset(self.file, 'w', format='NETCDF4_CLASSIC') 
        else:  
          fid = Dataset(self.file, 'w')
        fid.description = "Ouptut from SWOT simulator"
        try: fid.corresponding_grid=self.gridfile
        except: pass
## - Create dimensions
        fid.createDimension('time', numpy.shape(self.lon)[0])
        #fid.createDimension('time_nadir', numpy.shape(self.lon_nadir)[0])
        fid.createDimension('x_ac', numpy.shape(self.lon)[1])

## - Create and write variables
        vtime = fid.createVariable('time', 'f', ('time',))
        #vlon_nadir = fid.createVariable('lon_nadir', 'f4', ('time_nadir',))
        #vlat_nadir = fid.createVariable('lat_nadir', 'f4', ('time_nadir',))
        vlon_nadir = fid.createVariable('lon_nadir', 'f4', ('time',))
        vlat_nadir = fid.createVariable('lat_nadir', 'f4', ('time',))
        vlon = fid.createVariable('lon', 'f4', ('time','x_ac'))
        vlat = fid.createVariable('lat', 'f4', ('time','x_ac'))
        vx_al = fid.createVariable('x_al', 'f4', ('time',))
        vx_ac = fid.createVariable('x_ac', 'f4', ('x_ac',))
        vtime[:]=self.time
        vtime.units="days"
        vtime.long_name="Time from beginning of simulation"
        vlon[:,:]=self.lon
        vlon.units="deg"
        vlon.long_name="Longitude"
        vlat[:,:]=self.lat
        vlat.units="deg"
        vlat.long_name="Latitude"
        vlon_nadir[:]=self.lon_nadir
        vlon_nadir.units="deg"
        vlat_nadir[:]=self.lat_nadir
        vlat_nadir.units="deg"
        vx_al[:]=self.x_al
        vx_al.units="km"
        vx_al.long_name="Along track distance from the beginning of the pass"
        vx_ac[:]=self.x_ac
        vx_ac.units="km"
        vx_ac.long_name="Across track distance from nadir"
        longname={"pd": "Simulated path delay error due to wet tropo", "SSH_model":"SSH interpolated from model", "SSH_obs":"Observed SSH (SSH_model+errors)", "index":"Equivalent model output number in list of file", "pd_err_1b":"Residual path delay error after a 1-beam radiometer correction", "pd_err_2b":"Residual path delay error after a 2-beams radiometer correction", "roll_err":"Residual roll error", "phase_err": "Phase error", "timing_err": "Timing error", "bd_err": "Baseline dilation error", "karin_err": "Karin instrument random error", "nadir_err":"Nadir error"}
        unit={"pd":"m", "SSH_model":"m", "SSH_obs":"m","index":" ", "pd_err_1b":"m", "pd_err_2b": "m", "roll_err":"m", "phase_err":"m", "timing_err": "m", "bd_err":"m", "karin_err":"m", "nadir_err":"m"}
        for key, value in kwargs.items():
            #if not value is None:
            if value.any():
                if len(value.shape) == 1 :
                    var = fid.createVariable(str(key), 'f4', ('time',), fill_value=-1.36e9)
                    vmin=numpy.nanmin(value)
                    vmax=numpy.nanmax(value)
                    value[numpy.isnan(value)]=-1.36e9
                    value[value==0]=-1.36e9
                    var[:] = value
                if len(value.shape) == 2 :
                    var = fid.createVariable(str(key), 'f4', ('time', 'x_ac'), fill_value=-1.36e9)
                    vmin=numpy.nanmin(value)
                    vmax=numpy.nanmax(value)
                    value[numpy.isnan(value)]=-1.36e9
                    value[value==0]=-1.36e9
                    var[:,:] = value
                try: var.units = unit[str(key)]
                except: var.units=''
                try:    var.long_name = longname[str(key)]
                except: var.long_name=str(key)
                #try:    var.missing_value = p.model_nan
                #except: var.missing_value = 0. 
                #fid.setncattr('missing_value','-9999.f')

        fid.close()
        return None

    def load_swath(self,**kwargs):
        '''Load swath variables stored in Satellite grid file sgridfile. \n
        (longitude, latitude, number of days in a cycle, crossed distance 
        during a cycle, time, along track and across track position).'''

## - Open Netcdf file
        try : 
            fid = Dataset(self.file, 'r')
        except IOError: 
            print('There was an error opening the file '+self.file)
            sys.exit()
        #fid = Dataset(self.file, 'r')
        stime=[]; slon=[]; slat=[]; slon_nadir=[]; slat_nadir=[]; cycle=[]; x_al=[]; x_ac=[]
        listvar={'time':stime, 'lon': slon, 'lat' : slat, 'lon_nadir': slon_nadir, 'lat_nadir' : slat_nadir}

## - Read variables in listvar and return them
        for stringvar in listvar:
            var = fid.variables[stringvar]
            if len(var.shape) == 1 :
                listvar[stringvar] = numpy.array(fid.variables[stringvar][:]).squeeze()
            elif len(var.shape) == 2 :
                listvar[stringvar] = numpy.array(fid.variables[stringvar][:,:]).squeeze()
            #listvar[stringvar][listvar[stringvar]==var.fill_value] = numpy.nan
            exec('self.'+stringvar+' = listvar[stringvar]') #listvar[stringvar]=var

## - Read variables in arguments
        for key, value in kwargs.items():
            var = fid.variables[key]
            if len(var.shape) == 1 :
                value = numpy.array(fid.variables[key][:]).squeeze()
            elif len(var.shape) == 2 :
                value = numpy.array(fid.variables[key][:,:]).squeeze()
            elif len(var.shape) == 3 :
                value = numpy.array(fid.variables[key][:,:,:]).squeeze()
            #value[value == var.fill_value] = numpy.nan
            exec('self.'+key+' = value')
        try: self.corresponding_grid=fid.corresponding_grid
        except: pass
        fid.close()
        return None

class Sat_nadir():
    def __init__(self,
                file=None,
                lon=None,
                lat=None,
                time=None,
                cycle=None,
                al_cycle=None,
                x_al=None,
                timeshift=None):
        self.file=file
        self.lon=lon
        self.lat=lat
        self.time=time
        self.cycle=cycle
        self.x_al=x_al
        self.al_cycle=al_cycle
        self.timeshift=timeshift

    def load_orb(self,**kwargs):
        '''Load swath variables stored in Satellite grid file sgridfile. \n
        (longitude, latitude, number of days in a cycle, crossed distance
        during a cycle, time, along track and across track position).'''

## - Open Netcdf file
        try :
            fid = Dataset(self.file, 'r')
        except IOError:
            print('There was an error opening the file '+self.file)
            sys.exit()
        #fid = Dataset(self.file, 'r')
        stime=[]; slon=[]; slat=[]; tcycle=[]; x_al=[]; x_ac=[]
        listvar={'time':stime, 'lon': slon, 'lat' : slat}

## - Read variables in listvar and return them
        for stringvar in listvar:
            var = fid.variables[stringvar]
            if len(var.shape) == 1 :
                listvar[stringvar] = numpy.array(fid.variables[stringvar][:]).squeeze()
            elif len(var.shape) == 2 :
                listvar[stringvar] = numpy.array(fid.variables[stringvar][:,:]).squeeze()
            exec('self.'+stringvar+' = listvar[stringvar]') #listvar[stringvar]=var

## - Read variables in arguments
        for key, value in kwargs.items():
            var = fid.variables[key]
            if len(var.shape) == 1 :
                value = numpy.array(fid.variables[key][:]).squeeze()
            elif len(var.shape) == 2 :
                value = numpy.array(fid.variables[key][:,:]).squeeze()
            elif len(var.shape) == 3 :
                value = numpy.array(fid.variables[key][:,:,:]).squeeze()
            exec('self.'+key+' = value')
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
## - Open netcdf file in write mode
        if netcdf4: 
          fid = Dataset(self.file, 'w', format='NETCDF4_CLASSIC') 
        else:  
          fid = Dataset(self.file, 'w')
        fid.description = "Orbit computed by SWOT simulator"
        try: fid.corresponding_grid=self.gridfile
        except: pass
## - Create dimensions
        fid.createDimension('time', numpy.shape(self.lon)[0])
        fid.createDimension('cycle', 1)

## - Create and write variables
        vtime = fid.createVariable('time', 'f', ('time',))
        vlon = fid.createVariable('lon', 'f4', ('time',))
        vlat = fid.createVariable('lat', 'f4', ('time',))
        vxal = fid.createVariable('x_al', 'f4', ('time',))
        vcycle = fid.createVariable('ncycle', 'f4', ('cycle',))
        vtime[:]=self.time
        vtime.units="days"
        vtime.long_name="Time from beginning of simulation"
        vlon[:]=self.lon
        vlon.units="deg"
        vlon.long_name="Longitude"
        vlat[:]=self.lat
        vlat.units="deg"
        vlat.long_name="Latitude"
        vxal[:]=self.x_al
        vxal.units="km"
        vxal.long_name="Along track distance"
        vcycle[0]=self.cycle
        vcycle.units="days"
        vcycle.long_name="Number of days in a cycle"
        longname={"SSH_model":"SSH interpolated from model", "SSH_obs":"Observed SSH (SSH_model+errors)", "index":"Equivalent model output number in list of file", "nadir_err":"Nadir error", "pd_err_1b":"Residual path delay error after a 1-beam radiometer correction", "pd": "Simulated path delay error due to wet tropo"}
        unit={"SSH_model":"m", "SSH_obs":"m","index":" ", "nadir_err":"m", "pd_err_1b":"m","pd":"m","pd_err_2b":"m"}
        for key, value in kwargs.items():
            #if not value is None:
            if value.any():
                if len(value.shape) == 1 :
                    var = fid.createVariable(str(key), 'f4', ('time',))
                    var[:] = value
                    #if len(value.shape) == 2 :
                    #var = fid.createVariable(str(key), 'f4', ('x_al'))
                    #var[:,:] = value
 #                  elif len(value.shape) == 3 :
#                   fid.createDimension('time', numpy.shape(value)[2])
  #                 var = fid.createVariable(str(key), 'f4', ('x_al', 'x_ac', 'time'))
   #                var[:,:,:] = value
                try:    var.units = unit[str(key)]
                except: var.units=''
                try:    var.long_name = longname[str(key)]
                except: var.long_name=str(key)
                try:    var.missing_value = p.model_nan
                except: var.missing_value = 0.
                #fid.setncattr('missing_value','-9999.f')
        try: fid.corresponding_grid=fid.corresponding_grid
        except: pass
        fid.close()
        return None

    def write_orb(self, **kwargs):
        '''Write swath location in Satellite grid file sgridfile.\n
        Dimensions are  x_al (along track distance), x_ac (across
        track distance) and cycle (1). \n
        Variables are longitude, latitude, number of days in a cycle,
        distance crossed in a cycle, time, along track and across track
        distances are stored.'''
## - Open Netcdf file in write mode
        if netcdf4: 
          fid = Dataset(self.file, 'w', format='NETCDF4_CLASSIC') 
        else:  
          fid = Dataset(self.file, 'w' )
        fid.description = "Orbit computed from SWOT simulator"

## - Create dimensions
        #if (not os.path.isfile(self.file)):
        fid.createDimension('time', numpy.shape(self.lon)[0])
        fid.createDimension('cycle', 1)

## - Create and write Variables
        vtime = fid.createVariable('time', 'f', ('time',))
        vlon = fid.createVariable('lon', 'f4', ('time',))
        vlat = fid.createVariable('lat', 'f4', ('time',))
        vcycle = fid.createVariable('cycle', 'f4', ('cycle',))
        valcycle = fid.createVariable('al_cycle', 'f4', ('cycle',))
        vtimeshift = fid.createVariable('timeshift', 'f4', ('cycle',))
        vx_al = fid.createVariable('x_al', 'f4', ('time',))
        vtime[:]=self.time
        vtime.units="days"
        vlon[:]=self.lon
        vlon.units="deg"
        vlat[:]=self.lat
        vlat.units="deg"
        vcycle[:]=self.cycle
        vcycle.units="days"
        vcycle.long_name="Number of days during a cycle"
        valcycle[:]=self.al_cycle
        valcycle.units="km"
        valcycle.long_name=" Distance travelled during the pass"
        vtimeshift[:]=self.timeshift
        vtimeshift.units="km"
        vtimeshift.long_name="Shift time to match model time"
        vx_al[:]=self.x_al
        vx_al.units="km"
        vx_al.long_name="Along track distance from the beginning of the pass"
        fid.close()
        return None

class file_instr():
    '''Class to open file containing instrumentation errors. \n
    USAGE: file_instr(file=file name) \n
    Mandatory argument is the file name. \n'''
    def __init__(self,
                 file=None,):
        self.file=file

    def read_var(self, **kwargs):
        '''Read variables in instrumentation file. \n 
        Possible arguments are the variable names in 
        the instrument file.'''
        try : 
            fid = Dataset(self.file, 'r')
        except : 
            print('There was an error opening the file '+self.file)
            sys.exit()
        for key, value in kwargs.items():
            var = fid.variables[key]
            if len(var.shape) == 1 :
                value = numpy.array(fid.variables[key][:]).squeeze()
            else: 
                print('Wrong dimension in instrumentation file')
                sys.exit()
            exec('self.'+key+' = value')
        fid.close()
        return None

class file_karin():
    '''Class to open file containing SWH and corresponding karin noise. \n
    USAGE: file_karin(file=file name) \n
    Mandatory argument is the file name. \n'''
    def __init__(self,
                 file=None,):
        self.file=file

    def read_karin(self, swh):
        '''Read and interpolate karin noise. \n
        Possible arguments is SWH. \n
        '''
        try : 
            fid = Dataset(self.file, 'r')
        except IOError: 
            print('There was an error opening the file '+self.file)
            sys.exit()
        var=fid.variables['height_sdt']
        if len(var.shape) == 2:
            hsdt=numpy.array(fid.variables['height_sdt'][:,:]).squeeze()
        else: 
            print('Wrong dimension in height_sdt in Karin noise file')
            sys.exit()
        var=fid.variables['cross_track']
        if len(var.shape) == 1:
            self.x_ac=numpy.array(fid.variables['cross_track'][:]).squeeze()
        else: 
            print( 'Wrong dimension in x_ac variable in Karin noise file')
            sys.exit()
        var=fid.variables['SWH']
        if len(var.shape) == 1:
            swh_tmp=numpy.array(fid.variables['SWH'][:]).squeeze()
        else: 
            print('Wrong dimension in SWH variable in Karin noise file')
            sys.exit()
        i=numpy.argmin(abs(swh_tmp-swh))
        if swh_tmp[i]>swh:
            i+=1
        if numpy.max(swh_tmp)<=swh:
            self.hsdt=hsdt[-1,:]
            print('WARNING: swh='+str(swh)+' is greater than the maximum value in '+self.file+', therefore swh is set to the file maximum value: swh='+str(numpy.max(swh_tmp)))
        else:
            rswh=swh-swh_tmp[i]
            self.hsdt=hsdt[i,:]*(1-rswh)+rswh*hsdt[i+1,:]
        fid.close()
        return None

  
class NEMO():
    '''Class to read NEMO data \n
    USAGE is NEMO(file=name of file ,var= variable name,
    lon=longitude name, lat=latitude name, depth= depth name,
    time=time name).\n
    Argument file is mandatory, other arguments have default
    values var='sossheig', lon='nav_lon', lat='nav_lat', depth='depth',
    time='time. \n'''
    def __init__(self,
                file=None,
                var='sossheig',
                lon='nav_lon',
                lat='nav_lat',
                time='time',
                depth='depth',
                ):
        self.nvar = var
        self.nlon = lon
        self.nlat = lat
        self.ntime=time
        self.nfile = file
        self.ndepth=depth
        try: self.model_nan=p.model_nan
        except: self.model_nan=0. ; p.model_nan=0.
    
    def read_var(self, index=None):
        '''Read variables from NEMO file \n 
        Argument is index=index to load part of the variable.'''
        try :SSH_factor=p.SSH_factor
        except: SSH_factor=1. ; p.SSH_factor=SSH_factor
        self.vvar = read_var(self.nfile, self.nvar, index=index,time=0, depth=0, model_nan=self.model_nan)*SSH_factor
        return None

    def read_coordinates(self, index=None):
        '''Read coordinates from NEMO file \n
        Argument is index=index to load part of the variable.'''
        if p.grid=='regular':
          lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat, twoD=False)
        else:
          lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat)
        
        self.vlat=lat
        self.vlon=(lon+360)%360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from NEMO file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if numpy.min(self.vlon)<1. and numpy.max(self.vlon)>359.:
            self.vlon[numpy.where(self.vlon>180.)]=self.vlon[numpy.where(self.vlon>180.)]-360
            lon1=(numpy.min(self.vlon)+360)%360
            lon2=(numpy.max(self.vlon)+360)%360
            if lon1==lon2: 
              lon1=0 
              lon2=360
        else:
            lon1=numpy.min(self.vlon)
            lon2=numpy.max(self.vlon)
        return [lon1, lon2, numpy.min(self.vlat), numpy.max(self.vlat)]

class ROMS():
    '''Class to read ROMS data \n
    USAGE is ROMS(file=name of file ,var= variable name,
    lon=longitude name, lat=latitude name, depth= depth name, 
    time=time name).\n
    Argument file is mandatory, other arguments have default
    values var='rho', lon='x_rho', lat='y_rho', depth='depth',
    time='time. \n
    Variable units is specified in params file and default value
    is True (coordinates in degree). \n
    If units is False (coordinates in km), specify left
    low corner of the domain (lon0, lat0) in params file.'''
    def __init__(self,
                file=None,
                var='rho',
                depth='depth',
                time='time', 
                lon='x_rho', 
                lat='y_rho',
                ):
        self.nvar = var
        self.nlon = lon
        self.nlat = lat
        self.ntime=time
        self.nfile = file
        self.ndepth = depth
        try : self.model_nan=p.model_nan
        except: self.model_nan=0.; p.model_nan=0.

    def read_var(self, index=None):
        '''Read variables from ROMS file\n
        Argument is index=index to load part of the variable.'''
        try: SSH_factor=p.SSH_factor
        except: SSH_factor=1. ; p.SSH_factor=1.
        self.vvar = read_var(self.nfile, self.nvar, index=index, time=0, depth=0, model_nan=self.model_nan)*SSH_factor
        return None


    def read_coordinates(self, index=None):
        '''Read coordinates from ROMS file \n
        Argument is index=index to load part of the variable.'''
        if p.grid=='regular':
          lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat,twoD=False)
        else:
          lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat)
        self.vlat=lat
        self.vlon=(lon+360)%360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from ROMS file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if numpy.min(self.vlon)<1. and numpy.max(self.vlon)>359.:
            self.vlon[numpy.where(self.vlon>180.)]=self.vlon[numpy.where(self.vlon>180.)]-360
            lon1=(numpy.min(self.vlon)+360)%360
            lon2=(numpy.max(self.vlon)+360)%360
        else:
            lon1=numpy.min(self.vlon)
            lon2=numpy.max(self.vlon)
        return [lon1, lon2, numpy.min(self.vlat), numpy.max(self.vlat)]



class NETCDF_MODEL():
    '''Class to read any netcdf data.\n 
    USAGE is NETCDF_MODEL(file=name of file ,var= variable name,
    lon=variable longitude, lat=variable latitude, units=).\n
    Argument file is mandatory, arguments var, lon, lat
    are specified in params file. \n
    '''
    def __init__(self,
                file=None,
                var=p.var,
                lon=p.lon,
                lat=p.lat,
                depth=0,
                time=0,
                ):
        self.nvar = var
        self.nlon = lon
        self.nlat = lat
        self.nfile = file
        self.depth=depth
        self.time=time
        try: self.model_nan=p.model_nan
        except: self.model_nan=0. ; p.model_nan=0.

    def read_var(self, index=None):
        '''Read variables from netcdf file \n
        Argument is index=index to load part of the variable.'''
        try: SSH_factor=p.SSH_factor
        except: SSH_factor=1. ; p.SSH_factor=1.
        self.vvar = read_var(self.nfile, self.nvar, index=index, time=self.time, depth=self.depth, model_nan=self.model_nan)*SSH_factor
        #self.vvar[numpy.where(numpy.isnan(self.vvar))]=0
        return None

    def read_coordinates(self, index=None):
        '''Read coordinates from netcdf file \n
        Argument is index=index to load part of the variable.'''
        if p.grid=='regular':
            lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat, twoD=False)
        else:
            lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat)
        self.vlat=lat
        self.vlon=(lon+360)%360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from netcdf file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if numpy.min(self.vlon)<1. and numpy.max(self.vlon)>359.:
            self.vlon[numpy.where(self.vlon>180.)]=self.vlon[numpy.where(self.vlon>180.)]-360
            lon1=(numpy.min(self.vlon)+360)%360
            lon2=(numpy.max(self.vlon)+360)%360
        else:
            lon1=numpy.min(self.vlon)
            lon2=numpy.max(self.vlon)
        return [lon1, lon2, numpy.min(self.vlat), numpy.max(self.vlat)]

class CLS_MODEL():
    '''Class to read any netcdf data.\n
    USAGE is NETCDF_MODEL(file=name of file ,var= variable name,
    lon=variable longitude, lat=variable latitude, units=).\n
    Argument file is mandatory, arguments var, lon, lat
    are specified in params file. \n
    '''
    def __init__(self,
                file=None,
                var=p.var,
                lon=p.lon,
                lat=p.lat,
                depth=0,
                time=0,
                ):
        self.nvar = var
        self.nlon = lon
        self.nlat = lat
        self.nfile = file
        self.depth=depth
        self.time=time
        try: self.model_nan=p.model_nan
        except: self.model_nan=0. ; p.model_nan=0.

    def read_var(self, index=None):
        '''Read variables from netcdf file \n
        Argument is index=index to load part of the variable.'''
        try: SSH_factor=p.SSH_factor
        except: SSH_factor=1. ; p.SSH_factor=1.
        self.vvar = numpy.transpose(read_var(self.nfile, self.nvar, index=index, time=self.time, depth=self.depth, model_nan=self.model_nan)*SSH_factor)
        #self.vvar[numpy.where(numpy.isnan(self.vvar))]=0
        return None

    def read_coordinates(self, index=None):
        '''Read coordinates from netcdf file \n
        Argument is index=index to load part of the variable.'''
        if p.grid=='regular':
            lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat, twoD=False)
        else:
            lon, lat = read_coordinates(self.nfile, self.nlon, self.nlat)
        self.vlat=lat
        self.vlon=(lon+360)%360
            #lon=numpy.transpose(lon)
            #lat=numpy.transpose(lat)
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from netcdf file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if numpy.min(self.vlon)<1. and numpy.max(self.vlon)>359.:
            self.vlon[numpy.where(self.vlon>180.)]=self.vlon[numpy.where(self.vlon>180.)]-360
            lon1=(numpy.min(self.vlon)+360)%360
            lon2=(numpy.max(self.vlon)+360)%360
        else:
            lon1=numpy.min(self.vlon)
            lon2=numpy.max(self.vlon)
        return [lon1, lon2, numpy.min(self.vlat), numpy.max(self.vlat)]



class MITgcm():
    '''Class to read MITgcm binary data.\n
    USAGE is MITgcm(file=name of file ,var= variable name,
    lon=variable longitude, lat=variable latitude, units=).\n
    Argument file is mandatory, arguments var, lon, lat
    are specified in params file. \n
    '''
    def __init__(self,
                file=None,
                var='Eta',
                lon='XC',
                lat='YC',
                depth=0,
                time=0,
                ):
        self.nvar = var
        self.nlon = lon
        self.nlat = lat
        self.nfile = file
        self.depth=depth
        self.time=time
        self.model_nan = p.model_nan = 0

    def read_var(self, index=None):
        '''Read variables from binary file \n
        Argument is index=index to load part of the variable.'''
        try: SSH_factor=p.SSH_factor
        except: SSH_factor = p.SSH_factor = 1.
        self.vvar = numpy.fromfile(self.nfile,'>f4',shape=(p.Ny,p.Nx),mode='r')[:]*SSH_factor
        return None

    def read_coordinates(self, index=None):
        '''Read MITgcm output coordinates saved in XC.data and YC.data \n
        Argument is index=index to load part of the variable.'''
        lon = numpy.memmap(p.indatadir+'XC.data','>f4',shape=(p.Ny,p.Nx),mode='r')
        lat = numpy.memmap(p.indatadir+'YC.data','>f4',shape=(p.Ny,p.Nx),mode='r')
        self.vlat=lat[:,:]
        self.vlon=(lon[:,:]+360)%360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from netcdf file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if numpy.min(self.vlon)<1. and numpy.max(self.vlon)>359.:
            self.vlon[numpy.where(self.vlon>180.)]=self.vlon[numpy.where(self.vlon>180.)]-360
            lon1=(numpy.min(self.vlon)+360)%360
            lon2=(numpy.max(self.vlon)+360)%360
        else:
            lon1=numpy.min(self.vlon)
            lon2=numpy.max(self.vlon)
        return [lon1, lon2, numpy.min(self.vlat), numpy.max(self.vlat)]
