'''Modules to run the main program: 
Contains the following functions: 
load_error(param_file) : load or compute random coefficients to compute errors
\n
load_sgrid(sgridfile, p): load or compute SWOT grid
\n
load_ngrid(sgridfile, p): load or compute nadir grid
\n
select_modelbox(sgrid, model_data): Select region of interest with params file or model data
\n
create_SWOTlikedata(cycle, ntotfile, list_file, sgrid, ngrid, model_data, modeltime, err, errnad, p):
interpolate data on SWOT swath 
\n
CreateNadirlikedata(cycle, ntotfile, list_file, ngrid, model_data, modeltime,  errnad, p, progress_bar=False):
Interpolate data on nadir track
\n
def save_SWOT(cycle, sgrid, err, p, time=[], vindice=[], SSH_true=[]):
Save SWOT-like data
\n
save_Nadir(cycle, ngrid, errnad, err , p, time=[], vindice_nadir=[], SSH_true_nadir=[]):
Save Nadir-like data
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
#Copyright (c) 2002-2014, California Institute of Technology.
#All rights reserved. Based on Government Sponsored Research under contracts NAS7-1407 and/or NAS7-03001.
#
#-----------------------------------------------------------------------
'''
import os, shutil
from scipy import interpolate
import numpy
import glob, gc
from optparse import OptionParser
import sys
import time as ti
import swotsimulator
import swotsimulator.build_swath as build_swath
import swotsimulator.rw_data as rw_data
import swotsimulator.build_error as build_error
import swotsimulator.mod_tools as mod_tools
import swotsimulator.const as const

## - Define global variables for progress bars
istep=0
ntot=1
ifile=0

def load_error(p):
    '''Initialize random coefficients that are used to compute
    random errors following the specified spectrum. \n
    If a random coefficient file is specified, random coefficients
    are loaded from this file.
    '''
    import swotsimulator.build_error as build_error
    err=build_error.error(p)
    if p.nadir:
        errnad=build_error.errornadir(p)
    nhalfswath=int( (p.halfswath-p.halfgap)/p.delta_ac) +1
    if p.file_coeff:
        if os.path.isfile(p.file_coeff) and (not p.makesgrid):
            print('\n WARNING: Existing random coefficient file used')
            err.load_coeff(p)
        else:
            err.init_error(p, 2*nhalfswath)
            err.save_coeff(p,2*nhalfswath)
        if p.nadir:
            if os.path.isfile(p.file_coeff[:-3]+'_nadir.nc')  and (not p.makesgrid):
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

    ## Load SWOT swath file
    sgrid = rw_data.Sat_SWOT(file=sgridfile)
    cycle=0 ; x_al=[] ; x_ac=[] ; al_cycle=0 ; timeshift=0
    sgrid.load_swath(cycle=cycle,x_al=x_al, x_ac=x_ac, al_cycle=al_cycle, timeshift=timeshift)
    sgrid.loncirc=numpy.rad2deg(numpy.unwrap(sgrid.lon))
    ## Extract the pass number from the file name
    ipass=int(sgridfile[-6:-3])
    sgrid.ipass=ipass
    return sgrid

def load_ngrid(sgridfile, p):
    ipass=int(sgridfile[-6:-3])
    ## Load Nadir track file
    ngrid = rw_data.Sat_nadir(file=p.filesgrid+'nadir_p'+str(ipass).zfill(3)+'.nc')
    cycle=0 ; x_al=[] ; al_cycle=0 ; timeshift=0
    ngrid.load_orb(cycle=cycle,x_al=x_al, al_cycle=al_cycle, timeshift=timeshift)
    ngrid.loncirc=numpy.rad2deg(numpy.unwrap(ngrid.lon))
    ngrid.ipass=ipass
    return ngrid

def select_modelbox(sgrid, model_data):
    #mask=numpy.zeros((numpy.shape(model_data.vlon)))
    #nal=len(sgrid.x_al)
    #for kk in range(0,nal,10):
    #  dlon1=model_data.vlon-sgrid.lon_nadir[kk]
    #  dlon=numpy.minimum(numpy.mod(dlon1,360),numpy.mod(-dlon1,360))
    #  ddist=(((dlon)*numpy.cos(sgrid.lat_nadir[kk]*numpy.pi/180.)*const.deg2km)**2 + ((model_data.vlat-sgrid.lat_nadir[kk])*const.deg2km)**2 )
    #  mask[ddist<const.radius_interp**2]=1
    model_index=numpy.where(((numpy.min(sgrid.lon))<=model_data.vlon) & (model_data.vlon<=(numpy.max(sgrid.lon))) & ((numpy.min(sgrid.lat))<=model_data.vlat) & (model_data.vlat<=(numpy.max(sgrid.lat))))
    model_data.model_index=model_index
    lonmodel1d=model_data.vlon[model_index].ravel()
    latmodel1d=model_data.vlat[model_index].ravel()

    #pdb.set_trace()
    if p.grid=='regular':
      model_data.lon1d=lon1[numpy.where(((numpy.min(sgrid.lon))<=lon1) & (lon1<=(numpy.max(sgrid.lon))))]
      model_data.lat1d=lat1[numpy.where(((numpy.min(sgrid.lat))<=lat1) & (lat1<=(numpy.max(sgrid.lat))))]
    else:
        model_index=numpy.where(((numpy.min(sgrid.lon))<=model_data.vlon) & (model_data.vlon<=(numpy.max(sgrid.lon))) & ((numpy.min(sgrid.lat))<=model_data.vlat) & (model_data.vlat<=(numpy.max(sgrid.lat))))
        model_data.lon1d=model_data.vlon[model_index].ravel()
        model_data.lat1d=model_data.vlat[model_index].ravel()
        model_data.vlon=model_data.vlon[model_index].ravel()
        model_data.vlat=model_data.vlat[model_index].ravel()

    #pdb.set_trace()
    nx=len(lon1)
    ny=len(lat1)
    return model_index


def create_SWOTlikedata(cycle, ntotfile, list_file, sgrid, ngrid, model_data, modeltime, err, errnad, p):
    import swotsimulator.rw_data as rw_data
    import swotsimulator.build_error as build_error
    import swotsimulator.mod_tools as mod_tools
    import swotsimulator.const as const
    '''Create SWOT and nadir errors err and errnad, interpolate model SSH model_data on swath and nadir track,
    compute SWOT-like and nadir-like data for cycle, SWOT grid sgrid and ngrid. '''
## - Progress bar variables are global
    global istep
    global ntot
##   Initialiaze errors and SSH
    progress=0
    err.karin=numpy.zeros((numpy.shape(sgrid.lon)[0], numpy.shape(sgrid.lon)[1]))
    err.roll=numpy.zeros((numpy.shape(sgrid.lon)[0], numpy.shape(sgrid.lon)[1]))
    err.phase=numpy.zeros((numpy.shape(sgrid.lon)[0], numpy.shape(sgrid.lon)[1]))
    err.baseline_dilation=numpy.zeros((numpy.shape(sgrid.lon)[0], numpy.shape(sgrid.lon)[1]))
    err.timing=numpy.zeros((numpy.shape(sgrid.lon)[0], numpy.shape(sgrid.lon)[1]))
    err.wet_tropo1=numpy.zeros((numpy.shape(sgrid.lon)[0], numpy.shape(sgrid.lon)[1]))
    err.wet_tropo2=numpy.zeros((numpy.shape(sgrid.lon)[0], numpy.shape(sgrid.lon)[1]))
    err.ssb=numpy.zeros((numpy.shape(sgrid.lon)[0], numpy.shape(sgrid.lon)[1]))
    err.wt=numpy.zeros((numpy.shape(sgrid.lon)[0], numpy.shape(sgrid.lon)[1]))
    SSH_true=numpy.zeros((numpy.shape(sgrid.lon)[0], numpy.shape(sgrid.lon)[1])) #, p.nstep))
    if p.nadir:
         err.wet_tropo1nadir=numpy.zeros((numpy.shape(ngrid.lon)[0]))
         err.wet_tropo2nadir=numpy.zeros((numpy.shape(ngrid.lon)[0]))
         err.wtnadir=numpy.zeros((numpy.shape(ngrid.lon)[0]))
         errnad.nadir=numpy.zeros((numpy.shape(ngrid.lon)[0]))
         SSH_true_nadir=numpy.zeros((numpy.shape(ngrid.lon)[0]))
         vindice_nadir=numpy.zeros(numpy.shape(ngrid.lon))*numpy.nan
    date1=cycle*sgrid.cycle
    vindice=numpy.zeros(numpy.shape(SSH_true))*numpy.nan
## Definition of the time in the model
    time=sgrid.time+date1
## Look for satellite data that are beween step-p.timesetp/2 end setp+p.step/2
    if p.file_input:
        index_filemodel=numpy.where(((time[-1]-sgrid.timeshift)>=(modeltime-p.timestep/2.)) & ((time[0]-sgrid.timeshift)<(modeltime+p.timestep/2.)) )#[0]
        ## At each step, look for the corresponding time in the satellite data
        for ifile in index_filemodel[0]:
          progress=mod_tools.update_progress(float(istep)/float(ntot*ntotfile), 'pass: '+str(sgrid.ipass), 'model file: '+list_file[ifile] +', cycle:'+str(cycle+1))
          ## If there are satellite data, Get true SSH from model
          if numpy.shape(index_filemodel)[1]>0:
            ## number of file to be processed used in the progress bar
            ntot=ntot+numpy.shape(index_filemodel)[1]-1
            #if numpy.shape(index)[1]>1:
            ## Select part of the track that corresponds to the time of the model (+-timestep/2)
            ind_time=numpy.where(((time-sgrid.timeshift)>=(modeltime[ifile]-p.timestep/2.)) & ((time-sgrid.timeshift)<(modeltime[ifile]+p.timestep/2.)) )
            if p.nadir: ind_nadir_time=numpy.where(((time-ngrid.timeshift)>=(modeltime[ifile]-p.timestep/2.)) & ((time-ngrid.timeshift)<(modeltime[ifile]+p.timestep/2.)) )
            ## Load data from this model file
            model_step=eval('rw_data.'+model_data.model+'(file=p.indatadir+os.sep+list_file[ifile], var=p.var)')
            if p.grid=='regular':
              model_step.read_var()
              SSH_model=model_step.vvar[model_data.model_index_lat, :]
              SSH_model=SSH_model[:,model_data.model_index_lon]
            else: 
              model_step.read_var(index=model_data.model_index)
              SSH_model=model_step.vvar
## - Interpolate Model data on a SWOT grid and/or along the nadir track
            ## if grid is regular, use interpolate.RectBivariateSpline to interpolate
            if p.grid=='regular' and len(numpy.shape(model_data.vlon))==1:
            #########################TODO
              #To be moved to routine rw_data
              indsorted=numpy.argsort(model_data.vlon)
              model_data.vlon=model_data.vlon[indsorted]
              SSH_model=SSH_model[:,indsorted]
                    
              ## Flatten satellite grid and select part of the track corresponding to the model time
              lonswot=sgrid.lon[ind_time[0],:].flatten()
              lonnadir=ngrid.lon[ind_nadir_time[0]].flatten()
              latswot=sgrid.lat[ind_time[0],:].flatten()
              latnadir=ngrid.lat[ind_nadir_time[0]].flatten()
              Teval=interpolate.RectBivariateSpline(model_data.vlat,model_data.vlon,numpy.isnan(SSH_model), kx=1, ky=1, s=0).ev(sgrid.lat[ind_time[0],:].ravel(),sgrid.lon[ind_time[0],:].ravel())
              SSH_model_mask=+SSH_model
              SSH_model_mask[numpy.isnan(SSH_model_mask)]=0.
              SSH_true_ind_time=interpolate.RectBivariateSpline(model_data.vlat,model_data.vlon,SSH_model_mask, kx=1, ky=1, s=0).ev(sgrid.lat[ind_time[0],:].ravel(),sgrid.lon[ind_time[0],:].ravel())
              SSH_true_ind_time[Teval>0]=numpy.nan
              nal,nac=numpy.shape(sgrid.lon[ind_time[0],:])
              SSH_true[ind_time[0], :]=SSH_true_ind_time.reshape(nal,nac)
              if p.nadir:
                Teval=interpolate.RectBivariateSpline(model_data.vlat,model_data.vlon,numpy.isnan(SSH_model), kx=1, ky=1, s=0).ev(ngrid.lon[ind_nadir_time[0]].ravel(),ngrid.lat[ind_nadir_time[0]].ravel())
                SSH_model_mask=+SSH_model
                SSH_model_mask[numpy.isnan(SSH_model_mask)]=0.
                SSH_true_ind_time=interpolate.RectBivariateSpline(model_data.vlat,model_data.vlon,SSH_model_mask, kx=1, ky=1, s=0).ev(ngrid.lon[ind_nadir_time[0]].ravel(),ngrid.lat[ind_nadir_time[0]].ravel())
                SSH_true_ind_time[Teval>0]=numpy.nan
                SSH_true_nadir[ind_nadir_time[0]]=SSH_true_ind_time
            else:
            ## Grid is irregular, interpolation can be done using pyresample module if it is installed or griddata function from scipy. 
            ## Note that griddata is slower than pyresample functions. 
              try: 
                import pyresample as pr
                model_data.vlon=pr.utils.wrap_longitudes(model_data.vlon)
                sgrid.lon=pr.utils.wrap_longitudes(sgrid.lon)
                if len(numpy.shape(model_data.vlon))<=1:
                  model_data.vlon, model_data.vlat=numpy.meshgrid(model_data.vlon, model_data.vlat)
                swath_def=pr.geometry.SwathDefinition(lons=model_data.vlon, lats=model_data.vlat)
                grid_def=pr.geometry.SwathDefinition(lons=sgrid.lon, lats=sgrid.lat)
                if p.interpolation=='nearest':
                 SSH_true[ind_time[0], :]=pr.kd_tree.resample_nearest(swath_def, SSH_model,grid_def, radius_of_influence=max(p.delta_al, p.delta_ac)*10**3, epsilon=100)
                else:
                 SSH_true[ind_time[0], :]=pr.kd_tree.resample_gauss(swath_def, SSH_model, grid_def, radius_of_influence=3*max(p.delta_al, p.delta_ac)*10**3,sigmas=max(p.delta_al, p.delta_ac)*10**3,fill_value=None)
                if p.nadir:
                  ngrid.lon=pr.utils.wrap_longitudes(ngrid.lon)
                  ngrid_def=pr.geometry.SwathDefinition(lons=ngrid.lon, lats=ngrid.lat)
                  if p.interpolation=='nearest':
                    SSH_true_nadir[ind_nadir_time[0]]=pr.kd_tree.resample_nearest(swath_def, SSH_model,ngrid_def, radius_of_influence=max(p.delta_al, p.delta_ac)*10**3, epsilon=100)
                  else:
                    SSH_true_nadir[ind_nadir_time[0]]=pr.kd_tree.resample_gauss(swath_def, SSH_model, ngrid_def, radius_of_influence=3*max(p.delta_al, p.delta_ac)*10**3,sigmas=max(p.delta_al, p.delta_ac)*10**3,fill_value=None) 
              except: 
                SSH_true[ind_time[0], :]=interpolate.griddata((model_data.vlon.ravel(), model_data.vlat.ravel()), SSH_model.ravel(), (sgrid.lon[ind_time[0],:], sgrid.lat[ind_time[0],:]), method=p.interpolation)
                if p.nadir: SSH_true_nadir[ind_nadir_time[0]]=interpolate.griddata((model_data.vlon.ravel(), model_data.vlat.ravel()), SSH_model.ravel(), (ngrid.lon[ind_nadir_time[0]], ngrid.lat[ind_nadir_time[0]]), method=p.interpolation)
                if p.interpolation=='nearest':
                  if modelbox[0]>modelbox[1]:
                    SSH_true[numpy.where(((sgrid.lon<modelbox[0]) & (sgrid.lon>modelbox[1])) | (sgrid.lat<modelbox[2]) | (sgrid.lat>modelbox[3]))] = numpy.nan
                    if p.nadir: SSH_true_nadir[numpy.where(((ngrid.lon<modelbox[0]) & (ngrid.lon>modelbox[1])) | (ngrid.lat<modelbox[2]) | (ngrid.lat>modelbox[3]))] = numpy.nan
                  else:
                    SSH_true[numpy.where((sgrid.lon<modelbox[0]) | (sgrid.lon>modelbox[1]) | (sgrid.lat<modelbox[2]) | (sgrid.lat>modelbox[3]))] = numpy.nan
                    if p.nadir: SSH_true_nadir[numpy.where((ngrid.lon<modelbox[0]) | (ngrid.lon>modelbox[1]) | (ngrid.lat<modelbox[2]) | (ngrid.lat>modelbox[3]))] = numpy.nan
            vindice[ind_time[0], :]=ifile
            if p.nadir: vindice_nadir[ind_nadir_time[0]]=ifile
            del ind_time, SSH_model, model_step, ind_nadir_time
          istep+=1
    else:
        istep+=1
        progress=mod_tools.update_progress(float(istep)/float(ntotfile*ntot), 'pass: '+str(sgrid.ipass), 'no model file provided'+', cycle:'+str(cycle+1))
    err.make_error(sgrid, cycle, SSH_true,p) #, ind[0])
    err.make_SSH_error(SSH_true,p)
    if p.nadir:
        errnad.make_error(ngrid, cycle, SSH_true_nadir,p) #, ind[0])
        if p.nbeam==1: errnad.SSH=SSH_true_nadir+errnad.nadir+err.wet_tropo1nadir
        else: errnad.SSH=SSH_true_nadir+errnad.nadir+err.wet_tropo2nadir
    #if p.file_input: del ind_time, SSH_model, model_step
    return SSH_true, SSH_true_nadir, vindice, vindice_nadir, time, progress

def createNadirlikedata(cycle, ntotfile, list_file, ngrid, model_data, modeltime,  errnad, p, progress_bar=False):
    import swotsimulator.rw_data as rw_data
    import swotsimulator.build_error as build_error
    import swotsimulator.mod_tools as mod_tools
    import swotsimulator.const as const

## - Progress bar variables are global
    global istep
    global ntot
    global ifile
    errnad.wet_tropo1nadir=numpy.zeros((numpy.shape(ngrid.lon)[0]))
    errnad.wt=numpy.zeros((numpy.shape(ngrid.lon)[0]))
    errnad.nadir=numpy.zeros((numpy.shape(ngrid.lon)[0]))
    SSH_true=numpy.zeros((numpy.shape(ngrid.lon)[0]))
    vindice=numpy.zeros(numpy.shape(ngrid.lon))*numpy.nan
    ## Definition of the time in the model
    time=ngrid.time+date1
    ## Look for satellite data that are beween step-p.timesetp/2 end setp+p.step/2
    if p.file_input:
        index_filemodel=numpy.where(((time[-1]-ngrid.timeshift)>=(modeltime-p.timestep/2.)) & ((time[0]-ngrid.timeshift)<(modeltime+p.timestep/2.)) )#[0]
        ## At each step, look for the corresponding time in the satellite data
        for ifile in index_filemodel[0]:
            if progress_bar:
                progress=mod_tools.update_progress(float(istep)/float(ntotfile*ntot), 'pass: '+str(ngrid.ipass), 'model file: '+list_file[ifile] +', cycle:'+str(cycle+1))
                ## If there are satellite data, Get true SSH from model
            if numpy.shape(index_filemodel)[1]>0:
                ntot=ntot+numpy.shape(index_filemodel)[1]-1
                ind_nadir_time=numpy.where(((time-ngrid.timeshift)>=(modeltime[ifile]-p.timestep/2.)) & ((time-ngrid.timeshift)<(modeltime[ifile]+p.timestep/2.)) )
                model_step=eval('rw_data.'+model+'(file=p.indatadir+os.sep+list_file[ifile], var=p.var)')
                if p.grid=='regular':
                  model_step.read_var()
                  SSH_model=model_step.vvar[model_data.model_index_lat, :]
                  SSH_model=SSH_model[:,model_data.model_index_lon]
                else:
                  model_step.read_var(index=model_data.model_index)
                  SSH_model=model_step.vvar
##########REWRITE INTERPOLATION
            #############CLEM###################
            #pdb.set_trace()
            #SSH_model=model_step.vvar.reshape(nx,ny)
## - Interpolate Model data on a SWOT grid and/or along the nadir track
            ## if grid is regular, use interpolate.RectBivariateSpline to interpolate
            if p.grid=='regular' and len(numpy.shape(model_data.vlon))==1:
            #########################TODO
              #To be moved to routine rw_data
              indsorted=numpy.argsort(model_data.vlon)
              model_data.vlon=model_data.vlon[indsorted]
              SSH_model=SSH_model[:,indsorted]

              lonnadir=ngrid.lon[ind_nadir_time[0]].flatten()
              latnadir=ngrid.lat[ind_nadir_time[0]].flatten()

              Teval=interpolate.RectBivariateSpline(model_data.vlat,model_data.vlon,numpy.isnan(SSH_model), kx=1, ky=1, s=0).ev(ngrid.lon[ind_nadir_time[0]].ravel(),ngrid.lat[ind_nadir_time[0]].ravel())
              SSH_model_mask=+SSH_model
              SSH_model_mask[numpy.isnan(SSH_model_mask)]=0.
              SSH_true_ind_time=interpolate.RectBivariateSpline(model_data.vlat,model_data.vlon,SSH_model_mask, kx=1, ky=1, s=0).ev(ngrid.lon[ind_nadir_time[0]].ravel(),ngrid.lat[ind_nadir_time[0]].ravel())
              SSH_true_ind_time[Teval>0]=numpy.nan
              SSH_true_nadir[ind_nadir_time[0]]=SSH_true_ind_time
############MARKPOINT
            else:
            ## Grid is irregular, interpolation can be done using pyresample module if it is installed or griddata function from scipy.
            ## Note that griddata is slower than pyresample functions.
              try:
                import pyresample as pr
                ngrid.lon=pr.utils.wrap_longitudes(ngrid.lon)
                ngrid_def=pr.geometry.SwathDefinition(lons=ngrid.lon, lats=ngrid.lat)
                if p.interpolation=='nearest':
                  SSH_true_nadir[ind_nadir_time[0]]=pr.kd_tree.resample_nearest(swath_def, SSH_model,ngrid_def, radius_of_influence=max(p.delta_al, p.delta_ac)*10**3, epsilon=100)
                else:
                  SSH_true_nadir[ind_nadir_time[0]]=pr.kd_tree.resample_gauss(swath_def, SSH_model, ngrid_def, radius_of_influence=3*max(p.delta_al, p.delta_ac)*10**3,sigmas=max(p.delta_al, p.delta_ac)*10**3,fill_value=None)
              except:
                SSH_true_nadir[ind_nadir_time[0]]=interpolate.griddata((model_data.vlon.ravel(), model_data.vlat.ravel()), SSH_model.ravel(), (ngrid.lon[ind_nadir_time[0]], ngrid.lat[ind_nadir_time[0]]), method=p.interpolation)
                if p.interpolation=='nearest':
                  if modelbox[0]>modelbox[1]:
                    SSH_true_nadir[numpy.where(((ngrid.lon<modelbox[0]) & (ngrid.lon>modelbox[1])) | (ngrid.lat<modelbox[2]) | (ngrid.lat>modelbox[3]))] = numpy.nan
                  else:
                    SSH_true_nadir[numpy.where((ngrid.lon<modelbox[0]) | (ngrid.lon>modelbox[1]) | (ngrid.lat<modelbox[2]) | (ngrid.lat>modelbox[3]))] = numpy.nan
            vindice_nadir[ind_nadir_time[0]]=ifile
    errnad.make_error(ngrid, cycle, SSH_true_nadir,p) #, ind[0])
    errnad.SSH=errnad.nadir+err.wet_tropo1nadir
    del SSH_model, model_step, ind_nadir_time
    return SSH_true_nadir, vindice, time, progress

#def interp_regularmodel():

#def interp_irregularmodel():


def save_SWOT(cycle, sgrid, err, p, time=[], vindice=[], SSH_true=[]):
  file_output=p.file_output+'_c'+str(cycle+1).zfill(2)+'_p'+str(sgrid.ipass).zfill(3)+'.nc'
  OutputSWOT=rw_data.Sat_SWOT(file=file_output, lon=(sgrid.lon+360)%360, lat=sgrid.lat, time=time, x_ac=sgrid.x_ac, x_al=sgrid.x_al, cycle=sgrid.cycle, lon_nadir=(sgrid.lon_nadir+360)%360, lat_nadir=sgrid.lat_nadir)
  OutputSWOT.write_data(SSH_model=SSH_true, index=vindice, roll_err=err.roll, bd_err=err.baseline_dilation, phase_err=err.phase, ssb_err=err.ssb, karin_err=err.karin, pd_err_1b=err.wet_tropo1, pd_err_2b=err.wet_tropo2, pd=err.wt, timing_err=err.timing, SSH_obs=err.SSH)
  return None

def save_Nadir(cycle, ngrid, errnad, err , p, time=[], vindice_nadir=[], SSH_true_nadir=[]):
    file_outputnadir=p.file_output+'nadir_c'+str(cycle+1).zfill(2)+'_p'+str(ngrid.ipass).zfill(3)+'.nc'
    OutputNadir=rw_data.Sat_nadir(file=file_outputnadir, lon=(ngrid.lon+360)%360, lat=ngrid.lat, time=time, x_al=ngrid.x_al, cycle=ngrid.cycle)
    OutputNadir.write_data(SSH_model=SSH_true_nadir, index=vindice_nadir, nadir_err=errnad.nadir, SSH_obs=errnad.SSH, pd_err_1b=err.wet_tropo1nadir,pd=err.wtnadir, pd_err_2b=err.wet_tropo2nadir)
    return None
    


