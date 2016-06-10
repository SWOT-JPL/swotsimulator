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
between model files and thus if several files are used in the SSH interpolation,
some discontinuities may be visible. \n

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

def run_nadir(file_param):

    if os.path.isfile(file_param):
#    basedir=os.path.dirname(swotsimulator.__file__)
        shutil.copyfile(file_param, 'params.py') #os.path.join(basedir,'params.py'))
    else: 
        print("Error: No such file: '%s'" % file_param)
        sys.exit()
    try: 
        import params as p
    except:
        if os.path.isfile('params.py'):
            print("There is a wrong entry in your params file")
            import params
        else: 
            print("Error: No params.py module found")
            sys.exit()
    import swotsimulator.build_swath as build_swath
    import swotsimulator.rw_data as rw_data
    import swotsimulator.build_error as build_error
    import swotsimulator.mod_tools as mod_tools

## - Initialize some parameters values
    try: p.shift_lon=p.shift_lon
    except: p.shift_lon=None
    try: p.shift_time=p.shift_time
    except: p.shift_time=None
    try: model=p.model
    except: model='NETCDF_MODEL'; p.model= model
    try: p.model_nan=p.model_nan
    except: p.model_nan=0.
    try:   p.SSH_factor=p.SSH_factor
    except: p.SSH_factor=1. #; p.SSH_factor=SSH_factor
    p.halfswath=0.


## - Read list of user model files """
    if p.file_input:
        list_file = [line.strip() for line in open(p.file_input)]

## - Read model input coordinates """
    if p.file_input: model_data= eval('rw_data.'+model+'(file=p.indatadir+os.sep+list_file[0])')
    if  p.modelbox:
        modelbox=numpy.array(p.modelbox, dtype='float')
        if modelbox[0]<0: modelbox[0]=modelbox[0]+360
        if modelbox[1]<0:modelbox[1]=modelbox[1]+360
    else: 
        try: modelbox=model_data.calc_box()
        except: 
            print('Modelbox could not be computed, a modelbox should be given if no model file is provided')
            sys.exit()
    if p.file_input: 
        model_data.read_coordinates()
        model_index=numpy.where(((modelbox[0]-1)<=model_data.vlon) & (model_data.vlon<=(modelbox[1]+1)) & ((modelbox[2]-1)<=model_data.vlat) & (model_data.vlat<=(modelbox[3]+1)))#[0]


## - Initialize random coefficients that are use to compute 
##   random errors following the specified spectrum
    errnad=build_error.errornadir(p)
    if p.file_coeff:
        if os.path.isfile(p.file_coeff[:-3]+'_nadir.nc') and (not p.makesgrid):
            print('\n WARNING: old random coefficient file are used')
            errnad.load_coeff(p)
        else:
            errnad.init_error(p)
            errnad.save_coeff(p)
    else: 
        errnad.init_error(p)
## - Compute interpolated SSH and errors for each pass, at each
##   cycle
    print('Compute interpolated SSH and errors:')
##   load all SWOT grid files (one for each pass)
##   model time step
    modeltime=numpy.arange(0,p.nstep*p.timestep, p.timestep) 
##   remove the grid from the list of model files
    if p.file_input: list_file.remove(list_file[0])
##   initialize progress bar variables
    istep=0; ntot=1
## - Loop on SWOT grid files
 # for ifile in p.filsat:
    istring=len(p.dir_setup)
    if not isinstance(p.filesat, list): 
        p.filesat=[p.filesat]
    for filesat in p.filesat:
##   load SWOT grid file
        ntmp, nfilesat=os.path.split(filesat[istring:-4])
        if p.makesgrid:
            print('\n Force creation of satellite grid')
            orb=build_swath.makeorbit(modelbox, p, orbitfile=filesat) 
            orb.file=p.outdatadir+os.sep+nfilesat+'_grid.nc'
            orb.write_orb()
        else:
            orb = rw_data.Sat_nadir(file=p.filesgrid+'.nc')
            cycle=0 ; x_al=[] ;  al_cycle=0 ; timeshift=0
            orb.load_orb(cycle=cycle,x_al=x_al, al_cycle=al_cycle, timeshift=timeshift)
        cycle=orb.cycle ; x_al=orb.x_al ; al_cycle=orb.al_cycle ; timeshift=orb.timeshift
##   Look for model indices corresponding to the localization of the pass
##   If the pass is not in the model region, continue to next loop
        if numpy.min(orb.lon)<1. and numpy.max(orb.lon)>359.:
            orb.lon[numpy.where(orb.lon>180.)]=orb.lon[numpy.where(orb.lon>180.)]-360
            if p.file_input: model_data.vlon[numpy.where(model_data.vlon>180.)]=model_data.vlon[numpy.where(model_data.vlon>180.)]-360
            if modelbox[0]>180.: modelbox[0]=modelbox[0]-360
            if modelbox[1]>180.: modelbox[1]=modelbox[1]-360
        if p.file_input: 
            model_index=numpy.where(((numpy.min(orb.lon))<=model_data.vlon) & (model_data.vlon<=(numpy.max(orb.lon))) & ((numpy.min(orb.lat))<=model_data.vlat) & (model_data.vlat<=(numpy.max(orb.lat))))
            lonmodel1d=model_data.vlon[model_index].ravel()
            latmodel1d=model_data.vlat[model_index].ravel()
            if not model_index[0].any():
                continue
##   Compute number of cycles needed to cover all nstep model timesteps
        rcycle=(p.timestep*p.nstep)/float(orb.cycle)
        ncycle=int(rcycle)
##   Switch to (-180,180) longitude range to interpolate near 360-0

## - Loop on all cycles
        for cycle in range(0,ncycle+1):
##   Initialiaze errors and SSH
            errnad.nadir=numpy.zeros((numpy.shape(orb.x_al)[0]))
            errnad.wet_tropo1=numpy.zeros((numpy.shape(orb.x_al)[0]))
            errnad.wt=numpy.zeros((numpy.shape(orb.x_al)[0]))
            date1=cycle*orb.cycle
            SSH_true=numpy.zeros((numpy.shape(orb.x_al)[0])) #, p.nstep))
            #SSH_obs=numpy.empty((numpy.shape(orb.x_al)[0]))
            vindice=numpy.zeros(numpy.shape(SSH_true))
## Definition of the time in the model 
            stime=orb.time+date1
## Look for satellite data that are beween step-p.timesetp/2 end setp+p.step/2
            if p.file_input: 
                index=numpy.where(((stime[-1]-timeshift)>=(modeltime-p.timestep/2.)) & ((stime[0]-timeshift)<(modeltime+p.timestep/2.)) )#[0]
## At each step, look for the corresponding time in the satellite data
                for ifile in index[0]:
                    ntot=1
                    mod_tools.update_progress(float(istep)/float(p.nstep*ntot*numpy.shape(p.filesat)[0]), 'Sat: '+ filesat[istring:], 'model file: '+list_file[ifile] +', cycle:'+str(cycle+1))
                    #mod_tools.update_progress(float(istep)/float(p.nstep*ntot), 'Sat: '+ filesat[istring:], 'model file: '+list_file[ifile] +', cycle:'+str(cycle+1))
## If there are satellite data, Get true SSH from model
                    if numpy.shape(index)[1]>0:
                        ntot=ntot+numpy.shape(index)[1]-1
                        #if numpy.shape(index)[1]>1:
                        ind=numpy.where(((stime-timeshift)>=(modeltime[ifile]-p.timestep/2.)) & ((stime-timeshift)<(modeltime[ifile]+p.timestep/2.)) )
                        exec('model_step=rw_data.'+model+'(file=p.indatadir+os.sep+list_file[ifile], var=p.var)')
                        model_step.read_var(index=model_index)
                        SSH_model=model_step.vvar
                        SSH_true[ind[0]]=interpolate.griddata((lonmodel1d, latmodel1d), SSH_model.ravel(), (orb.lon[ind[0]], orb.lat[ind[0]]), method=p.interpolation)
                        if p.interpolation=='nearest':
                            if modelbox[0]>modelbox[1]:
                                SSH_true[numpy.where(((orb.lon<modelbox[0]) & (orb.lon>modelbox[1])) | (orb.lat<modelbox[2]) | (orb.lat>modelbox[3]))] = numpy.nan 
                            else:
                                SSH_true[numpy.where((orb.lon<modelbox[0]) | (orb.lon>modelbox[1]) | (orb.lat<modelbox[2]) | (orb.lat>modelbox[3]))] = numpy.nan 
                        #SSH_true[numpy.where((orb.lon<modelbox[0]) | (orb.lon>modelbox[1]) | (orb.lat<modelbox[2]) | (orb.lat>modelbox[3]))] = numpy.nan 
                        vindice[ind[0]]=ifile
                        del ind, SSH_model, model_step
                    istep+=1
            else:
                istep+=1
                mod_tools.update_progress(float(istep)/float(rcycle*numpy.shape(p.filesat)[0]), 'Sat: '+ filesat[istring:], 'no model file provide '+', cycle:'+str(cycle+1))

            SSH_true[SSH_true==0]=numpy.nan
            errnad.make_error(orb, cycle, SSH_true, p) #, ind[0])
            if p.wet_tropo: errnad.SSH=SSH_true+errnad.nadir+errnad.wet_tropo1 #make_SSH_errornadir(SSH_true,p)
            else: errnad.SSH=SSH_true+errnad.nadir #make_SSH_errornadir(SSH_true,p)
## Save outputs in a netcdf file
            if vindice.any() or not p.file_input:
                file_output=p.file_output+'_c'+str(cycle+1).zfill(2)+'.nc'
                Output=rw_data.Sat_nadir(file=file_output, lon=(orb.lon+360)%360, lat=orb.lat, time=stime, x_al=orb.x_al, cycle=orb.cycle)
                Output.write_data(SSH_model=SSH_true, index=vindice, nadir_err=errnad.nadir, SSH_obs=errnad.SSH, pd_err_1b=errnad.wet_tropo1,pd=errnad.wt) #, ncycle=cycle)
        del stime
        if p.file_input: del index
        #orb.lon=(orb.lon+360)%360
        #if p.file_input: model_data.vlon=(model_data.vlon+360)%360
        #modelbox[0]=(modelbox[0]+360)%360
        #modelbox[1]=(modelbox[1]+360)%360
        del orb
## - Write Selected parameters in a txt file
    rw_data.write_params(p,p.outdatadir+os.sep+'orbit_simulator.output')
    print("\n Simulated orbit files have been written in " + p.outdatadir)
    print("----------------------------------------------------------")
