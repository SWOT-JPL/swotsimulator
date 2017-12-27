import numpy
import math
from scipy import interpolate
import swotsimulator.mod_tools as mod_tools
import swotsimulator.const as const
import swotsimulator.rw_data as rw_data
import os
import logging
import sys

# Define logger level for debug purposes
logger = logging.getLogger(__name__)


def makeorbit(modelbox, p, orbitfile='orbit_292.txt', filealtimeter=None):
    '''Computes the orbit nadir on a subdomain.
    The path of the satellite is given by the orbit file and the subdomain
    corresponds to the one in the model. Note that a subdomain can be manually
    added in the parameters file. \n
    Inputs are satellite orbit (p.filesat), subdomain (modelbox), Along track
    sampling, along track resolution). \n
    Outputs are Sat_Nadir object containing Nadir track (along track distance
    x_al, longitude lon and latitude lat,
    number of days in a cycle cycle, distance crossed in a cycle cycle_al,
    time, time shfit and time of pass passtime'''
    # npoints = 1
    # - Load SWOT orbit ground track
    logger.info('Load data from orbit file')
    if p.order_orbit_col is None:
        volon, volat, votime = numpy.loadtxt(orbitfile, usecols=(1, 2, 0),
                                             unpack=True)

    else:
        ncols = p.order_orbit_col
        volon, volat, votime = numpy.loadtxt(orbitfile, usecols=ncols,
                                             unpack=True)
        votime *= const.secinday
    if (volon > 360).any() or (numpy.abs(volat) > 90).any():
        logger.error('Error in orbit file or wrong order of column \n'
                     'Columns should be in the following order'
                     '(time, lon, lat)')
        sys.exit(1)
    # - If orbit is at low resolution, interpolate at 0.5 s resolution
    if numpy.mean(votime[1:] - votime[:-1]) > 0.5:
        x, y, z = mod_tools.spher2cart(volon, volat)
        time_hr = numpy.arange(0., votime[-1], 0.5)
        f = interpolate.interp1d(votime, x)
        x_hr = f(time_hr)
        f = interpolate.interp1d(votime, y)
        y_hr = f(time_hr)
        f = interpolate.interp1d(votime, z)
        z_hr = f(time_hr)
        lon_hr = numpy.zeros(len(x_hr)) + numpy.nan
        lat_hr = numpy.zeros(len(x_hr)) + numpy.nan
        for ii in range(len(x_hr)):
            lon_hr[ii], lat_hr[ii] = mod_tools.cart2spher(x_hr[ii],
                                                          y_hr[ii],
                                                          z_hr[ii])
        time_hr = time_hr / const.secinday
        ind = numpy.where((time_hr < const.tcycle))
        volon = lon_hr[ind]
        volat = lat_hr[ind]
        votime = time_hr[ind]

    # - Get number of points in orbit
    nop = numpy.shape(votime)[0]

    # - Get cycle period.
    tcycle = votime[nop-1] + votime[1] - votime[0]
    # shift time if the user needs to shift the time of the orbit
    try:
        pshift_time = p.shift_time
        if pshift_time is not None:
            shift_index = numpy.where(votime >= pshift_time)[0]
            volon = numpy.hstack([volon[shift_index[0]:],
                                 volon[:shift_index[0]]])
            volat = numpy.hstack([volat[shift_index[0]:],
                                 volat[:shift_index[0]]])
    except:
        p.shift_time = None
    # shift lon if the user needs to shift the localisation of the orbit
    try:
        pshift_lon = p.shift_lon
        if pshift_lon is not None:
            volon = volon + pshift_lon
    except:
        p.shift_lon = None
    volon = (volon + 360) % 360

    # - Rearrange orbit starting from pass 1
    # Detect the beginning of pass 1 in orbit txt file. By definition, it is
    # the first passage at southernmost latitude.
    dlat = numpy.roll(volat, 1) - volat
    ind = numpy.where((dlat < 0) & (numpy.roll(dlat, 1) >= 0))
    # Shift coordinates, so that the first point of the orbit is the beginning
    # of pass 1
    decal = ind[0][-1]
    # timeshift = votime[-1] - votime[decal]
    volon = numpy.hstack([volon[decal:], volon[:decal]])
    volat = numpy.hstack([volat[decal:], volat[:decal]])
    votime = numpy.hstack([votime[decal:], votime[:decal]])
    votime = (votime - votime[0]) % tcycle
    if votime[numpy.where(votime < 0)]:
        logger.warn('WARNING: there are negative times in your orbit')
    del ind
    # Compute the initial time of each pass
    dlat = numpy.roll(volat, 1) - volat
    ind = numpy.where(((dlat < 0) & (numpy.roll(dlat, 1) >= 0)) | ((dlat > 0)
                      & (numpy.roll(dlat, 1) <= 0)))
    # index=numpy.hstack([0,ind[0]-1])
    index = ind[0]
    passtime = votime[index]  # Times of pass

    # - Compute accumulated along-track distance, longitude, latitude
    #   in the subdomain chosen by the user or given by the model

    # Extract points in the domain and count the number of points (nop) in the
    # subdomain
    # modelbox=[lonmin lonmax latmin latmax] add margin around the domain
    # plus one point to compute Satdir
    matnpbox = numpy.zeros((nop))
    if modelbox[0] > modelbox[1]:
        matnpbox[numpy.where((((modelbox[0] - p.halfswath / (const.deg2km
                 * math.cos(modelbox[2]*math.pi/180.))) <= volon)
                 | (volon <= (modelbox[1] + p.halfswath/(const.deg2km
                              * math.cos(modelbox[3]*math.pi/180.)))))
                 & ((modelbox[2] - p.halfswath/const.deg2km) <= volat)
                 & ((modelbox[3] + p.halfswath/const.deg2km) >= volat))] = 1
    else:
        matnpbox[numpy.where(((modelbox[0] - p.halfswath / (const.deg2km
                 * math.cos(modelbox[2]*math.pi/180.))) <= volon)
                 & (volon <= (modelbox[1] + p.halfswath/(const.deg2km
                              * math.cos(modelbox[3]*math.pi/180.))))
                 & ((modelbox[2] - p.halfswath/const.deg2km) <= volat)
                 & ((modelbox[3] + p.halfswath/const.deg2km) >= volat))] = 1
    norp = int(numpy.sum(matnpbox))
    # Initialize total distance travelled by the satellite since the first
    # point of the cycle in the subdomain at low (orbital file) resolution
    x_al_lr = numpy.zeros((norp))
    lon_lr = numpy.zeros((norp))
    lat_lr = numpy.zeros((norp))
    stime_lr = numpy.zeros((norp))

    # Initialize vector with accumulated distance travelled by the satellite
    indp = 0
    distance = numpy.zeros((nop))
    # Compute and store distances and coordinates that are in the defined
    # subdomain
    logger.info('Compute SWOT nadir coordinate in the new domain')
    for i in range(0, nop - 1):
        mod_tools.update_progress(float(i) / float(nop-1), None, None)
        if abs(volon[i + 1] - volon[i]) > 1:
            if volon[i + 1] > 180.:
                volon[i + 1] = volon[i + 1] - 360
            if volon[i] > 180.:
                volon[i] = volon[i] - 360
        distance[i+1] = (distance[i] + numpy.sqrt(((volon[i+1]-volon[i])
                         * const.deg2km*numpy.cos(volat[i+1]
                         * 2*math.pi/360.))**2 + ((volat[i+1] - volat[i])
                                                  * const.deg2km)**2))
        volon[i + 1] = (volon[i + 1] + 360) % 360
        if matnpbox[i]:
            x_al_lr[indp] = distance[i]
            lon_lr[indp] = (volon[i] + 360) % 360
            lat_lr[indp] = volat[i]
            stime_lr[indp] = votime[i]
            indp += 1

    # - Interpolate orbit at delta_al km resolution (default is delta_al=1)

    # Detect gap in time in stime (to detect step in x_al, lon and lat)
    dstime = stime_lr[:] - numpy.roll(stime_lr[:], 1)
    ind = numpy.where(dstime > 3*(votime[1] - votime[0]))
    index = numpy.hstack([0, ind[0]])
    nindex = numpy.shape(index)[0]
    # Initialize along track distance, time and coordinates at delta_al
    # resolution
    if nindex > 1:
        dgap = numpy.zeros((nindex))
        for i in range(1, nindex):
            dgap[i] = x_al_lr[index[i]] - x_al_lr[max(index[i] - 1, 0)]
        Ninterp = (int((x_al_lr[-1] - x_al_lr[0] - sum(dgap))
                   / float(p.delta_al)) + 1)
        x_al = numpy.zeros((Ninterp))
        stime = numpy.zeros((Ninterp))
        lon = numpy.zeros((Ninterp))
        lat = numpy.zeros((Ninterp))
        imin = 0
        imax = 0
        for i in range(0, nindex - 1):
            imax = imin + int((x_al_lr[index[i+1]-1] - x_al_lr[index[i]])
                              / float(p.delta_al)) + 1
            if imax <= (imin + 1):
                x_al[imin] = x_al_lr[index[i]]
                stime[imin] = stime_lr[index[i]]
                lon[imin] = lon_lr[index[i]]
                lat[imin] = lat_lr[index[i]]
            else:
                slicei = slice(index[i], index[i + 1])
                x_al[imin: imax] = numpy.arange(x_al_lr[index[i]],
                                                x_al_lr[index[i+1] - 1],
                                                p.delta_al)
                stime[imin: imax] = numpy.interp(x_al[imin: imax],
                                                 x_al_lr[slicei],
                                                 stime_lr[slicei])
                loncirc = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(
                                        lon_lr[slicei])))
                # if numpy.min(lon_lr[index[i]:index[i+1]])<1.
                # and numpy.max(lon_lr[index[i]:index[i+1]])>359.:
                # lontmp=lon_lr[index[i]:index[i+1]]
                # lontmp[numpy.where(lontmp>180.)]=lontmp[numpy.where(
                # lontmp>180.)]-360.
                # lon[imin:imax]=numpy.interp(x_al[imin:imax],
                # x_al_lr[index[i]:index[i+1]], lontmp)
                #    lon[imin:imax]=(lon[imin:imax]+360)%360
                # else:
                #    lon[imin:imax]=numpy.interp(x_al[imin:imax],
                # x_al_lr[index[i]:index[i+1]], lon_lr[index[i]:index[i+1]])
                lon[imin: imax] = numpy.interp(x_al[imin: imax],
                                               x_al_lr[slicei],
                                               loncirc)
                lat[imin: imax] = numpy.interp(x_al[imin: imax],
                                               x_al_lr[slicei],
                                               lat_lr[slicei])
            imin = imax
        x_al[imin:] = numpy.arange(x_al_lr[index[-1]], x_al_lr[index[-1]]
                                   + (Ninterp - imin)*p.delta_al, p.delta_al)
        stime[imin:] = numpy.interp(x_al[imin:], x_al_lr[index[-1]:],
                                    stime_lr[index[-1]:])
        loncirc = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(
                                lon_lr[index[-1]:])))
        lon[imin:] = numpy.interp(x_al[imin:], x_al_lr[index[-1]:], loncirc)
        lat[imin:] = numpy.interp(x_al[imin:], x_al_lr[index[-1]:],
                                  lat_lr[index[-1]:])
    else:
        Ninterp = int((x_al_lr[-2] - x_al_lr[0]) / float(p.delta_al)) + 1
        x_al = numpy.zeros((Ninterp))
        stime = numpy.zeros((Ninterp))
        lon = numpy.zeros((Ninterp))
        lat = numpy.zeros((Ninterp))
        x_al = numpy.arange(x_al_lr[0], x_al_lr[-2], p.delta_al)
        stime = numpy.interp(x_al, x_al_lr[:-1], stime_lr[:-1])
        loncirc = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(lon_lr[:-1])))
        lon = numpy.interp(x_al, x_al_lr[:-1], loncirc)
        lat = numpy.interp(x_al, x_al_lr[:-1], lat_lr[:-1])
    lon = lon % 360
    nfile = '{}.nc'.format(orbitfile[:-4])
    orb = rw_data.Sat_nadir(nfile=nfile)
    orb.x_al = x_al
    orb.time = stime
    orb.lon = lon
    orb.lat = lat
    orb.cycle = tcycle
    orb.al_cycle = distance[-1]
    orb.passtime = numpy.sort(passtime)
    orb.timeshift = p.timeshift
    return orb


def orbit2swath(modelbox, p, orb):
    '''Computes the swath of SWOT satellites on a subdomain from an orbit.
    The path of the satellite is given by the orbit file and the subdomain
    corresponds to the one in the model. Note that a subdomain can be manually
    added in the parameters file. \n
    Inputs are satellite orbit (p.filesat), subdomain (modelbox), Swath
    parameters (half gap distance p.halfgap, half swath distance p.halfswath,
    along track
    resolution p.delta_al, across track resolution p.delta_ac). \n
    Outputs are netcdf files containing SWOT grid (along track distance x_al,
    across track distance from nadir x_ac, longitude lon and latitude lat,
    number of days in a cycle cycle, distance crossed in a cycle cycle_al,
    time'''
    ''' Compute orbit from Swath '''
    # - Load altimeter orbit
    npoints = 1
    x_al = orb.x_al
    stime = orb.time
    lon = orb.lon
    lat = orb.lat
    tcycle = orb.cycle
    al_cycle = orb.al_cycle
    passtime = orb.passtime
    # - Compute accross track distances from nadir
    # Number of points in half of the swath
    nhalfswath = int((p.halfswath-p.halfgap)/p.delta_ac) + 1
    # Across track distance from nadir
    x_ac = numpy.zeros((2*nhalfswath))
    for i in range(0, int(nhalfswath)):
        x_ac[i] = -(nhalfswath - i)*p.delta_ac - p.halfgap + p.delta_ac
        x_ac[i + nhalfswath] = i * p.delta_ac + p.halfgap

    # - Computation of SWOT grid and storage by passes
    logger.info('\n Compute SWOT grid')
    # Detect first pass that is in the subdomain
    ipass0 = 0
    # strpass = []
    # Loop on all passes after the first pass detected
    for ipass in range(ipass0, numpy.shape(passtime)[0]):
        # Detect indices corresponding to the pass
        if ipass == numpy.shape(passtime)[0]-1:
            ind = numpy.where((stime >= passtime[ipass]))[0]
        else:
            ind = numpy.where((stime >= passtime[ipass])
                              & (stime < passtime[ipass+1]))[0]
        nind = numpy.shape(ind)[0]
        # Compute swath grid if pass is in the subdomain
        if nind > 5:
            mod_tools.update_progress(float(ipass+1)
                                      / float(numpy.shape(passtime)[0]),
                                      'selected pass: ' + str(ipass+1), None)
            # Initialize SWOT grid, grid variables and Satellite
            # direction and Location
            filesgrid = '{}_p{:03d}.nc'.format(p.filesgrid,ipass + 1)
            sgrid = rw_data.Sat_SWOT(nfile=filesgrid)
            sgrid.x_al = x_al[ind]
            sgrid.x_ac = x_ac
            sgrid.cycle = tcycle
            sgrid.al_cycle = al_cycle
            sgrid.time = stime[ind]
            sgrid.lon = numpy.zeros((nind, 2*nhalfswath))
            sgrid.lat = numpy.zeros((nind, 2*nhalfswath))
            SatDir = numpy.zeros((int(nind/npoints), 3))
            SatLoc = numpy.zeros((int((nind)/npoints), 3))

            # Initialize Nadir track, grid variables
            filengrid = '{}nadir_p{:03d}.nc'.format(p.filesgrid,ipass + 1)
            ngrid = rw_data.Sat_nadir(nfile=filengrid)
            ngrid.x_al = x_al[ind]
            ngrid.cycle = tcycle
            ngrid.al_cycle = al_cycle
            ngrid.time = stime[ind]

            # Project in cartesian coordinates satellite ground location
            s2cart = mod_tools.spher2cart(lon[ind[0]: ind[-1]+1: npoints],
                                          lat[ind[0]: ind[-1]+1: npoints])
            SatLoc[:, 0], SatLoc[:, 1], SatLoc[:, 2] = s2cart
            # Compute satellite direction (SatLoc is periodic)
            SatDir[1: -1, 0] = ((SatLoc[2:, 0]-SatLoc[: -2, 0])
                                / numpy.sqrt(SatLoc[1: -1, 0]**2
                                + SatLoc[1: -1, 1]**2 + SatLoc[1: -1, 2]**2))
            SatDir[1: -1, 1] = ((SatLoc[2:, 1] - SatLoc[: -2, 1])
                                / numpy.sqrt(SatLoc[1: -1, 0]**2
                                + SatLoc[1: -1, 1]**2 + SatLoc[1: -1, 2]**2))
            SatDir[1: -1, 2] = ((SatLoc[2:, 2] - SatLoc[: -2, 2])
                                / numpy.sqrt(SatLoc[1: -1, 0]**2
                                + SatLoc[1: -1, 1]**2 + SatLoc[1: -1, 2]**2))
            SatDir[-1] = SatDir[-2]
            SatDir[0] = SatDir[1]
            # Rotate from earth center around satellite direction to compute
            # swath points of angles between the borders of the swath in left
            # and right swath
            for i in range(0, nind, npoints):
                for j in range(0, int(nhalfswath)):
                    R = mod_tools.rotationmat3D(float((j*p.delta_ac+p.halfgap)
                                                / (const.Rearth*10**-3)),
                                                SatDir[int(i/npoints), :])
                    ObsLoc = numpy.dot(R, SatLoc[int(i/npoints)])
                    cs = mod_tools.cart2spher(ObsLoc[0], ObsLoc[1], ObsLoc[2])
                    sgrid.lon[i, nhalfswath+j], sgrid.lat[i, nhalfswath+j] = cs
                    ObsLoc = numpy.dot(numpy.transpose(R),
                                       SatLoc[int(i/npoints)])
                    sgrid.lon[i, nhalfswath-j-1], sgrid.lat[i, nhalfswath-j-1] = cs
                    if npoints > p.delta_al:
                        if i >= npoints:
                            sgrid.lon[i-npoints: i, nhalfswath+j] = numpy.arange(sgrid.lon[i-npoints, nhalfswath+j], sgrid.lon[i, nhalfswath+j], (sgrid.lon[i, nhalfswath+j]-sgrid.lon[i-npoints, nhalfswath+j])/npoints)
                            sgrid.lat[i-npoints: i, nhalfswath+j] = numpy.arange(sgrid.lat[i-npoints, nhalfswath+j], sgrid.lat[i, nhalfswath+j], (sgrid.lat[i, nhalfswath+j]-sgrid.lat[i-npoints, nhalfswath+j])/npoints)
            # if npoints>p.delta_al:
            # print 'interp not coded'
            # for j in range(0, 2*int(nhalfswath+1)):
            # sgrid.lon[:,j]=numpy.arange(sgrid.lon[0,j], sgrid.lon[ind[-1],j],
            # (sgrid.lon[-1,j]-sgrid.lon[0,j])/npoints)
            # sgrid.lat[:,j]=numpy.arange(sgrid.lat[0,j], sgrid.lat[-1,j],
            # (sgrid.lat[-1,j]-sgrid.lat[0,j])/npoints)
            # Save Sgrid object
            sgrid.timeshift = orb.timeshift
            ngrid.timeshift = orb.timeshift
            ngrid.lon = (lon[ind] + 360) % 360
            ngrid.lat = lat[ind]
            sgrid.lon_nadir = (lon[ind] + 360) % 360
            sgrid.lat_nadir = lat[ind]
            # Remove grid file if it exists and save it
            if os.path.exists(filesgrid):
                os.remove(filesgrid)
            sgrid.write_swath()
            if p.nadir:
                if os.path.exists(filengrid):
                    os.remove(filengrid)
                ngrid.write_orb()
    mod_tools.update_progress(1,  'All swaths have been processed', ' ')
    return None
