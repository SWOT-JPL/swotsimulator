''' Spectral  and algebra tools  for swot simulator. \n
Contains the following functions:
- gen_signal1d: compute 1d random error from spectrum \n
- gen_signal2d: compute 2d random error from spectrum \n
- rotationmat3d: rotate data  \n
- spher2cart: convert spherical to cartesian coordinates \n
- cart2spher: convert cartesian to spherical coordinates \n
- update_progress: Progress bar'''
import numpy
import math
import sys


def gen_coeff_signal1d(f, PS, nc):
    '''Generate nc random coefficient from a spectrum PS
    with frequencies f. \n
    Return Amplitude, phase and frequency of nc realisations'''

    f1 = min(f)
    f2 = max(f)
    logf = numpy.log10(f)
    logf1 = numpy.log10(f1)
    logf2 = numpy.log10(f2)
    # ''' Compute nc random vectors in [logf1, logf2] '''
    logfr = (logf2 - logf1)*numpy.random.random(nc) + logf1
    fr = 10**(logfr)
    # ''' Compute amplitude for random vectors '''
    logPS = numpy.log10(PS)
    logPSr = numpy.interp(logfr, logf, logPS)
    PSr = 10**(logPSr)
    A = numpy.sqrt(0.5 * PSr * fr * ((f2/f1)**(1./nc)-1))
    # ''' Compute nc phases in (0,2*pi)'''
    phi = 2 * math.pi * numpy.random.random(nc)
    return A, phi, fr


def gen_coeff_signal2d(f, PS, nc):
    '''Generate nc random coefficient from a spectrum PS
    with frequencies f. \n
    Inputs are: frequency [f], spectrum [PS], number of realisation [nc]
    Return Amplitude, phase and frequency in 2D (frx, fry) of nc realisations'''

    # ''' Compute nc random vectors in an annular
    # (radius is between (min(f), max(f)) '''
    f1 = min(f)
    f2 = max(f)
    logf = numpy.log10(f)
    fr = (f2 - f1)*numpy.random.random(nc) + f1
    logfr = numpy.log10(fr)
    dir = 2. * math.pi*numpy.random.random(nc)
    frx = fr * numpy.cos(dir)
    fry = fr * numpy.sin(dir)

    # ''' Compute coefficients corresponding to random vectors '''
    logPS = numpy.log10(PS)
    logPSr = numpy.interp(numpy.log10(fr), logf, logPS)
    PSr = 10.**logPSr
    A = numpy.sqrt(0.5*PSr*2*math.pi*(f2 - f1)/(nc))

    # ''' Compute nc phases in (0,2*pi)'''
    phi = 2*math.pi*numpy.random.random(nc)

    return A, phi, frx, fry


def rotationmat3D(theta, axis):
    ''' Creates a rotation matrix: Slow method. \n
    Inputs are rotation angle theta and rotation axis axis.
    The rotation matrix correspond to a rotation of angle theta
    with respect to axis axis. \n
    Return the rotation matrix.'''
    mat = numpy.eye(3, 3)
    axis = axis / math.sqrt(numpy.dot(axis, axis))
    a = math.cos(theta/2.)
    b, c, d = -axis*math.sin(theta/2.)

    return numpy.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                       [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
                       [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]])


def sign(a):
    return (a > 0) - (a < 0)


def spher2cart(lon, lat):
    ''' Convert spherical coordinates to cartesian coordinates.\n
    Inputs are longitude, latitude. \n
    Return x, y, z'''
    x = numpy.cos(lon*math.pi/180.) * numpy.cos(lat*math.pi/180.)
    y = numpy.sin(lon*math.pi/180.) * numpy.cos(lat*math.pi/180.)
    z = numpy.sin(lat*math.pi/180.)
    return x, y, z


def cart2sphervect(x, y, z):
    ''' Convert cartesian coordinates to spherical coordinates. \n
    Inputs are cartiesian coordinates x, y, z. \n
    Return lon, lat. '''
    norm = numpy.sqrt(x*x + y*y + z*z)
    lat = numpy.arcsin(z/norm) * 180./math.pi
    lon = numpy.arctan(y/x) % (2*math.pi)
    if (x < 0).any():
        lon[x < 0] = (numpy.arctan(y [x < 0]/x[x < 0]) % (2*math.pi)
                      + x[x < 0] / x[x < 0]*math.pi)
    lon = lon * 180/math.pi
    return lon % 360, lat

def cart2spher(x, y, z):
    ''' Convert cartesian coordinates to spherical coordinates. \n
    Inputs are cartiesian coordinates x, y, z. \n
    Return lon, lat. '''
    norm = numpy.sqrt(x*x + y*y + z*z)
    lat = numpy.arcsin(z/norm) * 180./math.pi
    if (x < 0):
        lon = (numpy.arctan(y/x) % (2*math.pi)
                      + max((x < 0) - (x > 0), 0)*math.pi)
        # + max(-sign(x),0)*math.pi
    else:
        lon = numpy.arctan(y/x) % (2*math.pi)
    # lon=(lon+2*math.pi)%2*math.pi
    lon = lon * 180/math.pi
    return lon % 360, lat


def update_progress(progress, arg1, arg2):
    '''Creation of a progress bar: print on screen the progress of the run'''
    barLength = 30  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    if arg1 and arg2:
        text = "\rPercent: [{0}] {1}%, {2}, {3}".format("#"*block + "-"*(barLength-block), "%.2f" % (progress*100), arg1 + ', ' + arg2, status)
    elif arg1:
        text = "\rPercent: [{0}] {1}%, {2}, {3}".format("#"*block + "-"*(barLength-block), "%.2f" % (progress*100), arg1, status)
    else:
        text = "\rPercent: [{0}] {1}%, {2} ".format("#"*block + "-"*(barLength-block), "%.2f" % (progress*100), status)
    sys.stdout.write(text)
    sys.stdout.flush()
    return progress
