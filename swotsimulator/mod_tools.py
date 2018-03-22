''' Spectral  and algebra tools  for swot simulator. \n
Contains the following functions:
- load_python_file: load and parse parameter file \n
- gen_signal1d: compute 1d random error from spectrum \n
- gen_signal2d: compute 2d random error from spectrum \n
- rotationmat3d: rotate data  \n
- spher2cart: convert spherical to cartesian coordinates \n
- cart2spher: convert cartesian to spherical coordinates \n
- update_progress: Progress bar'''
import numpy
import math
import logging
import sys
import os

# Define logger level for debug purposes
logger = logging.getLogger(__name__)


def load_python_file(file_path):
    """Load a file and parse it as a Python module."""
    if not os.path.exists(file_path):
        raise IOError('File not found: {}'.format(file_path))

    full_path = os.path.abspath(file_path)
    python_filename = os.path.basename(full_path)
    module_name, _ = os.path.splitext(python_filename)
    module_dir = os.path.dirname(full_path)
    if module_dir not in sys.path:
        sys.path.append(module_dir)

    module = __import__(module_name, globals(), locals(), [], 0)
    return module


def initialize_parameters(p):
    p.shift_lon = getattr(p, 'shift_lon', None)
    p.shift_time = getattr(p, 'p.shift_time', None)
    p.timeshift = getattr(p, 'p.timeshift', 0)
    if p.shift_time is None:
        p.timeshift = 0
    model = getattr(p, 'model', 'NETCDF_MODEL')
    p.model = model
    p.model_nan = getattr(p, 'model_nan', 0)
    p.SSH_factor = getattr(p, 'SSH_factor', 1.)
    p.nadir = getattr(p, 'nadir', True)
    p.grid = getattr(p, 'grid', 'regular')
    p.order_orbit_col = getattr(p, 'order_orbit_col', None)
    p.ice_mask = getattr(p, 'ice_mask', True)
    p.save_variables = getattr(p, 'save_variables', 'classic')
    p.savesignal = getattr(p, 'savesignal', False)
    return None


def check_path(p):
    if os.path.isdir(p.dir_setup) is False:
        logger.error('Data directory {} not found'.format(p.dir_setup))
        sys.exit(1)
    if os.path.isdir(p.indatadir) is False:
        logger.error('Input directory {} not found'.format(p.indatadir))
        sys.exit(1)
    if os.path.isdir(p.outdatadir) is False:
        logger.warn('Output directory {} did not exist and was '
                    'created'.format(p.dir_setup))
        os.makedirs(p.outdatadir)
    filesat_path = os.path.join(p.dir_setup, p.filesat)
    if os.path.isfile(filesat_path) is False:
        logger.error('Orbit file {} not found'.format(filesat_path))
        sys.exit(1)
    if p.file_input is not None:
        if os.path.isfile(p.file_input) is False:
            logger.error('Model file list {} not found'.format(p.file_input))
            sys.exit(1)
    return None


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


def gen_rcoeff_signal1d(f, PS, lambda_min, lambda_max, npseudoper, repeat):
    '''Generate nc random coefficient from a spectrum PS
    with frequencies f. \n
    Return Amplitude, phase and frequency of nc realisations'''

    logffl = numpy.arange(numpy.log10(1. / lambda_max),
                          numpy.log10(1. / lambda_min + 1. / lambda_max),
                          numpy.log10(1 + 1. / npseudoper))
    ffl = 10**(logffl)
    logf = numpy.log10(f)
    logPS = numpy.log10(PS)
    logPSl = numpy.interp(logffl, logf, logPS)
    PSl = 10**(logPSl)
    A = numpy.sqrt(2 * PSl * (ffl / npseudoper))
    phi = []
    for k in range(len(ffl)):
        phi.append(2 * math.pi * numpy.random.random(int(2 * repeat * ffl[k]
                                                     / npseudoper + 3)))
    return A, phi


def gen_signal1d(xx, A, phi, lambda_min, lambda_max, npseudoper):
    '''Generate 1d random noise signal from coefficent computed using
     gen_rcoeff_signal1d. \n
    Return The random noise signal'''
    S = numpy.full(numpy.shape(xx), 0.)
    logffl = numpy.arange(numpy.log10(1. / lambda_max),
                          numpy.log10(1. / lambda_min + 1. / lambda_max),
                          numpy.log10(1 + 1. / npseudoper))
    ffl = 10**(logffl)
    for k in range(len(ffl)):
        ka = 2 * (xx * ffl[k] / npseudoper).astype(int) + 1
        Cka = numpy.abs(numpy.sin(2 * math.pi * xx * ffl[k] / npseudoper / 2.))
        kb = (2 * ((xx + npseudoper / 2. / ffl[k]) * ffl[k]
              / npseudoper).astype(int))
        Ckb = numpy.abs(numpy.sin(2 * math.pi * (xx + npseudoper / 2 / ffl[k])
                        * ffl[k] / npseudoper / 2))
        S = (S + A[k] * numpy.cos(2 * math.pi * ffl[k] * xx + phi[k][ka]) * Cka
             + A[k] * numpy.cos(2 * math.pi * ffl[k] * xx + phi[k][kb]) * Ckb)
    return S


def gen_coeff_signal2d(f, PS, nc):
    '''Generate nc random coefficient from a spectrum PS
    with frequencies f. \n
    Inputs are: frequency [f], spectrum [PS], number of realisation [nc]
    Return Amplitude, phase and frequency in 2D (frx, fry) of nc
    realisations'''

    # ''' Compute nc random vectors in an annular
    # (radius is between (min(f), max(f)) '''
    f1 = min(f)
    f2 = max(f)
    logf = numpy.log10(f)
    fr = (f2 - f1)*numpy.random.random(nc) + f1
    # logfr = numpy.log10(fr)
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
    # mat = numpy.eye(3, 3)
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
        lon[x < 0] = (numpy.arctan(y[x < 0]/x[x < 0]) % (2*math.pi)
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
               + max(-numpy.sign(x), 0)*math.pi)
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
