import numpy
from scipy import signal
from math import pi
from scipy.fftpack import fft

def psd1d(hh = None,dx = 1.,tap = 0.05, detrend = True):


  hh = hh-numpy.mean(hh)
  nx = numpy.shape(hh)[0]

  if detrend:
    hh = signal.detrend(hh)

  if tap>0:  
    ntaper = numpy.int(tap * nx + 0.5)
    taper = numpy.zeros(nx)+1.
    taper[:ntaper] = numpy.cos(numpy.arange(ntaper)/(ntaper-1.)*pi/2+3*pi/2)
    taper[-ntaper:] = numpy.cos(-numpy.arange(-ntaper+1,1)/(ntaper-1.)*pi/2+3*pi/2)
    hh = hh*taper

  ss = fft(hh)
  if nx % 2 != 0 : nx-=1

  ff = numpy.arange(1,nx/2-1)/(nx*dx)

  PSD = 2*dx/(nx)*numpy.abs(ss[1:int(nx/2)-1])**2


  return ff, PSD
