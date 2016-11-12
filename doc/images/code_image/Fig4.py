''' Plotting Karin noise along track (no location on a map) '''
import matplotlib.pylab as plt
import numpy
import matplotlib as mpl
import swotsimulator.rw_data as rw_data
from netCDF4 import Dataset
import params as p

filegrid = p.indatadir + 'OREGON_grd.nc'
fileswot = p.outdatadir + 'OREGON_swot292_c01_p024.nc'
vmin = -0.01
vmax = 0.01
fig = plt.figure(figsize=(30,10))
tloc=0.11
tfont=24
data = rw_data.Sat_SWOT(file=fileswot)
data.load_swath(karin_err=[], x_ac=[], x_al=[])
x_al, x_ac = numpy.meshgrid(data.x_al, data.x_ac)
x_al = x_al - numpy.min(numpy.min(x_al))
SSH = data.karin_err
stitle = 'Karin noise along track'
nac = numpy.shape(data.lon)[1]
SSH[abs(SSH) > 1000] = numpy.nan
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
SSH[SSH == 0] = numpy.nan
SSH = numpy.ma.array(SSH, mask=numpy.isnan(SSH))
print(numpy.shape(x_al), numpy.shape(x_ac), numpy.shape(SSH))
col = plt.pcolormesh(x_al[(nac/2), :], x_ac[:int(nac/2), 0],
        numpy.transpose(SSH[:,:int(nac/2)]), norm=norm)
plt.xlabel('along track (km)')
plt.ylabel('across track (km)')
plt.colorbar()
plt.title(stitle, y=-tloc, fontsize=tfont) #size[1])
plt.savefig('Fig4.png')
