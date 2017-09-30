''' Statistical wet troposphere as a function of across track and along track variables '''
import matplotlib.pylab as plt
import numpy
import matplotlib as mpl
import swotsimulator.rw_data as rw_data
import numpy
import params as p

filegrid = p.indatadir + 'OREGON_grd.nc'
fileswot = p.outdatadir + 'OREGON_swot292_c01_p024.nc'
vmin = -0.05
vmax = 0.051
fig, axes = plt.subplots(figsize=(30, 10), nrows=1, ncols=1, sharex=True,
                         sharey=False)
tloc=0.11
tfont=24
data = rw_data.Sat_SWOT(file=fileswot)
data.load_swath(pd=[], x_ac=[], x_al=[])
x_al, x_ac = numpy.meshgrid(data.x_al, data.x_ac)
x_al = x_al - numpy.min(numpy.min(x_al))
ax = axes
SSH = data.pd
stitle = 'Wet troposphere'
nac = numpy.shape(data.lon)[1]
SSH[abs(SSH) > 1000] = numpy.nan
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
SSH[SSH == 0] = numpy.nan
SSH = numpy.ma.array(SSH, mask=numpy.isnan(SSH))
print(numpy.shape(x_al), numpy.shape(x_ac), numpy.shape(SSH))
col = ax.pcolormesh(x_al[(nac/2), :], x_ac[:int(nac/2), 0],
        numpy.transpose(SSH[:,:int(nac/2)]), norm=norm)
col = ax.pcolormesh(x_al[(nac/2), :], x_ac[int(nac/2):, 0],
        numpy.transpose(SSH[:,int(nac/2):]), norm=norm)
#axis('tight')
ax.set_xlabel('along track (km)')
ax.set_ylabel('across track (km)')
ax.set_title(stitle, y=-tloc, fontsize=tfont) #size[1])
plt.colorbar(col)
plt.axis('tight')
plt.savefig('Fig15.png')
