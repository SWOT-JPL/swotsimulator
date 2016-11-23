''' Wet trop error corrected with a 1 beam and a 2 beam radiometer as a function of across track and along track variables '''
import matplotlib.pylab as plt
import numpy
import matplotlib as mpl
import swotsimulator.rw_data as rw_data
import numpy
import params as p

filegrid = p.indatadir + 'OREGON_grd.nc'
fileswot = p.outdatadir + 'OREGON_swot292_c01_p024.nc'
vmin = -0.01
vmax = 0.01
fig, axes = plt.subplots(figsize=(30, 20), nrows=2, ncols=1, sharex=True,
                         sharey=False)
tloc=0.11
tfont=24
data = rw_data.Sat_SWOT(file=fileswot)
data.load_swath(pd_err_2b=[], pd=[], x_ac=[], x_al=[])
x_al, x_ac = numpy.meshgrid(data.x_al, data.x_ac)
x_al = x_al - numpy.min(numpy.min(x_al))
ax = axes[0]
SSH = data.pd_err_2b
stitle = '(a) 1 beam residual error along swath'
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
ax.set_ylabel('across track (km)')
ax.set_title(stitle, y=-tloc, fontsize=tfont) #size[1])
ax.axis('tight')
#plt.colorbar(col)
ax = axes[1]
SSH = data.pd_err_2b
stitle = '(a) 2 beams residual error along swath'
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
ax.set_ylabel('across track (km)')
ax.set_xlabel('along track (km)')
ax.set_title(stitle, y=-tloc, fontsize=tfont) #size[1])
ax.axis('tight')
cax, kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
#plt.colorbar(col)
plt.colorbar(col, cax=cax, **kw)
plt.savefig('Fig17.png')
