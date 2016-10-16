''' Plot map of model SSH, model SSH colocated on one SWOT pass, SWOT-like SSH on this pass'''
import matplotlib.pylab as plt
import numpy
isBasemap=True
try: from mpl_toolkits.basemap import Basemap
except: isBasemap=False
import matplotlib as mpl
import swotsimulator.rw_data as rw_data
from netCDF4 import Dataset
import params as p

modelbox = [-130., -123., 42., 48.]
filemodel = p.indatadir + 'OREGON_0101.nc'
filegrid = p.indatadir + 'OREGON_grd.nc'
fileswot = p.outdatadir + 'OREGON_swot292_c01_p024.nc'
proj='cyl'
vmin = -0.15
vmax = 0.15
#fig = plt.figure(figsize=(35,10))
fig, axes = plt.subplots(figsize=(36, 7.5), nrows=1, ncols=3, sharex=False,
#fig, axes = plt.subplots(nrows=1, ncols=3, sharex=False,
                         sharey=True)
#plt.clf()
size = fig.get_size_inches()*fig.dpi # get fig size in pixels
indice = 0
tloc=0.09
tfont=24
for indice in range(numpy.shape(axes.flat)[0]):
    ax = axes[indice]
    print(ax, axes)
    if isBasemap is True:
        m = Basemap(llcrnrlon=modelbox[0],
                    llcrnrlat=modelbox[2],
                    urcrnrlon=modelbox[1],
                    urcrnrlat=modelbox[3],
                    resolution='h',
                    projection=proj,
                    lon_0=(modelbox[1]-modelbox[0])/2,
                    lat_0=(modelbox[3]-modelbox[2])/2,
                    ax=ax)
        m.drawcoastlines()
        m.fillcontinents(color='#FFE4B5', lake_color = 'aqua')
        m.drawmeridians(numpy.arange(int(modelbox[0]), int(modelbox[1]),
                        int((modelbox[1]-modelbox[0])/4.)),
                        labels = [0, 0, 0, -1])
        if indice == 0:
            m.drawparallels(numpy.arange(int(modelbox[2]), int(modelbox[3]),
                            int((modelbox[3]-modelbox[2])/4.)),
                            labels = [-1, 0, 0, 0])
        else:
            m.drawparallels(numpy.arange(int(modelbox[2]), int(modelbox[3]),
                            int((modelbox[3]-modelbox[2])/4.)),
                            labels = [0, 0, 0, 0])
    if indice == 0:
        datamodel = rw_data.NETCDF_MODEL(filemodel, var=p.var, lon=p.lon, lat=p.lat)
        datamodel.read_var()
        datagrid = rw_data.NETCDF_MODEL(filegrid, lon=p.lon, lat=p.lat)
        datagrid.read_coordinates()
        datagrid.vlon[datagrid.vlon > 180] = datagrid.vlon[datagrid.vlon > 180] - 360
        if isBasemap is True:
            x, y = m(datagrid.vlon[:, :], datagrid.vlat[:, :])
        else:
            x = datagrid.vlon
            y = datagrid.vlat
        SSH = datamodel.vvar
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        SSH[SSH == 0] = numpy.nan
        SSH = numpy.ma.array(SSH, mask=numpy.isnan(SSH))
        col = ax.pcolormesh(x[:, :], y[:, :], SSH[:, :], norm=norm)
        ax.set_title('(a)', y=-tloc, fontsize=tfont) #size[1])
        #ax.clim(vmin, vmax)
    else:
        data = rw_data.Sat_SWOT(file=fileswot)
        if indice == 1:
            data.load_swath(SSH_model=[])
            SSH = data.SSH_model
            stitle = '(b)'
        elif indice == 2:
            data.load_swath(SSH_obs=[])
            SSH = data.SSH_obs
            stitle = '(c)'
        nac = numpy.shape(data.lon)[1]
        data.lon[data.lon > 180] = data.lon[data.lon > 180] - 360
        if isBasemap is True:
            x, y = m(data.lon[:, :], data.lat[:, :])
        else:
            x = data.lon
            y = data.lat
        SSH[abs(SSH) > 1000] = numpy.nan
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        SSH[SSH == 0] = numpy.nan
        SSH = numpy.ma.array(SSH, mask=numpy.isnan(SSH))
        col = ax.pcolormesh(x[:, :int(nac/2)], y[:, :int(nac/2)],
                            SSH[:, :int(nac/2)], norm=norm)
        #ax.clim(vmin, vmax)
        col = ax.pcolormesh(x[:, int(nac/2):], y[:, int(nac/2):],
                            SSH[:, int(nac/2):], norm=norm)
        ax.set_title(stitle, y=-tloc, fontsize=tfont) #size[1])
        #ax.clim(vmin, vmax)
    #indice += 1
    print(indice)
cax, kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
plt.colorbar(col, cax=cax, **kw)
plt.savefig('Fig3.png')
