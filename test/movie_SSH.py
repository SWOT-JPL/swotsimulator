from matplotlib import pyplot
isBasemap=True
try:
    import cartopy
except ImportError:
    isBasemap=False
import matplotlib
import swotsimulator.rw_data as rw_data
import params as p
import glob
import numpy


proj = 'cyl'
listfile = glob.glob(p.file_output+'_c01_p*.nc')
if p.modelbox is not None:
  modelbox = p.modelbox
elif p.config == "OREGON":
    modelbox = [-130., -123., 42., 48.]
else:
    modelbox = [float(x) for x in input("specify your modelbox with format lower_lon, upper_lon, lower_lat, upper_lat: ").split(',')]

image = 0
for coordfile in listfile:
    fig = pyplot.figure(figsize=(12,12))
    ax1 = pyplot.subplot(2, 1, 1)
    ax2 = pyplot.subplot(2, 2, 1)
    pyplot.clf()
    #pyplot.ion()
    if isBasemap is True:
        projection = cartopy.crs.PlateCarree()
        transform = cartopy.crs.PlateCarree()
        ax1 = pyplot.subplot(2, 1, 1, projection=projection)
        ax2 = pyplot.subplot(2, 2, 1, projection=projection)
        if modelbox is None:
            projection = cartopy.crs.Orthographic(0, 0)
        ax1 = pyplot.axes(projection=projection)
        ax2 = pyplot.axes(projection=projection)
        if modelbox is not None:
            ax2.set_extent([modelbox[0], modelbox[1],  modelbox[2], modelbox[3]],
                           crs=transform)
            norder = 6
        else:
            ax1.set_global()
            ax2.set_global()
            norder = 1
        #ax.add_feature(cartopy.feature.OCEAN, zorder=norder)
        #ax.add_feature(cartopy.feature.LAND, zorder=norder, edgecolor='black')
        gl1 = ax1.gridlines(crs=transform, draw_labels=True, color='gray',
                            linestyle='--', alpha=0.5)
        gl2 = ax2.gridlines(crs=transform, draw_labels=True, color='gray',
                            linestyle='--', alpha=0.5)
        gl1.xlabels_top = False
        gl1.ylabels_left = False
        gl2.xlabels_top = False
        gl2.ylabels_left = False
    # Loop on files
    print(coordfile)
    data = rw_data.Sat_SWOT(nfile=coordfile)
    data.load_swath(ssh_model=[])
    nac = numpy.shape(data.lon)[1]
    if modelbox[0] < 0:
        data.lon[numpy.where(data.lon > 180)] = (data.lon[numpy.where(
                                            data.lon > 180)] - 360)
    norm = matplotlib.colors.Normalize(vmin=-0.1, vmax=0.1)
    SSH = data.ssh_model
    SSH2 = data.ssh_obs
    mask = ((SSH==0) | (abs(SSH) > 9999.))
    SSH = numpy.ma.array(SSH, mask=mask)
    _min = -1
    _max = 1
    if _max < _min:
        continue
    ax1.pcolormesh(data.lon[:, :int(nac/2)], data.lat[:, :int(nac/2)],
                      SSH[:, :int(nac/2)], cmap='jet',
                      transform=transform)
    pyplot.clim(_min, _max)
    ax1.pcolormesh(data.lon[:, int(nac/2):], data.lat[:, int(nac/2):],
                      SSH[:, int(nac/2):], cmap='jet',
                      transform=transform)
    pyplot.clim(_min, _max)
    ax2.pcolormesh(data.lon[:, :int(nac/2)], data.lat[:, :int(nac/2)],
                      SSH2[:, :int(nac/2)], cmap='jet',
                      transform=transform)
    pyplot.clim(_min, _max)
    ax2.pcolormesh(data.lon[:, int(nac/2):], data.lat[:, int(nac/2):],
                      SSH2[:, int(nac/2):], cmap='jet',
                      transform=transform)
    pyplot.clim(_min, _max)
    pyplot.colorbar()
    pyplot.title('SWOT like data for config {}'.format(p.config))
    pyplot.savefig('./pict/img_{:04d}'.format(image))
    image += 1
