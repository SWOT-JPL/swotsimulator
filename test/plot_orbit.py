from matplotlib import pyplot
isBasemap=True
try: from mpl_toolkits.basemap import Basemap
except: isBasemap=False
import matplotlib.cm as cm
import matplotlib as mpl
import swotsimulator.rw_data as rw_data
import params as p
import glob
import numpy


proj = 'cyl'
listfile = glob.glob(p.file_output+'_c01_p*.nc')
if p.modelbox is not None:
  modelbox = p.modelbox
elif p.config=="OREGON":
  modelbox = [-130., -123., 42., 48.]
else:
  modelbox = [float(x) for x in input("specify your modelbox with format lower_lon, upper_lon, lower_lat, upper_lat: ").split(',')]
  print(modelbox)
fig = pyplot.figure
pyplot.clf()
#pyplot.ion()

listfile = sorted(listfile)
if isBasemap is True:
    m = Basemap(llcrnrlon=modelbox[0], \
            llcrnrlat=modelbox[2], \
            urcrnrlon=modelbox[1], \
            urcrnrlat=modelbox[3], \
            resolution='h', \
            projection=proj, \
            lon_0=(modelbox[1]-modelbox[0])/2, \
            lat_0=(modelbox[3]-modelbox[2])/2)
    m.drawcoastlines()
    m.fillcontinents(color='#FFE4B5', lake_color='aqua')
    m.drawmeridians(numpy.arange(int(modelbox[0]),
                    int(modelbox[1]) + 0.1, (modelbox[1] - modelbox[0])/7.),
                    labels = [0, 0, 0, 2])
    m.drawparallels(numpy.arange(int(modelbox[2]),
                    int(modelbox[3]) + 0.1, (modelbox[3] - modelbox[2])/7.),
                    labels = [2, 0, 0, 0])
for coordfile in listfile:
    print(coordfile)
    data = rw_data.Sat_SWOT(nfile=coordfile)
    data.load_swath()
    if modelbox[0] < 0:
        data.lon_nadir[numpy.where(data.lon_nadir > 180)] = (data.lon_nadir[
                                          numpy.where(data.lon_nadir>180)]-360)
    if isBasemap is True:
        x,y = m(data.lon_nadir,data.lat_nadir)
    else:
        x = data.lon_nadir
        y = data.lat_nadir
    ipass = coordfile[-6:-3]
    if ((numpy.max(data.lon_nadir) > modelbox[0])
        and (numpy.min(data.lon_nadir) < modelbox[1])):
        C = pyplot.plot(x, y, '.', label='{}'.format(ipass))
    else:
        C = pyplot.plot(x, y, '.')
pyplot.legend(loc=3, mode="expand", ncol=3)
pyplot.savefig('{}_sat_pass.png'.format(p.config))
pyplot.show()
