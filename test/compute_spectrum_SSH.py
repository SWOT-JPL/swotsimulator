import glob
import numpy
import matplotlib.pylab as plt
import myspectools
import params as p
import swotsimulator.rw_data as rw_data


dirname = p.outdatadir
rootname = p.file_output
dx = p.delta_al #km
listfile = glob.glob(p.file_output+'_c*'+'_p*.nc')

nr = 0
distance_min = 400
f0 = numpy.linspace(1/distance_min*dx, 1/2-1/distance_min*dx,num = int(distance_min/dx/2.))
for coordfile in listfile:
    print(coordfile)
    data=rw_data.Sat_SWOT(file = coordfile)
    data.load_swath(SSH_obs = [], SSH_model = [])
    nal, nac = numpy.shape(data.SSH_obs)
    for iac in range(nac):
      if nal*dx>=distance_min:

        tap = 0.04
        ff, PSD_obs = myspectools.psd1d(hh=data.SSH_obs[:,iac],dx=dx, detrend=True, tap=tap)
        ff, PSD_model = myspectools.psd1d(hh=data.SSH_model[:,iac],dx=dx, detrend=True, tap=tap)

        try: 

          SS_obs = SS_obs + numpy.interp(f0, ff, PSD_obs)
          SS_model = SS_model + numpy.interp(f0, ff, PSD_model)
        except:
      
          SS_obs = numpy.interp(f0, ff, PSD_obs)
          SS_model = numpy.interp(f0,ff, PSD_model)
        nr+=1
  
SS_obs/=nr
SS_model/=nr

ff = f0
dff = ff[1]-ff[0]

plt.close()
plt.loglog(ff,SS_obs, color='k', label='SSH_obs'); plt.grid()
plt.loglog(ff,SS_model, color='red',lw=2, label='SSH_model')
#plt.axis([5e-3,0.25,1e-4,1e3])
plt.xlabel(u'cy/km')
plt.ylabel(u'm²/(cy/km)')
plt.legend()
plt.title('SSH spectra for SWOT-like data and model data interpolated on the swath')
plt.savefig('SSH_spectra.png')
plt.show()

